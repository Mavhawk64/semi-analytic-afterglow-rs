use crate::modules::constants::{CCGS, PII, QCGS, RM_ELEC, THIRD, XH, XH_BAR, XME};
use crate::modules::globals::Globals;

#[derive(Debug, Clone, PartialEq)]
pub enum PhotonSynError {
    ElectronGammaArrayTooShort { len: usize, required_min_len: usize },
    ElectronDistributionTooShort { len: usize, required_min_len: usize },
    PhotonEnergyArrayTooShort { len: usize, required_min_len: usize },
    SynXArrayTooShort { len: usize, required_min_len: usize },
    SynFArrayTooShort { len: usize, required_min_len: usize },
}

const SQRT_3: f64 = 1.732_050_807_568_877_3;

/// Calculates photon production by an electron distribution due to
/// synchrotron emission.
/// Computes j_nu: luminosity per volume, per unit frequency, per solid angle,
/// dP/[dnu dOmega dVol], in units of ergs/(sec cm^3 ster Hz), which is NOT
/// converted to a flux here
///
/// # Arguments
///
/// * `num_dp` - number of momentum bins in the distribution of particles
/// * `elec_gam_array` - Lorentz factor boundary values of dn(p) distribution
/// * `dn_pf` - particle distribution. This is the number of particles in each bin,
///   NOT number of particles per energy (dn/dE). It is furthermore a density n,
///   not a pure number count N as used in other versions of this subroutine
/// * `b_sq` - square of magnetic field strength
/// * `phot_en_array` - array of energy values for photons
/// * `i_params` - array of up to 7 integer parameters, used for debugging purposes
/// * `globals` - runtime global state containing synchrotron lookup arrays
///
/// # Returns
///
/// Synchrotron luminosity per volume, per unit frequency, per solid angle
pub fn photon_syn(
    num_dp: usize,
    elec_gam_array: &[f64],
    dn_pf: &[f64],
    b_sq: f64,
    phot_en_array: &[f64],
    _i_params: &[i32; 7],
    globals: &Globals,
) -> Result<Vec<f64>, PhotonSynError> {
    if elec_gam_array.len() <= num_dp {
        return Err(PhotonSynError::ElectronGammaArrayTooShort {
            len: elec_gam_array.len(),
            required_min_len: num_dp + 1,
        });
    }

    if dn_pf.len() < num_dp {
        return Err(PhotonSynError::ElectronDistributionTooShort {
            len: dn_pf.len(),
            required_min_len: num_dp,
        });
    }

    let num_ph = globals.num_ph;

    if phot_en_array.len() < num_ph {
        return Err(PhotonSynError::PhotonEnergyArrayTooShort {
            len: phot_en_array.len(),
            required_min_len: num_ph,
        });
    }

    if globals.syn_x.len() < num_ph {
        return Err(PhotonSynError::SynXArrayTooShort {
            len: globals.syn_x.len(),
            required_min_len: num_ph,
        });
    }

    if globals.syn_f.len() < num_ph {
        return Err(PhotonSynError::SynFArrayTooShort {
            len: globals.syn_f.len(),
            required_min_len: num_ph,
        });
    }

    // Compute the magnetic field strength at this location, and the prefactor
    // for synchrotron emission (Rybicki & Lightman, eq. 6.18).
    // Units of syn_prefac are power per unit frequency *per electron* (cgs
    // units), and the factor of sin(alpha) is not present -- it's handled by
    // the "isotropy correction" in calc_syn_Fx.
    let b_loc = b_sq.sqrt();
    let syn_prefac = SQRT_3 * QCGS.powi(3) * b_loc / (2.0 * PII * RM_ELEC);

    // Initialize syn_jnu and omega_phot_array prior to the main loop; the latter
    // array is computed here to save divisions later
    let mut syn_jnu = vec![1.0e-99; num_ph];

    let omega_phot_array: Vec<f64> = phot_en_array
        .iter()
        .take(num_ph)
        .map(|phot_en| phot_en / XH_BAR)
        .collect(); // from E = h_bar*omega

    //-------------------------------------------------------------------------
    // Now loop over momenta (i.e. energy) of electrons and calculate the
    // synchrotron emission.
    //-------------------------------------------------------------------------

    for l_mom in 1..=num_dp {
        // Find the number of electrons in this momentum bin; skip calculation if
        // the bin is empty
        let den_elec = dn_pf[l_mom - 1];

        if den_elec <= 1.0e-60 {
            continue;
        }

        // Compute the Lorentz factor of electrons in this energy bin; if the
        // energy is too low (~3 MeV), assume they don't contribute any
        // synchrotron emission
        let gam_elec = (elec_gam_array[l_mom - 1] * elec_gam_array[l_mom]).sqrt();

        if gam_elec < 6.0 {
            continue;
        }

        // Equation 6.17c of Rybicki & Lightman, without the factor of sin(alpha)
        let o_o_omega_c = (2.0 * XME * CCGS) / (3.0 * gam_elec.powi(2) * QCGS * b_loc);

        //-----------------------------------------------------------------------
        // Calculate the F factor in Rybicki & Lightman 6.18, given by Equation
        // 6.31c
        //-----------------------------------------------------------------------

        let mut i_x_lo = 0usize;

        for m_ph in 1..=num_ph {
            // Angular frequency of photons at this energy, and corresponding value
            // of x for the calculation.
            // Again, a check to skip the calculation if we don't expect any
            // emission at this energy. We can skip all higher energies as well,
            // since they will have even larger values of xxx.
            let xxx = omega_phot_array[m_ph - 1] * o_o_omega_c;

            if xxx >= 30.0 {
                break;
            }

            // Find the appropriate value of syn_x for this photon energy
            for i_x in i_x_lo..=(num_ph - 2) {
                if globals.syn_x[i_x + 1] > xxx {
                    i_x_lo = i_x;
                    break;
                }
            }

            // Interpolate between values of syn_F to get x_F. If xxx < syn_x(1),
            // it falls in the power-law part of the function F(x)
            let x_f = if xxx < globals.syn_x[1] {
                (xxx * globals.o_o_sone).powf(THIRD) * globals.syn_f[1]
            } else {
                globals.syn_f[i_x_lo]
                    + (xxx - globals.syn_x[i_x_lo])
                        / (globals.syn_x[i_x_lo + 1] - globals.syn_x[i_x_lo])
                        * (globals.syn_f[i_x_lo + 1] - globals.syn_f[i_x_lo])
            };

            // The electron spectra passed to photon_syn contain the number of
            // particles per cc in each momentum bin. Therefore, the returned
            // photon spectra have units of [erg/sec-cm^3].
            // NOTE: Eq. (6.18) from Rybicki & Lightman, [xnum_elec * syn_prefac *
            // x_F], is energy production rate per frequency, dP/dw. Since
            // omega_phot is E/hbar, dw = dE/hbar, and so omega_phot/dw = E/dE.
            // Then [dP/dw * omega_phot] is equal to [dP/dE * E], or dP/d(lnE).
            // Only include emission if it's sufficiently positive.
            let tmp_add = den_elec * omega_phot_array[m_ph - 1] * syn_prefac * x_f;

            if tmp_add > 1.0e-60 {
                syn_jnu[m_ph - 1] += tmp_add;
            }

        }

        //-----------------------------------------------------------------------
        // F factor computed for all photon energies
        //-----------------------------------------------------------------------
    }

    //-------------------------------------------------------------------------
    // Loop over electrons finished
    //-------------------------------------------------------------------------

    // Per Rybicki & Lightman, Eq. 1.16, j_nu = P_nu / 4pi. But we don't have
    // P_nu exactly: we have dP/d(lnE). So divide by photon energy, multiply
    // by Planck's constant, and divide by 4pi to convert from dP/d(lnE) to
    // j_nu
    for m_ph in 0..num_ph {
        if syn_jnu[m_ph] > 1.0e-55 {
            syn_jnu[m_ph] = syn_jnu[m_ph] * XH / (phot_en_array[m_ph] * 4.0 * PII);
        }
    }

    Ok(syn_jnu)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_globals(num_ph: usize) -> Globals {
        Globals {
            num_ph,
            syn_x: vec![1.0e-4, 1.0e-2, 1.0, 10.0, 100.0],
            syn_f: vec![0.1, 0.5, 1.0, 0.2, 0.01],
            o_o_sone: 1.0,
            ..Default::default()
        }
    }

    fn base_arrays(num_dp: usize, num_ph: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let elec_gam_array: Vec<f64> = (0..=num_dp).map(|i| 10.0 + i as f64).collect();

        let dn_pf: Vec<f64> = vec![1.0; num_dp];

        let phot_en_array: Vec<f64> = (0..num_ph).map(|i| 1.0e-12 * (i as f64 + 1.0)).collect();

        (elec_gam_array, dn_pf, phot_en_array)
    }

    #[test]
    fn returns_error_when_electron_gamma_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        );

        assert!(matches!(
            result,
            Err(PhotonSynError::ElectronGammaArrayTooShort {
                len: 4,
                required_min_len: 5
            })
        ));
    }

    #[test]
    fn returns_error_when_distribution_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp - 1];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        );

        assert!(matches!(
            result,
            Err(PhotonSynError::ElectronDistributionTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_photon_energy_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph - 1];
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        );

        assert!(matches!(
            result,
            Err(PhotonSynError::PhotonEnergyArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_syn_x_array_is_too_short() {
        let num_dp = 4;

        let mut globals = base_globals(4);
        globals.syn_x = vec![1.0e-4, 1.0e-2, 1.0];

        let (elec_gam_array, dn_pf, phot_en_array) = base_arrays(num_dp, globals.num_ph);
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        );

        assert!(matches!(
            result,
            Err(PhotonSynError::SynXArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_syn_f_array_is_too_short() {
        let num_dp = 4;

        let mut globals = base_globals(4);
        globals.syn_f = vec![0.1, 0.5, 1.0];

        let (elec_gam_array, dn_pf, phot_en_array) = base_arrays(num_dp, globals.num_ph);
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        );

        assert!(matches!(
            result,
            Err(PhotonSynError::SynFArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_spectrum_with_num_ph_entries() {
        let num_dp = 4;
        let globals = base_globals(4);

        let (elec_gam_array, dn_pf, phot_en_array) = base_arrays(num_dp, globals.num_ph);
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        )
        .expect("photon_syn should return a spectrum");

        assert_eq!(result.len(), globals.num_ph);
    }

    #[test]
    fn returned_spectrum_is_finite() {
        let num_dp = 4;
        let globals = base_globals(4);

        let (elec_gam_array, dn_pf, phot_en_array) = base_arrays(num_dp, globals.num_ph);
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        )
        .expect("photon_syn should return a spectrum");

        assert!(result.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn empty_electron_distribution_returns_floor_values() {
        let num_dp = 4;
        let globals = base_globals(4);

        let (elec_gam_array, _dn_pf, phot_en_array) = base_arrays(num_dp, globals.num_ph);
        let dn_pf = vec![0.0; num_dp];
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        )
        .expect("photon_syn should return a spectrum");

        assert!(
            result
                .iter()
                .all(|value| (*value - 1.0e-99).abs() < 1.0e-120)
        );
    }

    #[test]
    fn low_energy_electrons_return_floor_values() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_syn(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            1.0,
            &phot_en_array,
            &i_params,
            &globals,
        )
        .expect("photon_syn should return a spectrum");

        assert!(
            result
                .iter()
                .all(|value| (*value - 1.0e-99).abs() < 1.0e-120)
        );
    }
}
