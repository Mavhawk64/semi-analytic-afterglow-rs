use crate::modules::constants::{CCGS, O_O_RME, PII, QCGS, RM_ELEC, XH};
use crate::modules::globals::Globals;

pub struct PhotonSscInput<'a> {
    pub num_dp: usize,
    pub elec_gam_array: &'a [f64],
    pub dn_pf: &'a [f64],
    pub phot_en_array: &'a [f64],
    pub o_o_alpha_array: &'a [f64],
    pub syn_jnu: &'a [f64],
    pub i_params: &'a [i32; 7],
    pub r_params: &'a [f64; 7],
    pub globals: &'a Globals,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PhotonSscError {
    ElectronGammaArrayTooShort { len: usize, required_min_len: usize },
    ElectronDistributionTooShort { len: usize, required_min_len: usize },
    PhotonEnergyArrayTooShort { len: usize, required_min_len: usize },
    OneOverAlphaArrayTooShort { len: usize, required_min_len: usize },
    SynJnuArrayTooShort { len: usize, required_min_len: usize },
}

/// Computes synchrotron self-Compton volume emissivity.
///
/// # Arguments
///
/// * `input` - struct containing the input parameters for the computation
///   * `num_dp` - number of momentum bins in the distribution of particles
///   * `i_params` - array (length 7) of up to 7 integer parameters to be used
///     in evaluating the function
///   * `r_params` - array (length 7) of up to 7 floating point parameters to
///     be used in evaluating the function
///   * `dn_pf` - particle distribution. This is the number of particles in each
///     bin, NOT number of particles per energy (`dn/dE`). It is furthermore a
///     density `n`, not a pure number count `N` as used in other versions of
///     this subroutine
///   * `elec_gam_array` - Lorentz factor boundary values of the `dn(p)`
///     distribution
///   * `phot_en_array` - array of energy values for photons
///   * `o_o_alpha_array` - array of `1 / alpha`, where alpha is photon energy
///     in units of electron rest mass
///   * `syn_jnu` - input synchrotron volume emissivity
///
/// # Returns
///
/// * `ssc_jnu` - output synchrotron self-Compton volume emissivity
pub fn photon_ssc(input: PhotonSscInput<'_>) -> Result<Vec<f64>, PhotonSscError> {
    let PhotonSscInput {
        num_dp,
        elec_gam_array,
        dn_pf,
        phot_en_array,
        o_o_alpha_array,
        syn_jnu,
        i_params,
        r_params,
        globals,
    } = input;

    let num_ph = globals.num_ph;

    // Do some rusty error checking that isn't in the Fortran code.
    if elec_gam_array.len() <= num_dp {
        return Err(PhotonSscError::ElectronGammaArrayTooShort {
            len: elec_gam_array.len(),
            required_min_len: num_dp + 1,
        });
    }

    if dn_pf.len() < num_dp {
        return Err(PhotonSscError::ElectronDistributionTooShort {
            len: dn_pf.len(),
            required_min_len: num_dp,
        });
    }

    if phot_en_array.len() < num_ph {
        return Err(PhotonSscError::PhotonEnergyArrayTooShort {
            len: phot_en_array.len(),
            required_min_len: num_ph,
        });
    }

    if o_o_alpha_array.len() < num_ph {
        return Err(PhotonSscError::OneOverAlphaArrayTooShort {
            len: o_o_alpha_array.len(),
            required_min_len: num_ph,
        });
    }

    if syn_jnu.len() < num_ph {
        return Err(PhotonSscError::SynJnuArrayTooShort {
            len: syn_jnu.len(),
            required_min_len: num_ph,
        });
    }

    let mut ssc_jnu = vec![1.0e-99; num_ph];
    let mut max_so_far: f64 = 1.0e-99;

    // Convert syn_jnu to photon density
    let ds = r_params[0];
    let dlognu = (phot_en_array[2] / phot_en_array[1]).ln();

    let mut phot_num_den = vec![1.0e-99; num_ph];

    for i in 0..num_ph {
        if syn_jnu[i] > 1.0e-55 {
            phot_num_den[i] = syn_jnu[i] * (4.0 * PII * ds * dlognu) / (XH * CCGS);
        }
    }

    // Compute photon energy in units of electron rest mass
    let phot_en_rm: Vec<f64> = phot_en_array
        .iter()
        .take(num_ph)
        .map(|phot_en| phot_en * O_O_RME)
        .collect();

    // Invert electron Lorentz factors to save divisions later
    let mut o_o_gam_array = vec![0.0; num_dp];

    for j_el in 0..num_dp {
        // Fortran:
        //   o_o_gam_array(1:num_dp) = 1.d0 / elec_gam_array(1:num_dp)
        //
        // Rust:
        //   j_el = 0 corresponds to Fortran j_el = 1, so use elec_gam_array[j_el + 1].
        o_o_gam_array[j_el] = 1.0 / elec_gam_array[j_el + 1];
    }

    let dndtda_prefac = 2.0 * PII * (QCGS.powi(2) * O_O_RME).powi(2) * CCGS;

    //-------------------------------------------------------------------------
    // Loop over outgoing photons from high to low energy, since SSC emission
    // skews toward high energy
    //-------------------------------------------------------------------------
    for i_out in (0..num_ph).rev() {
        // Make our lives clearer later on by defining this here
        let alpha_out = phot_en_rm[i_out];

        //-----------------------------------------------------------------------
        // Loop over electron energies
        //-----------------------------------------------------------------------
        for j_el in 0..num_dp {
            // Skip computation if there are no electrons at this energy
            if dn_pf[j_el] < 1.0e-55 {
                continue;
            }

            // If outgoing photon would have more energy than pre-scatter electron
            // did, skip it
            if alpha_out > elec_gam_array[j_el + 1] - 1.0 {
                continue;
            }

            // Initialize the dN_dtdalph array
            let mut dndtdalpha_temp = vec![1.0e-99; num_ph];

            // Grab the correct value of gamma
            let gam_elec = elec_gam_array[j_el + 1];

            // Save some computation in the next loop
            let q_pp_prefac =
                alpha_out / (4.0 * gam_elec.powi(2) * (1.0 - alpha_out * o_o_gam_array[j_el]));

            //---------------------------------------------------------------------
            // Loop over incoming photons. We only care about cases when photons
            // *gain* energy from the scattering. Fortunately, the arrays for
            // photon energy are identical for incoming and outgoing photons, so
            // we only need to compare indices. That leads to the following
            // limit on k_in
            //---------------------------------------------------------------------
            for k_in in 0..=i_out {
                // If there were no photons at this energy, don't bother continuing
                if phot_num_den[k_in] < 1.0e-55 {
                    continue;
                }

                // We need two things to find dN/(dt dalpha): q'' and alpha*gamma.
                // Compute those here
                let q_pp = q_pp_prefac * o_o_alpha_array[k_in];
                let alpha_gam = phot_en_rm[k_in] * gam_elec;

                // There is a physical limit on q_pp due to the energetics of
                // scattering
                if q_pp > 1.0 {
                    continue;
                }

                // Compute dN/(dt dalpha)
                dndtdalpha_temp[k_in] = 2.0 * q_pp * q_pp.ln()
                    + (1.0 + 2.0 * q_pp) * (1.0 - q_pp)
                    + 8.0 * (1.0 - q_pp) * (alpha_gam * q_pp).powi(2)
                        / (1.0 + 4.0 * alpha_gam * q_pp);
            }

            // Take care of prefactors to get the final value of dN/dtdalpha
            let mut weighted_sum = 0.0;

            for k_in in 0..=i_out {
                weighted_sum += dndtdalpha_temp[k_in] * o_o_alpha_array[k_in] * phot_num_den[k_in];
            }

            let dn_dtdalpha = dndtda_prefac * o_o_gam_array[j_el].powi(2) * weighted_sum;

            //---------------------------------------------------------------------
            // Loop over incoming photons
            //---------------------------------------------------------------------

            // Update ssc_jnu, ignoring the prefactor for now
            ssc_jnu[i_out] += dn_dtdalpha * dn_pf[j_el];
        }

        //-----------------------------------------------------------------------
        // Loop over electron energies
        //-----------------------------------------------------------------------

        // Incorporate the prefactor into ssc_jnu and put a floor on reasonable
        // photon production. If we've reached low enough energies that SSC
        // emission stops, exit the loop
        //
        // *However*, if we haven't even reached energies with emission yet
        // (because electron energies or photon energies are too low), don't
        // exit the loop
        //
        // Final prefactor here is needed to convert d2N/dtda into d2N/dtdnu
        // (the factor h/[m_e c^2]), to convert that into dP/dnu (the factor of
        // photon energy, since E*dN/dt = dP), and to convert dP/dnu into j_nu
        // (the factor 1/4pi).
        max_so_far = max_so_far.max(ssc_jnu[i_out]);

        if ssc_jnu[i_out] > 1.0e-55 {
            ssc_jnu[i_out] = ssc_jnu[i_out] * phot_en_array[i_out] * XH / (4.0 * PII * RM_ELEC);
        } else if max_so_far > 1.0e-55 {
            break;
        }
    }

    //-------------------------------------------------------------------------
    // Loop over outgoing photons
    //-------------------------------------------------------------------------

    let _ = i_params; // not sure what this is for, but it was in the fortran code, so I set it to _ to not have an error but keep the same function signature.

    Ok(ssc_jnu)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_globals(num_ph: usize) -> Globals {
        Globals {
            num_ph,
            ..Default::default()
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn base_input<'a>(
        num_dp: usize,
        elec_gam_array: &'a [f64],
        dn_pf: &'a [f64],
        phot_en_array: &'a [f64],
        o_o_alpha_array: &'a [f64],
        syn_jnu: &'a [f64],
        i_params: &'a [i32; 7],
        r_params: &'a [f64; 7],
        globals: &'a Globals,
    ) -> PhotonSscInput<'a> {
        PhotonSscInput {
            num_dp,
            elec_gam_array,
            dn_pf,
            phot_en_array,
            o_o_alpha_array,
            syn_jnu,
            i_params,
            r_params,
            globals,
        }
    }

    #[test]
    fn returns_error_when_electron_gamma_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonSscError::ElectronGammaArrayTooShort {
                len: 4,
                required_min_len: 5
            })
        ));
    }

    #[test]
    fn returns_error_when_electron_distribution_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp - 1];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonSscError::ElectronDistributionTooShort {
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
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonSscError::PhotonEnergyArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_o_o_alpha_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_alpha_array = vec![1.0; globals.num_ph - 1];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonSscError::OneOverAlphaArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_syn_jnu_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph - 1];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonSscError::SynJnuArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_spectrum_with_num_ph_entries() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![1.0, 10.0, 20.0, 30.0, 40.0];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12, 2.0e-12, 4.0e-12, 8.0e-12];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ))
        .expect("photon_ssc should return a spectrum");

        assert_eq!(result.len(), globals.num_ph);
    }

    #[test]
    fn returned_spectrum_is_finite() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![1.0, 10.0, 20.0, 30.0, 40.0];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12, 2.0e-12, 4.0e-12, 8.0e-12];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ))
        .expect("photon_ssc should return a spectrum");

        assert!(result.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn empty_electron_distribution_returns_floor_values() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![1.0, 10.0, 20.0, 30.0, 40.0];
        let dn_pf = vec![0.0; num_dp];
        let phot_en_array = vec![1.0e-12, 2.0e-12, 4.0e-12, 8.0e-12];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![1.0e-40; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ))
        .expect("photon_ssc should return a spectrum");

        assert!(
            result
                .iter()
                .all(|value| (*value - 1.0e-99).abs() < 1.0e-120)
        );
    }

    #[test]
    fn empty_synchrotron_input_returns_floor_values() {
        let num_dp = 4;
        let globals = base_globals(4);

        let elec_gam_array = vec![1.0, 10.0, 20.0, 30.0, 40.0];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12, 2.0e-12, 4.0e-12, 8.0e-12];
        let o_o_alpha_array = vec![1.0; globals.num_ph];
        let syn_jnu = vec![0.0; globals.num_ph];
        let i_params = [0; 7];
        let r_params = [1.0; 7];

        let result = photon_ssc(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_alpha_array,
            &syn_jnu,
            &i_params,
            &r_params,
            &globals,
        ))
        .expect("photon_ssc should return a spectrum");

        assert!(
            result
                .iter()
                .all(|value| (*value - 1.0e-99).abs() < 1.0e-120)
        );
    }
}
