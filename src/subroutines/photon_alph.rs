use crate::modules::constants::{CCGS, O_O_RME, PII, QCGS, RM_ELEC, THIRD, XH, XME};
use crate::modules::globals::Globals;
use crate::subroutines::calc_ssa_x_f::CalcSsaXFError;
use crate::subroutines::calc_ssa_x_f::calc_ssa_x_f;

const SQRT_3: f64 = 1.732_050_807_568_877_3;

pub struct PhotonAlphInput<'a> {
    pub num_dp: usize,
    pub elec_gam_array: &'a [f64],
    pub dn_pf: &'a [f64],
    pub b_sq: f64,
    pub phot_en_array: &'a [f64],
    pub o_o_nu_sq: &'a [f64],
    pub i_params: &'a [i32; 7],
    pub globals: &'a Globals,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PhotonAlphError {
    ElectronGammaArrayTooShort { len: usize, required_min_len: usize },
    ElectronDistributionTooShort { len: usize, required_min_len: usize },
    PhotonEnergyArrayTooShort { len: usize, required_min_len: usize },
    OneOverNuSqArrayTooShort { len: usize, required_min_len: usize },
    SynXArrayTooShort { len: usize, required_min_len: usize },
    SynFArrayTooShort { len: usize, required_min_len: usize },
    NegativeAlpha { m_ph: usize, alpha_nu: f64 },
    CalcSsaXF(CalcSsaXFError),
}

impl From<CalcSsaXFError> for PhotonAlphError {
    fn from(err: CalcSsaXFError) -> Self {
        PhotonAlphError::CalcSsaXF(err)
    }
}

/// Computes the absorption coefficient `alpha_nu` for synchrotron self-
/// absorption, using Rybicki & Lightman's Eq. 6.50.
/// All quantities are assumed to be in the local rest frame of the plasma,
/// and the electron distribution is assumed to be isotropic.
///
/// # Arguments
///
/// * `input` - struct containing the input parameters for the computation
///   * `num_dp` - number of momentum bins in the distribution of particles
///   * `b_sq` - square of magnetic field strength
///   * `i_params` - array (length 7) of up to 7 integer parameters,
///     used for debugging purposes
///   * `dn_pf` - particle distribution. This is the density of particles
///     in each bin, NOT number of particles per energy (`dn/dE`).
///     It is also not a pure number count `N` as used in other versions
///     of this subroutine
///   * `elec_gam_array` - Lorentz factor boundary values of the `dn(p)`
///     distribution
///   * `phot_en_array` - array of energy values for photons
///   * `o_o_nu_sq` - `1 / nu^2` for the energies given in
///     `phot_en_array`. Saves repeated divisions in the loops
///
/// # Returns
///
/// * `alpha_nu` - absorption coefficient at the frequencies given by
///   entries in `phot_en_array`
pub fn photon_alph(input: PhotonAlphInput<'_>) -> Result<Vec<f64>, PhotonAlphError> {
    let PhotonAlphInput {
        num_dp,
        elec_gam_array,
        dn_pf,
        b_sq,
        phot_en_array,
        o_o_nu_sq,
        i_params,
        globals,
    } = input;

    let num_ph = globals.num_ph;

    // Do some rusty error checking that isn't in the Fortran code.
    if elec_gam_array.len() <= num_dp {
        return Err(PhotonAlphError::ElectronGammaArrayTooShort {
            len: elec_gam_array.len(),
            required_min_len: num_dp + 1,
        });
    }

    if dn_pf.len() < num_dp {
        return Err(PhotonAlphError::ElectronDistributionTooShort {
            len: dn_pf.len(),
            required_min_len: num_dp,
        });
    }

    if phot_en_array.len() < num_ph {
        return Err(PhotonAlphError::PhotonEnergyArrayTooShort {
            len: phot_en_array.len(),
            required_min_len: num_ph,
        });
    }

    if o_o_nu_sq.len() < num_ph {
        return Err(PhotonAlphError::OneOverNuSqArrayTooShort {
            len: o_o_nu_sq.len(),
            required_min_len: num_ph,
        });
    }

    if globals.syn_x.len() <= num_ph {
        return Err(PhotonAlphError::SynXArrayTooShort {
            len: globals.syn_x.len(),
            required_min_len: num_ph + 1,
        });
    }

    if globals.syn_f.len() <= num_ph {
        return Err(PhotonAlphError::SynFArrayTooShort {
            len: globals.syn_f.len(),
            required_min_len: num_ph + 1,
        });
    }

    // Compute the magnetic field strength at this location, and the prefactor
    // for synchrotron emission (Rybicki & Lightman, eq. 6.18).
    // Units of p_fac power per unit frequency *per electron* (cgs units), and
    // the factor of sin(alpha) is not present -- it's handled by the
    // "isotropy correction" in calc_syn_Fx.
    // o_o_nuc uses Rybicki & Lightman, Eq. 6.17c, without the electron Lorentz
    // factors; those will be included in the loop over momentum
    let b_loc = b_sq.sqrt();
    let syn_prefac = SQRT_3 * QCGS.powi(3) * b_loc / (2.0 * PII * RM_ELEC);
    let o_o_nuc = 4.0 * PII * XME * CCGS / (3.0 * QCGS * b_loc);

    // alph_prefac, the collection of terms in RL79's Eq. 6.50 that don't
    // depend on electron energy
    let mut alph_prefac = vec![1.0e-99; num_ph];

    for m_ph in 0..num_ph {
        alph_prefac[m_ph] = -CCGS.powi(2) / (8.0 * PII) * syn_prefac * o_o_nu_sq[m_ph];
    }

    // Initialize alpha_nu prior to the main loop
    let mut alpha_nu = vec![1.0e-99; num_ph];

    //-------------------------------------------------------------------------
    // If the electron distribution is sufficiently delta-function-like, then
    // we can skip the computation completely and compute alpha_nu
    // analytically
    // Delta function is defined here as 1 or 2 filled bins in the rebinning;
    // with 3+ bins in a rebinned distribution, there is at least one bin
    // whose 2nd derivative will be computed correctly
    //-------------------------------------------------------------------------

    // Same behavior as Fortran maxloc(dn_pf, 1): first max wins.
    //
    // Fortran:
    //   i_max is 1-based for dn_pf(1:num_dp)
    //
    // Rust:
    //   i_max is 0-based for dn_pf[0..num_dp]
    let i_max = dn_pf
        .iter()
        .enumerate()
        .fold(
            (0usize, f64::NEG_INFINITY),
            |(best_i, best_value), (i, &value)| {
                // if value is STRICTLY greater,
                // update value and index (gives us first max if there's a tie)
                if value > best_value {
                    (i, value)
                } else {
                    (best_i, best_value)
                }
            },
        )
        .0;

    let filled_bin_count = dn_pf.iter().filter(|value| **value > 1.0e-55).count();

    if filled_bin_count <= 2 {
        // Fortran i_max is 1-based; Rust i_max is 0-based.
        // Fortran elec_gam_array(i_max) maps to Rust elec_gam_array[i_max + 1],
        // because elec_gam_array is a boundary array corresponding to 0:num_dp.
        let gamma_at_max = elec_gam_array[i_max + 1];

        // If electrons are sufficiently low in energy, they neither emit nor
        // absorb synchrotron radiation
        if gamma_at_max < 6.0 {
            return Ok(alpha_nu);
        }

        let mut i_x_lo = 0usize;
        let den_tot: f64 = dn_pf.iter().sum();

        for m_ph in 0..num_ph {
            // Compute x (RL79 Eq. 6.31c) for this photon energy
            let xxx = phot_en_array[m_ph] * o_o_nuc / (XH * gamma_at_max.powi(2));

            // If we are too far in the exponential tail of syn_F, ignore this
            // photon energy and move to the next
            if xxx >= 30.0 {
                alpha_nu[m_ph] = 1.0e-99;
                continue;
            }

            // Find the appropriate value of syn_x for this photon energy
            for i_x in i_x_lo..num_ph {
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

            // Combine x_F with other factors to get alpha_nu for this photon
            // energy, recalling that the prefactor for a delta function is 4/3
            // rather than (p+2) as originally used here.
            alpha_nu[m_ph] = 4.0 * THIRD * den_tot * syn_prefac * x_f
                / (8.0 * PII * XME * gamma_at_max * (phot_en_array[m_ph] / XH).powi(2));

            if alpha_nu[m_ph] < 1.0e-55 {
                alpha_nu[m_ph] = 1.0e-99;
            }
        }

        return Ok(alpha_nu);
    }

    //-------------------------------------------------------------------------
    // Delta-function possibility dealt with
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // Precompute the array pertaining to the electron distribution, since that
    // won't change over the integration
    // We'll need the electron energies during the integration, so switch from
    // Lorentz factor to energy here
    // We'll also need the number of electrons per unit energy, i.e. dN/dE
    //-------------------------------------------------------------------------

    // Array of electron energies
    let mut elec_en_array = vec![0.0; num_dp];

    let mut dnde = vec![0.0; num_dp];
    let mut d2nde2 = vec![0.0; num_dp];
    let mut diff_term = vec![0.0; num_dp];

    for l_mom in 0..num_dp {
        elec_en_array[l_mom] = (elec_gam_array[l_mom] * elec_gam_array[l_mom + 1]).sqrt() * RM_ELEC;

        dnde[l_mom] = dn_pf[l_mom] * O_O_RME / (elec_gam_array[l_mom + 1] - elec_gam_array[l_mom]);
    }

    d2nde2[0] = (dnde[1] - dnde[0]) / (elec_en_array[1] - elec_en_array[0]);

    for l_mom in 1..num_dp - 1 {
        d2nde2[l_mom] = (dnde[l_mom + 1] - dnde[l_mom - 1])
            / (elec_en_array[l_mom + 1] - elec_en_array[l_mom - 1]);
    }

    d2nde2[num_dp - 1] = (dnde[num_dp - 1] - dnde[num_dp - 2])
        / (elec_en_array[num_dp - 1] - elec_en_array[num_dp - 2]);

    for l_mom in 0..num_dp {
        diff_term[l_mom] = d2nde2[l_mom] - 2.0 * dnde[l_mom] / elec_en_array[l_mom];
    }

    //-------------------------------------------------------------------------
    // Array precomputation finished
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // Loop over photon energies, finding alpha_nu at each frequency
    //-------------------------------------------------------------------------

    for m_ph in 0..num_ph {
        //-----------------------------------------------------------------------
        // Integrate over the electron distribution to get alpha_nu, using the
        // trapezoidal rule.
        // For each boundary between bins, take the following steps
        // (1) Compute xxx: the photon frequency divided by the critical
        // frequency, or (E_nu/h)*(1/nu_c). Recall that we've already
        // computed most of 1/nu_c, and just need to incorporate the electron
        // Lorentz factor
        // (2) Find d^2n/dE^2, the derivative of dn/dE at that point
        // (3) Find x_F, the factor corresponding to syncrotron power radiated
        // at the current frequency by an electron at the current energy
        // (4) Combine those quantities to get the integrand at that location
        //
        // Do the low end of the first bin explicitly outside the loop so that
        // we can do *exactly* the same thing every time through
        //-----------------------------------------------------------------------

        let xxx_array: Vec<f64> = elec_gam_array
            .iter()
            .take(num_dp + 1)
            .map(|gamma| phot_en_array[m_ph] * o_o_nuc / (XH * gamma.powi(2)))
            .collect();

        let xxx_lo = xxx_array[1];

        let result = calc_ssa_x_f(xxx_lo, true, 117, globals)?;
        let x_f_lo = result.x_f;
        let mut n_x_out = result.n_x_out;
        let mut integrand_lo = x_f_lo * diff_term[0];

        let mut integral = 0.0;

        for l_mom in 1..=num_dp - 1 {
            // If electrons are sufficiently low in energy, they neither emit nor
            // absorb synchrotron radiation
            if elec_gam_array[l_mom] < 6.0 {
                integrand_lo = 1.0e-99;
                continue;
            }

            // Compute xxx_hi
            let xxx_hi = xxx_array[l_mom + 1];

            // Compute x_F_hi
            let n_x_in = n_x_out;
            let result = calc_ssa_x_f(xxx_hi, false, n_x_in, globals)?;
            let x_f_hi = result.x_f;
            let new_n_x_out = result.n_x_out;
            n_x_out = new_n_x_out;

            // Compute integrand_hi
            let integrand_hi = x_f_hi * diff_term[l_mom];

            // Compute contribution to integral for this bin
            integral += 0.5
                * (integrand_lo + integrand_hi)
                * (elec_en_array[l_mom] - elec_en_array[l_mom - 1]);

            // Set integrand_lo for the next cycle
            integrand_lo = integrand_hi;
        }

        //-----------------------------------------------------------------------
        // Integration finished
        //-----------------------------------------------------------------------

        // Take the integral and re-incorporate the prefactor. Also, if the
        // absorption coefficient isn't sufficiently positive, set it to roughly
        // zero. And since the absorption coefficient decreases with increasing
        // photon energy, we can also skip the rest of the computation
        alpha_nu[m_ph] = alph_prefac[m_ph] * integral;

        // Just a sanity check in case something unexpected happens
        if alpha_nu[m_ph] < 0.0 {
            return Err(PhotonAlphError::NegativeAlpha {
                m_ph,
                alpha_nu: alpha_nu[m_ph],
            });
        }

        if alpha_nu[m_ph] < 1.0e-55 {
            alpha_nu[m_ph] = 1.0e-99;
            break;
        }
    }

    //-------------------------------------------------------------------------
    // alpha_nu computed for all photon frequencies
    //-------------------------------------------------------------------------

    let _ = i_params;

    Ok(alpha_nu)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_globals(num_ph: usize) -> Globals {
        let mut syn_x = vec![0.0; 119];
        let mut syn_f = vec![0.0; 119];

        for i in 0..119 {
            syn_x[i] = i as f64;
            syn_f[i] = 1.0 + i as f64;
        }

        Globals {
            num_ph,
            syn_x,
            syn_f,
            o_o_sone: 1.0,
            ..Default::default()
        }
    }

    fn base_input<'a>(
        num_dp: usize,
        elec_gam_array: &'a [f64],
        dn_pf: &'a [f64],
        phot_en_array: &'a [f64],
        o_o_nu_sq: &'a [f64],
        i_params: &'a [i32; 7],
        globals: &'a Globals,
    ) -> PhotonAlphInput<'a> {
        PhotonAlphInput {
            num_dp,
            elec_gam_array,
            dn_pf,
            b_sq: 1.0,
            phot_en_array,
            o_o_nu_sq,
            i_params,
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
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonAlphError::ElectronGammaArrayTooShort {
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
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonAlphError::ElectronDistributionTooShort {
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
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonAlphError::PhotonEnergyArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_o_o_nu_sq_array_is_too_short() {
        let num_dp = 4;
        let globals = base_globals(4);
        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_nu_sq = vec![1.0; globals.num_ph - 1];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonAlphError::OneOverNuSqArrayTooShort {
                len: 3,
                required_min_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_syn_x_array_is_too_short() {
        let num_dp = 4;
        let mut globals = base_globals(4);
        globals.syn_x.truncate(globals.num_ph);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonAlphError::SynXArrayTooShort {
                len: 4,
                required_min_len: 5
            })
        ));
    }

    #[test]
    fn returns_error_when_syn_f_array_is_too_short() {
        let num_dp = 4;
        let mut globals = base_globals(4);
        globals.syn_f.truncate(globals.num_ph);

        let elec_gam_array = vec![10.0; num_dp + 1];
        let dn_pf = vec![1.0; num_dp];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ));

        assert!(matches!(
            result,
            Err(PhotonAlphError::SynFArrayTooShort {
                len: 4,
                required_min_len: 5
            })
        ));
    }

    #[test]
    fn delta_like_low_energy_distribution_returns_floor_values() {
        let num_dp = 4;
        let globals = base_globals(4);
        let elec_gam_array = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let dn_pf = vec![1.0, 0.0, 0.0, 0.0];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ))
        .expect("photon_alph should succeed");

        assert_eq!(result.len(), globals.num_ph);
        assert!(
            result
                .iter()
                .all(|value| (*value - 1.0e-99).abs() < 1.0e-120)
        );
    }

    #[test]
    fn delta_like_distribution_returns_finite_values() {
        let num_dp = 4;
        let globals = base_globals(4);
        let elec_gam_array = vec![1.0, 10.0, 20.0, 30.0, 40.0];
        let dn_pf = vec![1.0, 0.0, 0.0, 0.0];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ))
        .expect("photon_alph should succeed");

        assert_eq!(result.len(), globals.num_ph);
        assert!(result.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn broad_distribution_returns_finite_values() {
        let num_dp = 4;
        let globals = base_globals(4);
        let elec_gam_array = vec![1.0, 10.0, 20.0, 30.0, 40.0];
        let dn_pf = vec![1.0, 2.0, 3.0, 4.0];
        let phot_en_array = vec![1.0e-12; globals.num_ph];
        let o_o_nu_sq = vec![1.0; globals.num_ph];
        let i_params = [0; 7];

        let result = photon_alph(base_input(
            num_dp,
            &elec_gam_array,
            &dn_pf,
            &phot_en_array,
            &o_o_nu_sq,
            &i_params,
            &globals,
        ))
        .expect("photon_alph should succeed");

        assert_eq!(result.len(), globals.num_ph);
        assert!(result.iter().all(|value| value.is_finite()));
    }
}
