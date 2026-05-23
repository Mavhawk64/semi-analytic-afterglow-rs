use crate::modules::constants::{CCGS, RM_ELEC, RM_PROT, XME, XMP};
use crate::modules::inputs::Inputs;
use crate::subroutines::calc_epse::calc_epse;

#[derive(Debug, Clone, PartialEq)]
pub enum SetElecDistError {
    MomentumArrayTooShort { len: usize, required_min_len: usize },
    ElectronGammaArrayTooShort { len: usize, required_min_len: usize },
    EmptyDistribution { reason: &'static str },
}

/// Create the before-losses electron distribution. For this, we need four
/// pieces of information:
///
/// 1. Spectral index of the non-thermal tail
/// 2. Location of the thermal peak
/// 3. Normalization of the non-thermal tail, as fraction of the total electron distribution
/// 4. Location of cooling cutoff
///
/// If we are only using a power law, then "the location of the thermal peak"
/// refers instead to the minimum energy of the power law.
///
/// Arguments:
/// * `gam_shock` - Lorentz factor of the shock when these particles were swept up
/// * `n_ext` - density into which the shock is expanding
/// * `p_pf_cgs` - array holding default (uncooled) momentum values
/// * `elec_gam_array` - list of Lorentz factors associated with p_pf_cgs
/// * `num_dp` - number of bins in particle distribution--not all of these will be filled
/// * `inputs` - struct containing user inputs
pub fn set_elec_dist(
    gam_shock: f64,
    n_ext: f64,
    p_pf_cgs: &[f64],
    elec_gam_array: &[f64],
    num_dp: usize,
    inputs: &Inputs,
) -> Result<Vec<f64>, SetElecDistError> {
    // Parameters for computing p_max
    const P_C: f64 = 6.35e5;
    const P_EPS: f64 = 0.0;
    const G_C: f64 = 0.532;
    const G_EPS: f64 = -0.338;
    const L_C: f64 = 2.12;
    const L_EPS: f64 = 0.0;
    const H_C: f64 = 1.01;
    const H_EPS: f64 = 0.0;
    const W_C: f64 = 5.26;
    const W_EPS: f64 = 0.11;

    // Check that the input arrays are long enough to hold the distribution.
    // Not in original Fortran code.
    if p_pf_cgs.len() <= num_dp {
        return Err(SetElecDistError::MomentumArrayTooShort {
            len: p_pf_cgs.len(),
            required_min_len: num_dp + 1,
        });
    }
    if elec_gam_array.len() <= num_dp {
        return Err(SetElecDistError::ElectronGammaArrayTooShort {
            len: elec_gam_array.len(),
            required_min_len: num_dp + 1,
        });
    }

    // Minimum Lorentz factor of all electrons (really, electron bins)
    let mut gam_min = 0.0;

    let mut mj_a = 0.0;
    let mut o_o_mja = 0.0;

    let mut i_cross = 1usize;

    // Maxwell-Juttner array
    let mut mj_array = vec![1.0e-99; num_dp];

    // Power-law array
    let mut pl_array = vec![1.0e-99; num_dp];

    // Thermal distribution values
    let mut th_dist = vec![1.0e-99; num_dp];

    // Non-thermal distribution values
    let mut nt_dist = vec![1.0e-99; num_dp];

    // (1) Spectral index; while this may eventually depend on shock speed,
    // for now just assume a fixed value
    // (Spectral index)
    let sigma: f64 = inputs.p_spec;

    // (2) Location of the thermal peak comes from bulk kinetic energy flux,
    // and the fraction thereof transferred to electrons
    // (Fraction of energy in electrons)
    let epse: f64 = calc_epse(gam_shock, inputs);

    if inputs.pl_only {
        // Granot & Sari (2002), Eqs. (A1), (A4), and (A6)
        gam_min =
            epse * (sigma - 2.0) / (sigma - 1.0) * gam_shock * RM_PROT / (2.0_f64.sqrt() * RM_ELEC);
    } else {
        mj_a = epse * gam_shock * XMP / (3.0 * XME);
        o_o_mja = 1.0 / mj_a;
    }

    // (3) Normalization of the non-thermal tail, which is set to 1.0 if only
    // a power-law is used
    let nt_norm = if inputs.pl_only { 1.0 } else { inputs.nt_norm };

    // (4) Location of cooling cutoff comes from fits to Monte Carlo sims.
    // See Warren et al. (2021ApJ...906...33W) for the origin of these
    // formulae. We assume here that eta_mfp = 1 -- this may be revisited
    // at a later date.
    let p_max = if inputs.do_pmax {
        let p_pk = P_C * (inputs.eps_b0 * n_ext).powf(P_EPS) * XMP * CCGS;
        let gb_pk = G_C * (inputs.eps_b0 * n_ext).powf(G_EPS);
        // Asymptotic power-law indices far below and above (respectively) of
        // the peak, as defined in Warren (2021, DOI: 10.3847/1538-4357/abc694)
        let l = L_C * (inputs.eps_b0 * n_ext).powf(L_EPS);
        let h = H_C * (inputs.eps_b0 * n_ext).powf(H_EPS);
        // Control parameter for the width of the peak, as defined in
        // Warren (2021, DOI: 10.3847/1538-4357/abc694)
        let w = W_C * (inputs.eps_b0 * n_ext).powf(W_EPS);
        // Product of Gamma (Lorentz factor) and beta
        let gb = (gam_shock.powi(2) - 1.0).sqrt();

        // Eq. (8) of Warren+ (2021)
        let x = gb / gb_pk;

        p_pk * x.powf(l) * ((1.0 + l / h) / (1.0 + (l / h) * x.powf(w))).powf((l + h) / w)
    } else {
        // Set p_max very high
        6.0e14 * XME * CCGS
    };

    // Do we have enough of a distribution to see a power law? Don't adjust
    // anything even if we don't, but throw up a warning
    if inputs.pl_only {
        if (1.0e1 * gam_min * XME * CCGS) > p_max {
            // Temporary variable for calculating gamma from momentum before throwing an error
            let gamma_mom = (1.0 + p_max.powi(2) / (XME * CCGS).powi(2)).sqrt();

            eprintln!("Warning in set_elec_dist: gam_min is close to MJ_a");
            eprintln!(
                "      gam_min = {gam_min:11.3e}; p_max = {p_max:11.3e} (gam = {gamma_mom:11.3e})"
            );
            eprintln!("      Gam_shock = {gam_shock:11.3e}");
        }
    } else if (1.0e1 * mj_a * XME * CCGS) > p_max {
        // Temporary variable for calculating gamma from momentum before throwing an error
        let gamma_mom = (1.0 + p_max.powi(2) / (XME * CCGS).powi(2)).sqrt();

        eprintln!("Warning in set_elec_dist: p_max is close to MJ_a");
        eprintln!("      MJ_a = {mj_a:11.3e}; p_max = {p_max:11.3e} (gam = {gamma_mom:11.3e})");
        eprintln!("      Gam_shock = {gam_shock:11.3e}");
    }

    let o_o_pmax = 1.0 / p_max;

    //--------------------------------------------------------------------------
    // If we are using a thermal + power-law distribution, find the crossover
    // momentum where we shift from thermal to non-thermal parts. Said
    // momentum is defined by three equations
    // (1) f_nt = \int_{p_x}^{\infty} C_1 p^{-\sigma}\exp[-(\frac{p}{p_{\max}})^2] dp
    // (2) f_th = 1 - f_nt = \int_0^{p_x} C_2 p^2 \exp[-\frac{\gamma(p)}{a}] dp
    // (3) C_1 p_x^{-\sigma} \exp[-(\frac{p_x}{p_{\max}})^2] = C_2 p_x^2 \exp[-\frac{\gamma(p_x)}{a}]
    // To identify the crossover momentum, reduce the above equations to a
    // single one and step through until the sign of the difference switches.
    //--------------------------------------------------------------------------
    if !inputs.pl_only {
        // TODO: Verify these are 0-indexed. The indices align with the fortran code...
        //
        // Compute the partial sums for LHS integral running from zone 1 to the
        // top of the electron distribution
        let mut mj_lo = 0.0;
        let mut mj_hi = p_pf_cgs[1].powi(2) * (-elec_gam_array[1] * o_o_mja).exp();

        mj_array[0] = nt_norm * 0.5 * (mj_hi + mj_lo) * (p_pf_cgs[1] - p_pf_cgs[0]);

        for l_mom in 2..=num_dp {
            mj_lo = mj_hi;

            mj_hi = if elec_gam_array[l_mom] * o_o_mja > 228.0 {
                0.0
            } else {
                p_pf_cgs[l_mom].powi(2) * (-elec_gam_array[l_mom] * o_o_mja).exp()
            };

            mj_array[l_mom - 1] = mj_array[l_mom - 2]
                + nt_norm * 0.5 * (mj_hi + mj_lo) * (p_pf_cgs[l_mom] - p_pf_cgs[l_mom - 1]);
        }

        // Compute partial sums for the RHS integral, running from zone num_dp-1
        // (since we are integrating to top of array) *down* to zone 1
        let mut pl_lo =
            p_pf_cgs[num_dp - 1].powf(-sigma) * (-(p_pf_cgs[num_dp - 1] * o_o_pmax).powi(2)).exp();
        let mut pl_hi =
            p_pf_cgs[num_dp].powf(-sigma) * (-(p_pf_cgs[num_dp] * o_o_pmax).powi(2)).exp();

        pl_array[num_dp - 2] =
            (1.0 - nt_norm) * 0.5 * (pl_hi + pl_lo) * (p_pf_cgs[num_dp] - p_pf_cgs[num_dp - 1]);

        for l_mom in (1..=num_dp - 2).rev() {
            pl_hi = pl_lo;
            pl_lo = p_pf_cgs[l_mom].powf(-sigma) * (-(p_pf_cgs[l_mom] * o_o_pmax).powi(2)).exp();

            pl_array[l_mom - 1] = pl_array[l_mom]
                + (1.0 - nt_norm) * 0.5 * (pl_hi + pl_lo) * (p_pf_cgs[l_mom + 1] - p_pf_cgs[l_mom]);
        }

        // Once we have the array of partial sums, include the prefactor
        for l_mom in 1..=num_dp - 1 {
            pl_array[l_mom - 1] = pl_array[l_mom - 1]
                * p_pf_cgs[l_mom].powf(2.0 + sigma)
                * ((p_pf_cgs[l_mom] * o_o_pmax).powi(2) - elec_gam_array[l_mom] * o_o_mja).exp();

            if mj_array[l_mom - 1] >= pl_array[l_mom - 1] {
                i_cross = l_mom;
                break;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Crossover momentum found, if necessary
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // If we are using a thermal + power-law distribution, loop over momentum
    // and set un-normalized thermal and non-thermal distributions. Otherwise
    // just set the non-thermal distribution
    //--------------------------------------------------------------------------

    if !inputs.pl_only {
        // Thermal distribution
        th_dist.fill(1.0e-99);

        for l_mom in 1..i_cross {
            let p_mid = (p_pf_cgs[l_mom] * p_pf_cgs[l_mom - 1]).sqrt();
            let gam_mid = (elec_gam_array[l_mom] * elec_gam_array[l_mom - 1]).sqrt();

            th_dist[l_mom - 1] = p_mid.powi(2)
                * (-gam_mid * o_o_mja).exp()
                * (p_pf_cgs[l_mom] - p_pf_cgs[l_mom - 1]);
        }

        // Nonthermal tail
        nt_dist.fill(1.0e-99);

        for l_mom in i_cross..=num_dp {
            let p_mid = (p_pf_cgs[l_mom] * p_pf_cgs[l_mom - 1]).sqrt();

            nt_dist[l_mom - 1] = p_mid.powf(-sigma)
                * (-(p_mid * o_o_pmax).powi(2)).exp()
                * (p_pf_cgs[l_mom] - p_pf_cgs[l_mom - 1]);

            // Exit loop if we're seeing really low values
            if nt_dist[l_mom - 1] < 1.0e-65 {
                break;
            }
        }
    } else {
        // Zero out the thermal distribution
        th_dist.fill(1.0e-99);

        // Create the nonthermal distribution
        nt_dist.fill(1.0e-99);

        for l_mom in 1..=num_dp {
            if gam_min < elec_gam_array[l_mom] {
                let p_mid = (p_pf_cgs[l_mom] * p_pf_cgs[l_mom - 1]).sqrt();

                nt_dist[l_mom - 1] = p_mid.powf(-sigma)
                    * (-(p_mid * o_o_pmax).powi(2)).exp()
                    * (p_pf_cgs[l_mom] - p_pf_cgs[l_mom - 1]);

                // Exit loop if we're seeing really low values
                if nt_dist[l_mom - 1] < 1.0e-65 {
                    break;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Un-normalized distribution(s) found
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Normalize the two distributions and combine them into one output array
    //--------------------------------------------------------------------------

    let th_dist_norm: f64 = if !inputs.pl_only {
        th_dist.iter().sum()
    } else {
        0.0
    };

    let nt_dist_norm: f64 = nt_dist.iter().sum();

    if !inputs.pl_only {
        let norm_fac = (1.0 - nt_norm) / th_dist_norm;

        for value in &mut th_dist {
            if *value > 1.0e-55 {
                *value *= norm_fac;
            }
        }
    }

    let norm_fac = nt_norm / nt_dist_norm;

    for value in &mut nt_dist {
        if *value > 1.0e-55 {
            *value *= norm_fac;
        }
    }

    let mut elec_dist = vec![1.0e-99; num_dp];

    for i in 0..num_dp {
        if th_dist[i] > 1.0e-55 || nt_dist[i] > 1.0e-55 {
            elec_dist[i] = th_dist[i] + nt_dist[i];
        }
    }

    Ok(elec_dist)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_arrays(num_dp: usize) -> (Vec<f64>, Vec<f64>) {
        let p_pf_cgs: Vec<f64> = (0..=num_dp)
            .map(|i| 1.0e-4 * 10.0_f64.powf(i as f64 * 0.25))
            .collect();

        let elec_gam_array: Vec<f64> = p_pf_cgs
            .iter()
            .map(|p| (1.0 + (p / (XME * CCGS)).powi(2)).sqrt())
            .collect();

        (p_pf_cgs, elec_gam_array)
    }

    fn base_inputs(pl_only: bool) -> Inputs {
        Inputs {
            pl_only,
            p_spec: 2.2,
            epse0: 0.1,
            nt_norm: 0.1,
            eps_b0: 0.01,
            do_pmax: false,
            ..Default::default()
        }
    }

    #[test]
    fn returns_error_when_momentum_array_is_too_short() {
        let inputs = base_inputs(true);
        let num_dp = 4;

        let p_pf_cgs = vec![1.0; num_dp];
        let elec_gam_array = vec![1.0; num_dp + 1];

        let result = set_elec_dist(10.0, 1.0, &p_pf_cgs, &elec_gam_array, num_dp, &inputs);

        assert!(matches!(
            result,
            Err(SetElecDistError::MomentumArrayTooShort {
                len: 4,
                required_min_len: 5
            })
        ));
    }

    #[test]
    fn returns_error_when_electron_gamma_array_is_too_short() {
        let inputs = base_inputs(true);
        let num_dp = 4;

        let p_pf_cgs = vec![1.0; num_dp + 1];
        let elec_gam_array = vec![1.0; num_dp];

        let result = set_elec_dist(10.0, 1.0, &p_pf_cgs, &elec_gam_array, num_dp, &inputs);

        assert!(matches!(
            result,
            Err(SetElecDistError::ElectronGammaArrayTooShort {
                len: 4,
                required_min_len: 5
            })
        ));
    }

    #[test]
    fn returns_power_law_distribution_with_expected_length() {
        let inputs = base_inputs(true);
        let num_dp = 64;
        let (p_pf_cgs, elec_gam_array) = simple_arrays(num_dp);

        let result = set_elec_dist(10.0, 1.0, &p_pf_cgs, &elec_gam_array, num_dp, &inputs)
            .expect("set_elec_dist should return a distribution");

        assert_eq!(result.len(), num_dp);
    }

    #[test]
    fn power_law_distribution_is_finite() {
        let inputs = base_inputs(true);
        let num_dp = 64;
        let (p_pf_cgs, elec_gam_array) = simple_arrays(num_dp);

        let result = set_elec_dist(10.0, 1.0, &p_pf_cgs, &elec_gam_array, num_dp, &inputs)
            .expect("set_elec_dist should return a distribution");

        assert!(result.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn power_law_distribution_has_nonzero_physical_values() {
        let inputs = base_inputs(true);
        let num_dp = 64;
        let (p_pf_cgs, elec_gam_array) = simple_arrays(num_dp);

        let result = set_elec_dist(10.0, 1.0, &p_pf_cgs, &elec_gam_array, num_dp, &inputs)
            .expect("set_elec_dist should return a distribution");

        assert!(result.iter().any(|value| *value > 1.0e-55));
    }

    #[test]
    fn power_law_distribution_is_normalized_close_to_one() {
        let inputs = base_inputs(true);
        let num_dp = 64;
        let (p_pf_cgs, elec_gam_array) = simple_arrays(num_dp);

        let result = set_elec_dist(10.0, 1.0, &p_pf_cgs, &elec_gam_array, num_dp, &inputs)
            .expect("set_elec_dist should return a distribution");

        let sum: f64 = result.iter().sum();

        assert!((sum - 1.0).abs() < 1.0e-10, "sum = {sum}");
    }
}
