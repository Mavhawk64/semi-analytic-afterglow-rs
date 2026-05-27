use crate::functions::fixed_rperp::fixed_rperp;
use crate::modules::inputs::Inputs;

/// Numerically computes the derivative of:
///
/// R_perp = R * sqrt(1 - mu^2)
///
/// with respect to mu using a centered finite difference method.
///
/// # Arguments
///
/// * `mu_input` - cosine of angle between observer and fluid parcel
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters to be used in evaluating the function
/// * `inputs` - reference to the Inputs struct containing various input parameters
pub fn d_rperp_dmu(
    mu_input: f64,
    i_params: &[i32; 7],
    r_params: &[f64; 7],
    inputs: &Inputs,
) -> f64 {
    // Compute R_perp just above and below the input value of mu
    // Rust note: fix a small edge case where mu_input is 0.0.
    let delta = mu_delta(mu_input);

    let mu_lo = mu_input - delta;
    let mu_hi = mu_input + delta;

    let rperp_lo = -fixed_rperp(mu_lo, i_params, r_params, inputs);
    let rperp_hi = -fixed_rperp(mu_hi, i_params, r_params, inputs);

    // Compute the derivative
    (rperp_hi - rperp_lo) / (mu_hi - mu_lo)
}

fn mu_delta(mu_input: f64) -> f64 {
    if mu_input == 0.0 {
        1.0e-6
    } else {
        mu_input.abs() * 1.0e-6
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_inputs(k_cbm: i32) -> Inputs {
        Inputs {
            k_cbm,
            do_jet_structure: false,
            e_iso: 1.0,
            ..Default::default()
        }
    }

    fn base_params(lambda: f64, ctz: f64, r_obs: f64, targ: f64) -> ([i32; 7], [f64; 7]) {
        let i_params = [0; 7];

        let mut r_params = [0.0; 7];
        r_params[0] = 0.1; // theta_obs
        r_params[1] = 0.2; // psi_fl
        r_params[2] = lambda;
        r_params[3] = ctz;
        r_params[4] = r_obs;
        r_params[5] = targ;

        (i_params, r_params)
    }

    #[test]
    fn mu_delta_handles_zero_without_returning_zero() {
        assert_eq!(mu_delta(0.0), 1.0e-6);
    }

    #[test]
    fn derivative_is_finite_for_typical_inputs() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0, 10.0, 5.0);

        let derivative = d_rperp_dmu(0.4, &i_params, &r_params, &inputs);

        assert!(derivative.is_finite());
    }

    #[test]
    fn derivative_changes_with_mu() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0, 10.0, 5.0);

        let derivative_1 = d_rperp_dmu(0.3, &i_params, &r_params, &inputs);
        let derivative_2 = d_rperp_dmu(0.7, &i_params, &r_params, &inputs);

        assert_ne!(derivative_1, derivative_2);
    }

    #[test]
    fn derivative_matches_manual_centered_difference() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0, 10.0, 5.0);

        let mu_input: f64 = 0.4;

        let delta = mu_delta(mu_input);

        let mu_lo = mu_input - delta;
        let mu_hi = mu_input + delta;

        let rperp_lo = -fixed_rperp(mu_lo, &i_params, &r_params, &inputs);
        let rperp_hi = -fixed_rperp(mu_hi, &i_params, &r_params, &inputs);

        let expected = (rperp_hi - rperp_lo) / (mu_hi - mu_lo);
        let actual = d_rperp_dmu(mu_input, &i_params, &r_params, &inputs);

        assert!((actual - expected).abs() < 1.0e-12);
    }

    #[test]
    fn derivative_is_finite_at_mu_zero() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0, 10.0, 5.0);

        let derivative = d_rperp_dmu(0.0, &i_params, &r_params, &inputs);

        assert!(derivative.is_finite());
    }
}
