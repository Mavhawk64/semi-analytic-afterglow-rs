use crate::functions::r_of_theta::r_of_theta;
use crate::modules::inputs::Inputs;

/// Computes the value of X-R_perp, where X is a target value we want to hit,
/// and R_perp = R(theta) * sqrt(1 - mu^2) is the perpendicular component of
/// the fluid particle's position relative to the line of sight.
///
/// This setup with X = 0 allows us to use this function to find the maximum
/// value of R_perp, since that minimizes X - R_perp for Brent's algorithm.
///
/// Returns fixed_rperp, negative of perpendicular component of radius.
///
/// # Arguments
///
/// * `mu_input` - cosine of angle between observer and fluid parcel
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters to be used in evaluating the function
/// * `inputs` - reference to the Inputs struct containing various input parameters
pub fn fixed_rperp(
    mu_input: f64,
    i_params: &[i32; 7],
    r_params: &[f64; 7],
    inputs: &Inputs,
) -> f64 {
    // Local variables

    // Get two variables from the input parameters
    //
    // Fortran used r_params(5) and r_params(6).
    // Rust is 0-indexed, so these become [4] and [5].
    let r_obs: f64 = r_params[4];
    let targ: f64 = r_params[5];

    // Compute R_of_theta
    let r_fl: f64 = if mu_input == 1.0 {
        r_obs
    } else {
        r_of_theta(mu_input, i_params, r_params, inputs)
    };

    // determine X-R_perp and return it
    targ - r_fl * (1.0 - mu_input.powi(2)).sqrt()
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
    fn returns_target_when_mu_is_one() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0, 10.0, 5.0);

        let result = fixed_rperp(1.0, &i_params, &r_params, &inputs);

        assert_eq!(result, 5.0);
    }

    #[test]
    fn computes_target_minus_perpendicular_radius_for_k2_case() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0, 10.0, 5.0);
        let mu_input = 0.4;

        let r_fl = r_of_theta(mu_input, &i_params, &r_params, &inputs);
        let expected = r_params[5] - r_fl * (1.0 - mu_input.powi(2)).sqrt();

        let result = fixed_rperp(mu_input, &i_params, &r_params, &inputs);

        assert!((result - expected).abs() < 1.0e-12);
    }
}
