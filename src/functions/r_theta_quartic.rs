/// Evaluates the quartic function R(theta) with the coefficients precomputed.
/// Needed only when the quartic formula returns senseless values in subroutine
/// R_of_theta.
///
/// Returns R(theta), the value of the quartic function.
///
/// # Arguments
///
/// * `r_input` - current value being tried for R
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used
///   in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters to
///   be used in evaluating the function
pub fn r_theta_quartic(r_input: f64, _i_params: &[i32; 7], r_params: &[f64; 7]) -> f64 {
    // Grab the coefficients from the inputs
    //
    // Fortran:
    //   r_params(1), r_params(2)
    //
    // Rust:
    //   r_params[0], r_params[1]
    let a = r_params[0];
    let b = r_params[1];

    r_input.powi(4) + a * r_input + b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates_quartic_with_precomputed_coefficients() {
        let i_params = [0; 7];

        let mut r_params = [0.0; 7];
        r_params[0] = 2.0; // a
        r_params[1] = -3.0; // b

        let result = r_theta_quartic(4.0, &i_params, &r_params);

        // R^4 + aR + b = 4^4 + 2*4 - 3 = 256 + 8 - 3 = 261
        assert_eq!(result, 261.0);
    }
}
