use crate::functions::e_of_theta::e_of_theta;
use crate::functions::r_theta_quartic::r_theta_quartic;
use crate::functions::theta_f::theta_f;
use crate::modules::constants::CCGS;
use crate::modules::inputs::Inputs;
use crate::subroutines::newtons_method::newtons_method;

/// External-facing function (does Newton's Method approach)
pub fn r_of_theta(mu_input: f64, i_params: &[i32; 7], r_params: &[f64; 7], inputs: &Inputs) -> f64 {
    r_of_theta_internal(mu_input, i_params, r_params, inputs, true)
}

/// Solves for R as a function of theta_fl. Takes different path if k = 0
/// compared to if k = 2, since the latter can be solved using the quadratic
/// formula.
///
/// Returns the difference between the right and left sides of the
/// equation -- obviously we want that to be 0.
/// # Arguments
///
/// * `mu_input` - current value being tried for mu
/// * `i_params` - array (length 7) of up to 7 integer parameters to be
///   used in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters
///   to be used in evaluating the function
/// * `inputs` - runtime input parameters, replacing the Fortran module-level `k_cbm`
/// * `do_newtons_method` - flag to determine whether to use Newton's method
pub fn r_of_theta_internal(
    mu_input: f64,
    i_params: &[i32; 7],
    r_params: &[f64; 7],
    inputs: &Inputs,
    do_newtons_method: bool,
) -> f64 {
    // Load in data from input arrays
    //
    // Fortran:
    //   r_params(1:4)
    //
    // Rust:
    //   r_params[0:3]
    let theta_obs: f64 = r_params[0];
    let psi_fl: f64 = r_params[1];
    let lambda: f64 = r_params[2];
    let ctz: f64 = r_params[3];

    // Compute theta_fl given mu, theta_obs, and psi_fl.
    // Also compute E_of_theta
    let theta_fl: f64 = theta_f(theta_obs, psi_fl, mu_input);
    let e_theta: f64 = e_of_theta(theta_fl, inputs);

    // If k = 2, then equation is quadratic and may be solved very easily --
    // especially since we know we need the positive root
    if inputs.k_cbm == 2 {
        let a: f64 = 1.0 / (4.0 * lambda * e_theta);
        let b: f64 = 1.0 - mu_input;
        let c: f64 = -ctz;

        let mut r_value: f64 = -b + (b.powi(2) - 4.0 * a * c).sqrt();
        r_value /= 2.0 * a;

        r_value

    // If k = 0, then the equation is quartic. This means much more work;
    // use Newton's method and numerically find the root rather than the
    // quartic formula because the latter has proven to be rather finicky.
    } else if !do_newtons_method {
        let a: f64 = 8.0 * e_theta * lambda * (1.0 - mu_input);
        let b: f64 = -8.0 * e_theta * lambda * ctz; // verify this is ct_obs / (1+z)
        let p: f64 = -b;
        let q: f64 = -a * a / 8.0;
        let m: f64 = (-q / 2.0 + (q * q / 4.0 + p * p * p / 27.0).sqrt()).cbrt()
            + (-q / 2.0 - (q * q / 4.0 + p * p * p / 27.0).sqrt()).cbrt();
        if a == 0.0 {
            // (-B)^(1/4)
            f64::powf(-b, 0.25)
        } else if a > 0.0 {
            // -sqrt(m/2) + sqrt(-(m/2 - A/(sqrt(8m))))
            -(m / 2.0).sqrt() + (-(m / 2.0 - a / (8.0 * m).sqrt())).sqrt()
        } else {
            // Just in case (should never happen)
            // sqrt(m/2) + sqrt(-(m/2 + A/(sqrt(8m))))
            (m / 2.0).sqrt() + (-(m / 2.0 + a / (8.0 * m).sqrt())).sqrt()
        }
    } else {
        // Coefficients of original quartic function,
        // R^4 + A*R + B = 0
        let a: f64 = (1.0 - mu_input) * 8.0 * lambda * e_theta;
        let b: f64 = -ctz * 8.0 * lambda * e_theta;

        let mut i_params_int = [0i32; 7];
        let mut r_params_int = [0.0f64; 7];

        // on-axis or not
        i_params_int[4] = i_params[4];

        r_params_int[0] = a;
        r_params_int[1] = b;

        newtons_method(
            r_theta_quartic,
            &i_params_int,
            &r_params_int,
            1.0e8 * CCGS,
            CCGS,
            1.0e-2 * f64::MAX.powf(0.25),
        )
        .expect("Newton's method failed in r_of_theta")
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

    fn base_params(lambda: f64, ctz: f64) -> ([i32; 7], [f64; 7]) {
        let i_params = [0; 7];

        let mut r_params = [0.0; 7];
        r_params[0] = 0.1; // theta_obs
        r_params[1] = 0.2; // psi_fl
        r_params[2] = lambda;
        r_params[3] = ctz;

        (i_params, r_params)
    }

    fn quartic_coefficients(mu_input: f64, r_params: &[f64; 7], inputs: &Inputs) -> (f64, f64) {
        let theta_fl = theta_f(r_params[0], r_params[1], mu_input);
        let e_theta = e_of_theta(theta_fl, inputs);

        let a = (1.0 - mu_input) * 8.0 * r_params[2] * e_theta;
        let b = -r_params[3] * 8.0 * r_params[2] * e_theta;

        (a, b)
    }

    fn quartic_residual(r: f64, a: f64, b: f64) -> f64 {
        r.powi(4) + a * r + b
    }

    fn relative_quartic_residual(r: f64, a: f64, b: f64) -> f64 {
        let residual = quartic_residual(r, a, b);
        let scale = r.powi(4).abs() + (a * r).abs() + b.abs();

        if scale == 0.0 {
            residual.abs()
        } else {
            residual.abs() / scale
        }
    }

    #[test]
    fn k2_quadratic_solution_satisfies_original_equation() {
        let inputs = base_inputs(2);
        let (i_params, r_params) = base_params(2.0, 3.0);
        let mu_input = 0.4;

        let r = r_of_theta(mu_input, &i_params, &r_params, &inputs);

        let theta_fl = theta_f(r_params[0], r_params[1], mu_input);
        let e_theta = e_of_theta(theta_fl, &inputs);
        let lambda = r_params[2];
        let ctz = r_params[3];

        let residual = r.powi(2) / (4.0 * lambda * e_theta) + (1.0 - mu_input) * r - ctz;

        assert!(residual.abs() < 1.0e-8, "r = {r}, residual = {residual}");
    }

    #[test]
    fn k0_newton_solution_satisfies_quartic_equation() {
        let inputs = base_inputs(0);
        let (i_params, r_params) = base_params(1.0e64, 1.0e18);
        let mu_input = 0.5;

        let r = r_of_theta_internal(mu_input, &i_params, &r_params, &inputs, true);

        let (a, b) = quartic_coefficients(mu_input, &r_params, &inputs);
        let relative_residual = relative_quartic_residual(r, a, b);

        assert!(
            relative_residual < 1.0e-8,
            "r = {r}, a = {a}, b = {b}, relative_residual = {relative_residual}"
        );
    }

    #[test]
    fn k0_quartic_formula_solution_satisfies_quartic_equation() {
        let inputs = base_inputs(0);
        let (i_params, r_params) = base_params(1.0, 1.0);
        let mu_input = 0.5;

        let r = r_of_theta_internal(mu_input, &i_params, &r_params, &inputs, false);

        let (a, b) = quartic_coefficients(mu_input, &r_params, &inputs);
        let relative_residual = relative_quartic_residual(r, a, b);

        assert!(r.is_finite(), "quartic formula returned non-finite r = {r}");
        assert!(
            relative_residual < 1.0e-8,
            "r = {r}, a = {a}, b = {b}, relative_residual = {relative_residual}"
        );
    }

    #[test]
    fn k0_newton_and_quartic_formula_are_close_for_safe_case() {
        let inputs = base_inputs(0);
        let (i_params, r_params) = base_params(1.0e64, 1.0e18);
        let mu_input = 0.5;

        let r_newton = r_of_theta_internal(mu_input, &i_params, &r_params, &inputs, true);
        let r_quartic = r_of_theta_internal(mu_input, &i_params, &r_params, &inputs, false);

        let relative_difference =
            (r_newton - r_quartic).abs() / r_newton.abs().max(r_quartic.abs());

        assert!(
            relative_difference < 1.0e-6,
            "r_newton = {r_newton}, r_quartic = {r_quartic}, relative_difference = {relative_difference}"
        );
    }

    use std::time::Instant;

    #[test]
    #[ignore = "timing comparison only"]
    fn compare_newton_vs_quartic_timing() {
        let inputs = base_inputs(0);
        let (i_params, r_params) = base_params(1.0e64, 1.0e18);

        let samples = 100_000;
        let mu_values: Vec<f64> = (0..samples)
            .map(|i| 0.1 + 0.8 * (i as f64 / samples as f64))
            .collect();

        let start = Instant::now();
        let mut newton_sum = 0.0;
        for &mu in &mu_values {
            newton_sum += r_of_theta_internal(mu, &i_params, &r_params, &inputs, true);
        }
        let newton_elapsed = start.elapsed();

        let start = Instant::now();
        let mut quartic_sum = 0.0;
        for &mu in &mu_values {
            quartic_sum += r_of_theta_internal(mu, &i_params, &r_params, &inputs, false);
        }
        let quartic_elapsed = start.elapsed();

        println!("Newton elapsed:  {:?}", newton_elapsed);
        println!("Quartic elapsed: {:?}", quartic_elapsed);
        println!("Newton sum:      {newton_sum}");
        println!("Quartic sum:     {quartic_sum}");
        println!(
            "Speedup: {:.3}x",
            newton_elapsed.as_secs_f64() / quartic_elapsed.as_secs_f64()
        );
    }
}
