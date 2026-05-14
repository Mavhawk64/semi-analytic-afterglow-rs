const MAX_ITRS: usize = 100;
const TOL: f64 = 1.0e-9;

#[derive(Debug, Clone, PartialEq)]
pub enum NewtonError {
    MaxIterations { max_itrs: usize, last_x: f64 },
    ZeroDerivative { x_curr: f64, xpd: f64 },
    NonFiniteValue { x_curr: f64 },
}

/// Uses Newton's method to find the root of a function. Function name is
/// passed to the subroutine. Note that the function has a very specific list
/// of arguments -- any additional inputs must be included from modules!
///
/// # Arguments
///
/// * `func_name` - name of the function we're trying to find the root of
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters to be used in evaluating the function
/// * `x_start` - initial guess for the location of the root
/// * `l_bound` - minimum allowed value for x
/// * `u_bound` - maximum allowed value for x
pub fn newtons_method<F>(
    func_name: F,
    i_params: &[i32; 7],
    r_params: &[f64; 7],
    x_start: f64,
    l_bound: f64,
    u_bound: f64,
) -> Result<f64, NewtonError>
where
    F: Fn(f64, &[i32; 7], &[f64; 7]) -> f64,
{
    let mut x_curr: f64 = x_start;
    let mut xpd: f64 = next_xpd(x_curr, u_bound);

    for _ in 1..=MAX_ITRS {
        // Compute f(x_n), f'(x_n)
        let f_x: f64 = func_name(x_curr, i_params, r_params);
        let f_xpd: f64 = func_name(xpd, i_params, r_params);

        if !f_x.is_finite() || !f_xpd.is_finite() {
            return Err(NewtonError::NonFiniteValue { x_curr });
        }

        let denominator = xpd - x_curr;
        if denominator == 0.0 {
            return Err(NewtonError::ZeroDerivative { x_curr, xpd });
        }

        // Dy/Dx
        let f_prime: f64 = (f_xpd - f_x) / denominator;
        if f_prime == 0.0 || !f_prime.is_finite() {
            return Err(NewtonError::ZeroDerivative { x_curr, xpd });
        }

        // Compute x_(n+1); make sure it doesn't fall outside the bounds set
        // by l_bound and u_bound
        let mut x_next: f64 = x_curr - f_x / f_prime;
        x_next = x_next.max(0.5 * (l_bound + x_curr));
        x_next = x_next.min(0.5 * (x_curr + u_bound));

        // If we can exit, do so...
        if (x_next - x_curr).abs() < TOL * x_curr.abs() {
            return Ok(x_curr);
        }

        // ...otherwise prepare for next time through loop
        x_curr = x_next;
        xpd = next_xpd(x_curr, u_bound);
    }

    // If the loop ended without finishing, tell us.
    Err(NewtonError::MaxIterations {
        max_itrs: MAX_ITRS,
        last_x: x_curr,
    })
}

fn next_xpd(x_curr: f64, u_bound: f64) -> f64 {
    if x_curr == 0.0 {
        // If x_curr is 0, can't just multiply by a factor to get
        // x_curr + delta.
        (1.0e-5 * u_bound).min(1.0e-5)
    } else {
        (1.0001 * x_curr).min(0.5 * (x_curr + u_bound))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Happy path:

    fn quadratic(x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        x * x - 4.0
    }

    #[test]
    fn finds_positive_root_of_simple_quadratic() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let root = newtons_method(quadratic, &i_params, &r_params, 3.0, 0.0, 10.0)
            .expect("Newton's method should find the root");

        assert!((root - 2.0).abs() < 1.0e-6);
    }

    // Sad paths:

    fn constant_function(_x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        1.0
    }

    fn non_finite_function(_x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        f64::NAN
    }

    fn exponential_no_root(x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        x.exp()
    }

    #[test]
    fn returns_zero_derivative_for_constant_function() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = newtons_method(constant_function, &i_params, &r_params, 3.0, 0.0, 10.0);

        assert!(matches!(result, Err(NewtonError::ZeroDerivative { .. })));
    }

    #[test]
    fn returns_non_finite_value_for_nan_function() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = newtons_method(non_finite_function, &i_params, &r_params, 3.0, 0.0, 10.0);

        assert!(matches!(result, Err(NewtonError::NonFiniteValue { .. })));
    }

    #[test]
    fn returns_zero_derivative_when_finite_difference_collapses() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = newtons_method(exponential_no_root, &i_params, &r_params, 1.0, 0.0, 10.0);

        assert!(matches!(result, Err(NewtonError::ZeroDerivative { .. })));
    }

    fn oscillatory_no_root(x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        2.0 + (1.0e6 * x).sin()
    }

    #[test]
    fn returns_max_iterations_for_oscillatory_function_without_root() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = newtons_method(oscillatory_no_root, &i_params, &r_params, 1.0, 0.0, 10.0);

        assert!(matches!(result, Err(NewtonError::MaxIterations { .. })));
    }
}
