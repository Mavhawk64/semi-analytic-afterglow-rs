const MAX_ITRS: usize = 100;
const TOL: f64 = 1.0e-6;
const CGOLD: f64 = 0.381966;
const EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrentResult {
    pub x_min: f64,
    // Additional information about the optimization process can be included here if desired.
    // Not originally in f90 code.
    pub iterations: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BrentError {
    MaxIterations { max_itrs: usize, last_x: f64 },
    InvalidBracket,
    NonFiniteValue { x: f64, fx: f64 },
}

/// Given a function func_name, and given a bracketing triplet of abscissas
///   ax, bx, and cx (such that bx lies between ax and cx, and f(bx) is less
///   than both f(ax) and f(cx)), this routine isolates the function's minimum
///   to a fractional precision of about tol (set as a parameter below) using
///   Brent's method.  The abscissa of the minimum is returned as x_min (we do
///   not return the function value at xmin, since we don't need it).
/// Note that tol should be at minimum 3e-8 for double-precision calculations.
///   Accuracy is limited to roughly the square root of machine precision, so
///   reducing tol means more function evaluations for marginal (or no) gains.
///
/// Function name is passed to the subroutine, and interface is written below.
///   Note that the function has a very specific list of arguments -- any
///   additional inputs must be included from modules!
///
/// Subroutine modified from brent in Numerical Recipes in Fortran 77,
///   2nd edition.
///
/// # Arguments
///
/// * `func_name` - name of the function we're trying to minimize
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters to be used in evaluating the function
/// * `ax` - x value that brackets the function minimum
/// * `bx` - x value that brackets the function minimum
/// * `cx` - x value that brackets the function minimum
///
/// # Returns
/// * `BrentResult` containing the x value of the minimum and other information
/// # Errors
/// * `BrentError::InvalidBracket` if the initial bracketing conditions are not met
/// * `BrentError::NonFiniteValue` if the function returns a non-finite value at any of the initial points or during the optimization process
/// * `BrentError::MaxIterations` if the maximum number of iterations is reached without convergence
///
/// # Notes on F90 -> Rust
/// * The original code would return a variable `istat` to indicate 1 = success, -1 = failure.
///   In Rust, we use a `Result` type to return either a successful `BrentResult` or a `BrentError`,
///   effectively encoding the success/failure state in the type system.
pub fn brents_minimum<F>(
    func_name: F,
    i_params: &[i32; 7],
    r_params: &[f64; 7],
    ax: f64,
    bx: f64,
    cx: f64,
) -> Result<BrentResult, BrentError>
where
    F: Fn(f64, &[i32; 7], &[f64; 7]) -> f64,
{
    // This is not in the original Fortran code,
    // but this is a good place to debug if our assumptions
    // are incorrect (we don't want that).
    // We'll return a BrentError if any of these conditions aren't met.
    let fa = func_name(ax, i_params, r_params);
    let fb = func_name(bx, i_params, r_params);
    let fc = func_name(cx, i_params, r_params);

    if !fa.is_finite() {
        return Err(BrentError::NonFiniteValue { x: ax, fx: fa });
    }
    if !fb.is_finite() {
        return Err(BrentError::NonFiniteValue { x: bx, fx: fb });
    }
    if !fc.is_finite() {
        return Err(BrentError::NonFiniteValue { x: cx, fx: fc });
    }

    if !((ax < bx && bx < cx) || (cx < bx && bx < ax)) || !(fb <= fa && fb <= fc) {
        return Err(BrentError::InvalidBracket);
    }

    // Initialize values. a and b must be in ascending order,
    // even though ax and cx need not be.
    let mut a = ax.min(cx);
    let mut b = ax.max(cx);

    let mut v = bx;
    let mut w = v;
    let mut x = v;

    // Distance moved on the step before last
    let mut e: f64 = 0.0;
    let mut d: f64 = 0.0; // set a temporary value for d (not done in f90)

    let mut fx = fb;
    let mut fv = fx;
    let mut fw = fx;

    for iteration in 1..=MAX_ITRS {
        let xm = 0.5 * (a + b);
        let tol1 = TOL * x.abs() + EPS;
        let tol2 = 2.0 * tol1;

        // Test for convergence here.
        if (x - xm).abs() <= tol2 - 0.5 * (b - a) {
            return Ok(BrentResult {
                x_min: x,
                iterations: iteration,
            });
        }

        let mut skip_gold = false;
        let mut u: f64;

        // If we're still here, construct a trial parabolic fit.
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;

            q = 2.0 * (q - r);

            if q > 0.0 {
                p = -p;
            }

            q = q.abs();

            let etemp = e;
            e = d;

            // The conditions below determine the acceptability of the parabolic
            // fit. If any is true, the fit fails and the golden section step is
            // needed.
            // Note: The logic is flipped from the f90 because we don't need an else statement.
            if !(p.abs() >= (0.5 * q * etemp).abs() || p <= q * (a - x) || p >= q * (b - x)) {
                // Can skip the golden section step.
                d = p / q;
                u = x + d;

                if (u - a) < tol2 || (b - u) < tol2 {
                    d = tol1.copysign(xm - x);
                }

                skip_gold = true;
            }
        }

        // Golden section step, if needed.
        if !skip_gold {
            if x >= xm {
                e = a - x;
            } else {
                e = b - x;
            }

            d = CGOLD * e;
        }

        // We now have d computed, either from the parabolic fit or from the
        // golden section.
        if d.abs() >= tol1 {
            u = x + d;
        } else {
            u = x + tol1.copysign(d);
        }

        // Perform the one function evaluation for this iteration, and do some
        // housekeeping afterward.
        let fu = func_name(u, i_params, r_params);
        // Same as before... error if bad (not in f90 code, but good for debugging).
        if !fu.is_finite() {
            return Err(BrentError::NonFiniteValue { x: u, fx: fu });
        }

        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }

            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }

            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;

                // NOTE:
                // The Fortran source has `fv = u` here, but `fv` stores a
                // function value, not an x-location. This should be `fv = fu`.
                fv = fu;
            }
        }
    }
    // Otherwise, we failed to converge within the maximum number of iterations.
    Err(BrentError::MaxIterations {
        max_itrs: MAX_ITRS,
        last_x: x,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parabola(x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        (x - 2.0).powi(2)
    }

    fn shifted_parabola(x: f64, _i_params: &[i32; 7], r_params: &[f64; 7]) -> f64 {
        let center = r_params[0];
        let offset = r_params[1];

        (x - center).powi(2) + offset
    }

    fn non_finite_function(_x: f64, _i_params: &[i32; 7], _r_params: &[f64; 7]) -> f64 {
        f64::NAN
    }

    #[test]
    fn finds_minimum_of_simple_parabola() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = brents_minimum(parabola, &i_params, &r_params, 0.0, 2.0, 4.0)
            .expect("Brent should find the minimum");

        assert!((result.x_min - 2.0).abs() < 1.0e-5);
    }

    #[test]
    fn finds_minimum_of_shifted_parabola_using_params() {
        let i_params = [0; 7];

        let mut r_params = [0.0; 7];
        r_params[0] = -1.5;
        r_params[1] = 7.0;

        let result = brents_minimum(shifted_parabola, &i_params, &r_params, -4.0, -1.5, 2.0)
            .expect("Brent should find the shifted minimum");

        assert!((result.x_min + 1.5).abs() < 1.0e-5);
    }

    #[test]
    fn accepts_reversed_ax_cx_order() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = brents_minimum(parabola, &i_params, &r_params, 4.0, 2.0, 0.0)
            .expect("Brent should accept reversed ax/cx");

        assert!((result.x_min - 2.0).abs() < 1.0e-5);
    }

    #[test]
    fn rejects_invalid_bracket_when_bx_not_between_ax_and_cx() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = brents_minimum(parabola, &i_params, &r_params, 0.0, 5.0, 4.0);

        assert!(matches!(result, Err(BrentError::InvalidBracket)));
    }

    #[test]
    fn rejects_invalid_bracket_when_bx_is_not_lower_than_endpoints() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        // f(3) = 1
        // f(0) = 4
        // f(2) = 0
        //
        // bx is NOT lower than both endpoints.
        let result = brents_minimum(parabola, &i_params, &r_params, 0.0, 3.0, 2.0);

        assert!(matches!(result, Err(BrentError::InvalidBracket)));
    }

    #[test]
    fn returns_non_finite_value_error() {
        let i_params = [0; 7];
        let r_params = [0.0; 7];

        let result = brents_minimum(non_finite_function, &i_params, &r_params, 0.0, 1.0, 2.0);

        assert!(matches!(result, Err(BrentError::NonFiniteValue { .. })));
    }
}
