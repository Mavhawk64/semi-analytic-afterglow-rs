use crate::subroutines::beschb::{BeschbResult, beschb};

const MAXIT: usize = 10_000;
const XMIN: f64 = 2.0;
const FPMIN: f64 = 1.0e-30;
const PI: f64 = std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BessikResult {
    pub ri: f64,
    pub rk: f64,
    pub rip: f64,
    pub rkp: f64,
}

/// Returns the modified Bessel functions ri = I_nu, rk = K_nu, and their
/// derivatives rip = I'_nu, rkp = K'_nu, for positive x and xnu (= nu) >= 0.
///
/// The relative accuracy is within one or two significant digits of the input
/// argument eps. FPMIN is a parameter set close to the machine's smallest
/// floating point number. All internal arithmetic is in double precision.
///
/// Uses beschb, which uses chebev.
///
/// Subroutine taken almost unmodified from Numerical Recipes in Fortran 77,
/// 2nd edition.
///
/// # Arguments
///
/// * `x` - positive input value
/// * `xnu` - order of the Bessel functions, must be non-negative
/// * `eps` - convergence tolerance
pub fn bessik(x: f64, xnu: f64, eps: f64) -> BessikResult {
    if x <= 0.0 || xnu < 0.0 {
        panic!("bad arguments in bessik: x = {x}, xnu = {xnu}");
    }

    // n is the number of downward recurrences of the I's and upward recurrences
    // of the K's. xmu lies between -1/2 and +1/2
    let nl = (xnu + 0.5).trunc() as usize;
    let xmu = xnu - nl as f64;
    let xmu2 = xmu * xmu;
    let xi = 1.0 / x;
    let xi2 = 2.0 * xi;
    let mut h = xnu * xi;

    // Evaluate CF1 by modified Lentz's method (section 5.2 in Numerical Recipes
    // in Fortran 77, 2nd ed.)
    if h < FPMIN {
        h = FPMIN;
    }

    let mut b = xi2 * xnu;
    let mut d = 0.0;
    let mut c = h;

    let mut converged_cf1 = false;
    for _ in 1..=MAXIT {
        b += xi2;
        d = 1.0 / (b + d);
        c = b + 1.0 / c;

        let del = c * d;
        h *= del;

        if (del - 1.0).abs() < eps {
            converged_cf1 = true;
            break;
        }
    }

    if !converged_cf1 {
        panic!("x too large in bessik; try asymptotic expansion");
    }

    // Initialize I_nu and I'_nu for downward recurrence, and store values for
    // later rescaling
    let mut ril = FPMIN;
    let mut ripl = h * ril;
    let ril1 = ril;
    let rip1 = ripl;
    let mut fact = xnu * xi;

    for _l in (1..=nl).rev() {
        let ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }

    // We now have unnormalized I_mu and I'_mu. Use the series
    let f = ripl / ril;

    let mut rkmu: f64;
    let mut rk1: f64;

    if x < XMIN {
        let x2 = 0.5 * x;
        let pimu = PI * xmu;

        let fact = if pimu.abs() < eps {
            1.0
        } else {
            pimu / pimu.sin()
        };

        let d_log = -x2.ln();
        let e_arg = xmu * d_log;

        let fact2 = if e_arg.abs() < eps {
            1.0
        } else {
            e_arg.sinh() / e_arg
        };

        // Chebyshev evaluation of Gamma_1 and Gamma_2
        let BeschbResult {
            gam1,
            gam2,
            gampl,
            gammi,
        } = beschb(xmu);

        let mut ff = fact * (gam1 * e_arg.cosh() + gam2 * fact2 * d_log); // f0
        let mut sum = ff;
        let e_exp = e_arg.exp();

        let mut p = 0.5 * e_exp / gampl; // p0
        let mut q = 0.5 / (e_exp * gammi); // q0

        let mut c_series = 1.0;
        let d_series = x2 * x2;
        let mut sum1 = p;

        let mut converged_series = false;
        for i in 1..=MAXIT {
            let i_f = i as f64;

            ff = (i_f * ff + p + q) / (i_f * i_f - xmu2);
            c_series = c_series * d_series / i_f;
            p /= i_f - xmu;
            q /= i_f + xmu;

            let del = c_series * ff;
            sum += del;

            let del1 = c_series * (p - i_f * ff);
            sum1 += del1;

            if del.abs() < sum.abs() * eps {
                converged_series = true;
                break;
            }
        }

        if !converged_series {
            // TODO: f90 had "bessk" instead of "bessik" here -- check that this is correct
            panic!("bessik series failed to converge");
        }

        rkmu = sum;
        rk1 = sum1 * xi2;
    } else {
        // Otherwise, evaluate CF2 by Steed's algorithm (also section 5.2 of
        // Numerical Recipes in Fortran 77, 2nd ed.). This is OK because there
        // can be no zero denominators
        b = 2.0 * (1.0 + x);
        d = 1.0 / b;

        let mut delh = d;
        h = delh;

        let mut q1 = 0.0; // Initializations for recurrence 6.7.35
        let mut q2 = 1.0;

        let a1 = 0.25 - xmu2;
        c = a1;

        let mut q = c; // First term in equation 6.7.34
        let mut a = -a1;
        let mut s = 1.0 + q * delh;

        let mut converged_cf2 = false;
        for i in 2..=MAXIT {
            let i_f = i as f64;

            a -= 2.0 * (i_f - 1.0);
            c = -a * c / i_f;

            let qnew = (q1 - b * q2) / a;
            q1 = q2;
            q2 = qnew;

            q += c * qnew;
            b += 2.0;

            d = 1.0 / (b + a * d);
            delh *= b * d - 1.0;
            h += delh;

            let dels = q * delh;
            s += dels;

            // Only need to test convergence of sum since CF2 itself converges
            // more quickly
            if (dels / s).abs() < eps {
                converged_cf2 = true;
                break;
            }
        }

        if !converged_cf2 {
            panic!("bessik: failure to converge in cf2");
        }

        h *= a1;

        if x > 350.0 {
            // TODO: Check this -- f90 only sets rk and returns early.
            return BessikResult {
                ri: 0.0,
                rk: 1.0e-150,
                rip: 0.0,
                rkp: 0.0,
            };
        }

        rkmu = (PI / (2.0 * x)).sqrt() * (-x).exp() / s; // Can omit the factor exp(-x) to scale
        rk1 = rkmu * (xmu + x + 0.5 - h) * xi; // all the returned functions by exp(x) for x >= XMIN
    }

    let rkmup = xmu * xi * rkmu - rk1;
    let rimu = xi / (f * rkmu - rkmup); // Get I_mu from the Wronskian

    let ri = (rimu * ril1) / ril; // Scale original I_nu and I'_nu
    let rip = (rimu * rip1) / ril;

    // Now use the upward recurrence of K_nu
    for i in 1..=nl {
        let i_f = i as f64;

        let rktemp = (xmu + i_f) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }

    let rk = rkmu;
    let rkp = xnu * xi * rkmu - rk1;

    BessikResult { ri, rk, rip, rkp }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, rel_tol: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        assert!(
            (actual - expected).abs() / scale < rel_tol,
            "actual = {actual}, expected = {expected}, rel_err = {}",
            (actual - expected).abs() / scale
        );
    }

    #[test]
    fn returns_finite_values_for_small_x_branch() {
        let result = bessik(1.0, 0.5, 1.0e-12);

        assert!(result.ri.is_finite());
        assert!(result.rk.is_finite());
        assert!(result.rip.is_finite());
        assert!(result.rkp.is_finite());
    }

    #[test]
    fn returns_finite_values_for_large_x_branch() {
        let result = bessik(3.0, 0.5, 1.0e-12);

        assert!(result.ri.is_finite());
        assert!(result.rk.is_finite());
        assert!(result.rip.is_finite());
        assert!(result.rkp.is_finite());
    }

    #[test]
    fn half_order_values_match_closed_form_at_x_one() {
        let x: f64 = 1.0;
        let result = bessik(x, 0.5, 1.0e-12);

        let expected_i = (2.0 / (std::f64::consts::PI * x)).sqrt() * x.sinh();
        let expected_k = (std::f64::consts::PI / (2.0 * x)).sqrt() * (-x).exp();

        assert_close(result.ri, expected_i, 1.0e-10);
        assert_close(result.rk, expected_k, 1.0e-10);
    }

    #[test]
    fn wronskian_identity_holds() {
        let x: f64 = 1.0;
        let xnu: f64 = 0.5;
        let result = bessik(x, xnu, 1.0e-12);

        // For modified Bessel functions:
        // I_nu(x) K'_nu(x) - I'_nu(x) K_nu(x) = -1/x
        let wronskian = result.ri * result.rkp - result.rip * result.rk;

        assert_close(wronskian, -1.0 / x, 1.0e-9);
    }

    #[test]
    fn returns_early_for_very_large_x() {
        let result = bessik(351.0, 0.5, 1.0e-12);

        assert_eq!(result.rk, 1.0e-150);
        assert_eq!(result.ri, 0.0);
        assert_eq!(result.rip, 0.0);
        assert_eq!(result.rkp, 0.0);
    }

    #[test]
    #[should_panic(expected = "bad arguments in bessik")]
    fn panics_for_nonpositive_x() {
        bessik(0.0, 0.5, 1.0e-12);
    }

    #[test]
    #[should_panic(expected = "bad arguments in bessik")]
    fn panics_for_negative_order() {
        bessik(1.0, -0.5, 1.0e-12);
    }
}
