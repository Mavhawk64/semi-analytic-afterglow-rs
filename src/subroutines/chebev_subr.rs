/// Chebyshev evaluation. All arguments are input. `c[0..m]` is an array of
/// Chebyshev coefficients, the first `m` elements of `c`, which is the output
/// of another subroutine (`chebft`, which must have previously been called
/// with the same value of `a` and `b`).
///
/// The Chebyshev polynomial is evaluated at a point `y` defined below, and the
/// result is returned.
///
/// Subroutine taken from Numerical Recipes in Fortran 77, 2nd edition.
///
/// # Arguments
///
/// * `a` - lower bound of the Chebyshev interval
/// * `b` - upper bound of the Chebyshev interval
/// * `c` - Chebyshev coefficients
/// * `x` - point to evaluate
pub fn chebev_subr(a: f64, b: f64, c: &[f64], x: f64) -> f64 {
    if (x - a) * (x - b) > 0.0 {
        panic!("x not in range in chebev");
    }

    let m = c.len();

    let mut d = 0.0;
    let mut dd = 0.0;

    let y = (2.0 * x - a - b) / (b - a);
    let y2 = 2.0 * y;

    // Fortran:
    //   do j = m, 2, -1
    //
    // Rust:
    //   j_one_based goes m, m-1, ..., 2
    for j_one_based in (2..=m).rev() {
        // store d in sv before updating d, so that we can use the old value of d
        let sv = d;
        d = y2 * d - dd + c[j_one_based - 1];
        dd = sv;
    }

    y * d - dd + 0.5 * c[0]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates_constant_series() {
        // With only c0 nonzero, chebev returns 0.5*c0.
        let c = [4.0];

        let result = chebev_subr(-1.0, 1.0, &c, 0.25);

        assert_eq!(result, 2.0);
    }

    #[test]
    fn evaluates_linear_chebyshev_series() {
        // f(y) = T1(y) = y
        // Numerical Recipes convention:
        // result = 0.5*c0 + c1*T1(y)
        let c = [0.0, 1.0];

        let result = chebev_subr(-1.0, 1.0, &c, 0.25);

        assert!((result - 0.25).abs() < 1.0e-14);
    }

    #[test]
    fn evaluates_quadratic_chebyshev_series() {
        // T2(y) = 2y^2 - 1
        let c = [0.0, 0.0, 1.0];

        let x: f64 = 0.25;
        let expected = 2.0 * x.powi(2) - 1.0;

        let result = chebev_subr(-1.0, 1.0, &c, x);

        assert!((result - expected).abs() < 1.0e-14);
    }

    #[test]
    #[should_panic(expected = "x not in range in chebev")]
    fn panics_when_x_is_outside_range() {
        let c = [1.0];

        chebev_subr(-1.0, 1.0, &c, 2.0);
    }
}
