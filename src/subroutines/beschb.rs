use crate::subroutines::chebev_subr::chebev_subr;

/// Output values from `beschb`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeschbResult {
    pub gam1: f64,
    pub gam2: f64,
    pub gampl: f64,
    pub gammi: f64,
}

/// Evaluates Gamma_1 and Gamma_2 by Chebyshev expansion for abs(x) <= 1/2.
/// Also returns 1/Gamma(1+x) and 1/Gamma(1-x).
///
/// Uses chebev_subr.
///
/// Subroutine taken almost unmodified from Numerical Recipes in Fortran 77,
/// 2nd edition.
///
/// # Arguments
///
/// * `x` - input value, expected to satisfy abs(x) <= 1/2
pub fn beschb(x: f64) -> BeschbResult {
    const NUSE1: usize = 7;
    const NUSE2: usize = 8;

    let c1 = [
        -1.142022680371172e0,
        6.516511267076e-3,
        3.08709017308e-4,
        -3.470626964e-6,
        6.943764e-9,
        3.6780e-11,
        -1.36e-13,
    ];

    let c2 = [
        1.843740587300906e0,
        -0.076852840844786e0,
        1.271927136655e-3,
        -4.971736704e-6,
        -3.3126120e-8,
        2.42310e-10,
        -1.70e-13,
        -1.0e-15,
    ];

    // Multiply x by 2 to convert range to -1 to 1, and then apply the
    // transformation for evaluating even Chebyshev series
    let xx = 8.0 * x * x - 1.0;

    let gam1 = chebev_subr(-1.0, 1.0, &c1[..NUSE1], xx);
    let gam2 = chebev_subr(-1.0, 1.0, &c2[..NUSE2], xx);

    let gampl = gam2 - x * gam1;
    let gammi = gam2 + x * gam1;

    BeschbResult {
        gam1,
        gam2,
        gampl,
        gammi,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_finite_values_at_zero() {
        let result = beschb(0.0);

        assert!(result.gam1.is_finite());
        assert!(result.gam2.is_finite());
        assert!(result.gampl.is_finite());
        assert!(result.gammi.is_finite());
    }

    #[test]
    fn gampl_and_gammi_match_definitions() {
        let x = 0.25;
        let result = beschb(x);

        let expected_gampl = result.gam2 - x * result.gam1;
        let expected_gammi = result.gam2 + x * result.gam1;

        assert!((result.gampl - expected_gampl).abs() < 1.0e-14);
        assert!((result.gammi - expected_gammi).abs() < 1.0e-14);
    }

    #[test]
    fn gampl_and_gammi_are_equal_when_x_is_zero() {
        let result = beschb(0.0);

        assert!((result.gampl - result.gammi).abs() < 1.0e-14);
    }

    #[test]
    fn values_are_finite_at_valid_range_edges() {
        for x in [-0.5, 0.5] {
            let result = beschb(x);

            assert!(result.gam1.is_finite(), "gam1 was not finite at x={x}");
            assert!(result.gam2.is_finite(), "gam2 was not finite at x={x}");
            assert!(result.gampl.is_finite(), "gampl was not finite at x={x}");
            assert!(result.gammi.is_finite(), "gammi was not finite at x={x}");
        }
    }

    #[test]
    fn destructuring_result_works_cleanly() {
        let BeschbResult {
            gam1,
            gam2,
            gampl,
            gammi,
        } = beschb(0.1);

        assert!(gam1.is_finite());
        assert!(gam2.is_finite());
        assert!(gampl.is_finite());
        assert!(gammi.is_finite());
    }
}
