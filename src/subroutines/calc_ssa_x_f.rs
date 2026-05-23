use crate::modules::constants::THIRD;
use crate::modules::globals::Globals;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalcSsaXFResult {
    pub x_f: f64,
    pub n_x_out: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CalcSsaXFError {
    SynXArrayTooShort { len: usize, required_min_len: usize },
    SynFArrayTooShort { len: usize, required_min_len: usize },
}

/// Computes F(x) from Rybicki & Lightman (1979), Eq. 6.18. Does this by
/// linear interpolation on a pre-computed array.
///
/// # Arguments
///
/// * `xxx` - value of `nu / nu_c` being used
/// * `first_call` - boolean variable telling us how to seek for the
///   appropriate entry of the `syn_x` array
/// * `n_x_in` - value of `n_x_out` the last time we called this function.
///   Ignored if `first_call` is true
/// * `globals` - runtime global state containing `syn_x`, `syn_f`, and
///   `o_o_sone`
pub fn calc_ssa_x_f(
    xxx: f64,
    first_call: bool,
    n_x_in: usize,
    globals: &Globals,
) -> Result<CalcSsaXFResult, CalcSsaXFError> {
    const N_X_LOW: usize = 1;
    const N_X_HIGH: usize = 118;
    const N_X_HIGH_TAIL: usize = 117;

    if globals.syn_x.len() <= N_X_HIGH {
        return Err(CalcSsaXFError::SynXArrayTooShort {
            len: globals.syn_x.len(),
            required_min_len: N_X_HIGH + 1,
        });
    }

    if globals.syn_f.len() <= N_X_HIGH {
        return Err(CalcSsaXFError::SynFArrayTooShort {
            len: globals.syn_f.len(),
            required_min_len: N_X_HIGH + 1,
        });
    }

    // If xxx is sufficiently low (or sufficiently high), we can calculate x_F differently.
    // If not, we will need to linearly interpolate after identifying the
    // entries in syn_x that bracket xxx
    if xxx <= globals.syn_x[1] {
        Ok(CalcSsaXFResult {
            x_f: (xxx * globals.o_o_sone).powf(THIRD) * globals.syn_f[1],
            n_x_out: 1,
        })
    } else if xxx >= 30.0 {
        Ok(CalcSsaXFResult {
            x_f: 1.0e-99,
            n_x_out: N_X_HIGH_TAIL,
        })
    } else {
        // If this is the first time we've called the subroutine, we need to
        //   search for the correct entries of the syn_x array.  Otherwise, we can
        //   use the value of n_x_in (from the last time we called this
        //   subroutine)
        let n_x = if first_call {
            let mut n_x_lo = N_X_LOW;
            let mut n_x_hi = N_X_HIGH;

            while (n_x_hi - n_x_lo) > 1 {
                // To calculate n_x = (n_x_lo+n_x_hi)/2, use bit shifting to the right
                // instead of a division by two. This *only* works because n_x_lo
                // and n_x_hi are guaranteed to be positive; otherwise the sign bit
                // would cause all sorts of problems
                let n_x = (n_x_lo + n_x_hi) >> 1;

                if globals.syn_x[n_x] > xxx {
                    n_x_hi = n_x;
                } else {
                    n_x_lo = n_x;
                }
            }

            n_x_lo
        } else {
            let mut found = N_X_LOW;

            for n_x in (N_X_LOW..=n_x_in).rev() {
                if globals.syn_x[n_x] <= xxx {
                    found = n_x;
                    break;
                }
            }

            found
        };

        // Interpolate between entries of syn_F to find x_F
        let x_f = (xxx - globals.syn_x[n_x]) / (globals.syn_x[n_x + 1] - globals.syn_x[n_x])
            * (globals.syn_f[n_x + 1] - globals.syn_f[n_x])
            + globals.syn_f[n_x];

        // Set the value of n_x_out so that we're in position the next time we
        // call this function
        Ok(CalcSsaXFResult { x_f, n_x_out: n_x })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn globals_for_test() -> Globals {
        let mut syn_x = vec![0.0; 119];
        let mut syn_f = vec![0.0; 119];

        for i in 0..119 {
            syn_x[i] = i as f64;
            syn_f[i] = 2.0 * i as f64;
        }

        Globals {
            syn_x,
            syn_f,
            o_o_sone: 1.0,
            ..Default::default()
        }
    }

    #[test]
    fn uses_low_x_power_law_branch() {
        let globals = globals_for_test();

        let result = calc_ssa_x_f(0.5, true, 117, &globals).expect("calc_ssa_x_f should succeed");

        let expected = 0.5_f64.powf(THIRD) * globals.syn_f[1];

        assert!((result.x_f - expected).abs() < 1.0e-12);
        assert_eq!(result.n_x_out, 1);
    }

    #[test]
    fn uses_high_x_tail_branch() {
        let globals = globals_for_test();

        let result = calc_ssa_x_f(30.0, true, 117, &globals).expect("calc_ssa_x_f should succeed");

        assert_eq!(result.x_f, 1.0e-99);
        assert_eq!(result.n_x_out, 117);
    }

    #[test]
    fn interpolates_on_first_call() {
        let globals = globals_for_test();

        let result = calc_ssa_x_f(10.5, true, 117, &globals).expect("calc_ssa_x_f should succeed");

        assert!((result.x_f - 21.0).abs() < 1.0e-12);
        assert_eq!(result.n_x_out, 10);
    }

    #[test]
    fn interpolates_using_prior_index_when_not_first_call() {
        let globals = globals_for_test();

        let result = calc_ssa_x_f(10.5, false, 15, &globals).expect("calc_ssa_x_f should succeed");

        assert!((result.x_f - 21.0).abs() < 1.0e-12);
        assert_eq!(result.n_x_out, 10);
    }

    #[test]
    fn returns_error_when_syn_x_is_too_short() {
        let mut globals = globals_for_test();
        globals.syn_x.truncate(118);

        let result = calc_ssa_x_f(10.5, true, 117, &globals);

        assert!(matches!(
            result,
            Err(CalcSsaXFError::SynXArrayTooShort {
                len: 118,
                required_min_len: 119
            })
        ));
    }

    #[test]
    fn returns_error_when_syn_f_is_too_short() {
        let mut globals = globals_for_test();
        globals.syn_f.truncate(118);

        let result = calc_ssa_x_f(10.5, true, 117, &globals);

        assert!(matches!(
            result,
            Err(CalcSsaXFError::SynFArrayTooShort {
                len: 118,
                required_min_len: 119
            })
        ));
    }
}
