use crate::subroutines::bessik::bessik;

/// Pre-calculates F(x) in Rybicki & Lightman's Eq. 6.31c. The arrays syn_F
/// and syn_x will be used in computing synchrotron emission.
///
/// # Arguments
///
/// * `n_ph` - number of photon bins
///
/// # Returns
///
/// * `(syn_f, syn_x)` - arrays containing F(x) values and associated x values.
///   Both arrays are length `n_ph + 1`, matching Fortran `dimension(0:n_ph)`.
pub fn calc_syn_fx(n_ph: usize) -> (Vec<f64>, Vec<f64>) {
    // Local variables

    // Set xnu for use in bessik
    let xnu = 5.0 / 3.0;

    // Initialize syn_x and syn_F
    let mut syn_x = vec![0.0; n_ph + 1];
    let mut syn_f = vec![0.0; n_ph + 1];

    //-------------------------------------------------------------------------
    // Fill syn_x by hand
    //-------------------------------------------------------------------------

    let syn_x_values = [
        0.0, 1.0e-6, 2.0e-6, 3.0e-6, 5.0e-6, 1.0e-5, 2.0e-5, 3.0e-5, 5.0e-5, 1.0e-4, 2.0e-4,
        3.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 3.0e-2, 4.0e-2, 5.0e-2,
        6.0e-2, 7.0e-2, 8.0e-2, 9.0e-2, 1.0e-1, 1.1e-1, 1.2e-1, 1.3e-1, 1.4e-1, 1.5e-1, 1.6e-1,
        1.7e-1, 1.8e-1, 1.9e-1, 2.0e-1, 2.1e-1, 2.2e-1, 2.3e-1, 2.4e-1, 2.5e-1, 2.6e-1, 2.7e-1,
        2.8e-1, 2.9e-1, 3.0e-1, 3.1e-1, 3.2e-1, 3.3e-1, 3.4e-1, 3.5e-1, 3.6e-1, 3.7e-1, 3.8e-1,
        3.9e-1, 4.0e-1, 4.1e-1, 4.2e-1, 4.3e-1, 4.4e-1, 4.5e-1, 4.6e-1, 4.7e-1, 4.8e-1, 4.9e-1,
        5.0e-1, 6.0e-1, 7.0e-1, 8.0e-1, 9.0e-1, 1.0e0, 1.1e0, 1.2e0, 1.3e0, 1.4e0, 1.5e0, 1.6e0,
        1.7e0, 1.8e0, 1.9e0, 2.0e0, 2.1e0, 2.2e0, 2.3e0, 2.4e0, 2.5e0, 2.6e0, 2.7e0, 2.8e0, 2.9e0,
        3.0e0, 4.0e0, 5.0e0, 6.0e0, 7.0e0, 8.0e0, 9.0e0, 1.0e1, 1.1e1, 1.2e1, 1.3e1, 1.4e1, 1.5e1,
        1.6e1, 1.7e1, 1.8e1, 1.9e1, 2.0e1, 2.1e1, 2.2e1, 2.3e1, 2.4e1, 2.5e1, 2.6e1, 2.7e1, 2.8e1,
        2.9e1, 3.0e1,
    ];

    let copy_len = syn_x_values.len().min(syn_x.len());

    syn_x[..copy_len].copy_from_slice(&syn_x_values[..copy_len]);

    //-------------------------------------------------------------------------
    // Now calculate F(x) at each value of syn_x
    //-------------------------------------------------------------------------

    let last_integrated_index = 117.min(n_ph);

    for j in 1..=last_integrated_index {
        let xxx_max: f64 = 30.0;
        let xxx: f64 = syn_x[j];
        let xxx_log = xxx.log10();

        let n_xxx = 200usize;
        let del_xxx = (xxx_max.log10() - xxx_log) / n_xxx as f64;
        let xxx_fac = 10.0_f64.powf(del_xxx);

        let mut sum_k = 0.0;
        let mut xx1 = xxx;
        let mut xx = xx1 / xxx_fac.sqrt(); // Both of these are initialized so they
        let mut xx2 = xxx; // have the right values in the loop

        for _i in 1..=n_xxx {
            xx2 *= xxx_fac; // = 10.0**(xxx_log + (i)*del_xxx)
            xx *= xxx_fac; // = sqrt(xx1*xx2)
            let del = xx2 - xx1;

            // Modified Bessel function of fractional order
            let bessik_result = bessik(xx, xnu, 1.0e-8);

            let rk = bessik_result.rk;

            // Here is the factor accounting for isotropy:
            sum_k += (1.0 - (xxx / xx).powi(2)).sqrt() * rk * del;

            // Use this if not accounting for isotropy:
            // sum_k += rk * del;

            xx1 = xx2;
        }

        syn_f[j] = xxx * sum_k;
    }

    if n_ph >= 118 {
        syn_f[118] = 0.0; // Set this one by hand at the end
    }

    //-------------------------------------------------------------------------
    // F(x) calculated
    //-------------------------------------------------------------------------

    (syn_f, syn_x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_arrays_with_n_ph_plus_one_entries() {
        let n_ph = 118;

        let (syn_f, syn_x) = calc_syn_fx(n_ph);

        assert_eq!(syn_f.len(), n_ph + 1);
        assert_eq!(syn_x.len(), n_ph + 1);
    }

    #[test]
    fn fills_expected_syn_x_values() {
        let (_syn_f, syn_x) = calc_syn_fx(118);

        assert_eq!(syn_x[0], 0.0);
        assert_eq!(syn_x[1], 1.0e-6);
        assert_eq!(syn_x[2], 2.0e-6);
        assert_eq!(syn_x[13], 1.0e-3);
        assert_eq!(syn_x[71], 1.0e0);
        assert_eq!(syn_x[72], 1.1e0);
        assert_eq!(syn_x[118], 3.0e1);
    }

    #[test]
    fn syn_x_is_monotonic_non_decreasing() {
        let (_syn_f, syn_x) = calc_syn_fx(118);

        for window in syn_x.windows(2) {
            assert!(window[1] >= window[0]);
        }
    }

    #[test]
    fn syn_f_values_are_finite() {
        let (syn_f, _syn_x) = calc_syn_fx(118);

        assert!(syn_f.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn computed_syn_f_values_are_nonnegative() {
        let (syn_f, _syn_x) = calc_syn_fx(118);

        assert!(syn_f.iter().all(|value| *value >= 0.0));
    }

    #[test]
    fn final_syn_f_entry_is_set_by_hand() {
        let (syn_f, _syn_x) = calc_syn_fx(118);

        assert_eq!(syn_f[118], 0.0);
    }
}
