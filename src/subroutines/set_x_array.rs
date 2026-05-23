/// Initializes the values of the x_array, which will be used to integrate
/// Eq. (A24). Because we will be using Simpson's rule, these are organized
/// into sets of three points with repeated boundary points.
///
/// # Returns
///
/// * `(num_x, x_array)` - number of Simpson intervals and the corresponding
///   array of `[x_left, x_mid, x_right]` values.
pub fn set_x_array() -> (usize, Vec<[f64; 3]>) {
    // Set x_array values
    let x_array = vec![
        [0.0, 0.2, 0.4],
        [0.4, 0.5, 0.6],
        [0.6, 0.7, 0.8],
        [0.8, 0.83, 0.86],
        [0.86, 0.89, 0.92],
        [0.92, 0.94, 0.96],
        [0.96, 0.98, 1.0],
    ];

    let num_x = x_array.len();

    (num_x, x_array)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_expected_number_of_intervals() {
        let (num_x, x_array) = set_x_array();

        assert_eq!(num_x, 7);
        assert_eq!(x_array.len(), 7);
    }

    #[test]
    fn first_and_last_intervals_match_fortran() {
        let (_num_x, x_array) = set_x_array();

        assert_eq!(x_array[0], [0.0, 0.2, 0.4]);
        assert_eq!(x_array[6], [0.96, 0.98, 1.0]);
    }

    #[test]
    fn each_interval_is_strictly_increasing() {
        let (_num_x, x_array) = set_x_array();

        for [left, mid, right] in x_array {
            assert!(left < mid);
            assert!(mid < right);
        }
    }

    #[test]
    fn adjacent_intervals_share_boundaries() {
        let (_num_x, x_array) = set_x_array();

        for window in x_array.windows(2) {
            assert_eq!(window[0][2], window[1][0]);
        }
    }

    #[test]
    fn all_values_are_between_zero_and_one_inclusive() {
        let (_num_x, x_array) = set_x_array();

        for interval in x_array {
            for value in interval {
                assert!((0.0..=1.0).contains(&value));
            }
        }
    }
}
