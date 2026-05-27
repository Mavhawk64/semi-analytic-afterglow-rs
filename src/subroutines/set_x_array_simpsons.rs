/// Initializes the values of the x_array, which will be used to integrate
/// Eq. (A24). Because we will be using Simpson's rule, these are organized
/// into sets of three points with repeated boundary points.
///
/// # Returns
///
/// * `(num_x, x_array)` - number of Simpson intervals and the corresponding
///   array of `[x_left, x_mid, x_right]` values.
pub fn set_x_array_simpsons() -> (usize, Vec<[f64; 3]>) {
    // Set x_array values
    //
    // NOTE:
    // The original Fortran first assigns a 10-interval grid, then immediately
    // overwrites it with the DEBUGLINE 15-interval grid. This function returns
    // the final effective Fortran values.
    // Alternate lower-resolution Simpson grid preserved from the original
    // Fortran source:
    //
    // let x_array = vec![
    //     [0.0, 0.12, 0.24],
    //     [0.24, 0.36, 0.48],
    //     [0.48, 0.56, 0.64],
    //     [0.64, 0.70, 0.76],
    //     [0.76, 0.80, 0.84],
    //     [0.84, 0.87, 0.90],
    //     [0.90, 0.92, 0.94],
    //     [0.94, 0.95, 0.96],
    //     [0.96, 0.97, 0.98],
    //     [0.98, 0.99, 1.0],
    // ];
    //
    // let num_x = x_array.len();
    let x_array = vec![
        [0.0, 0.08, 0.16],
        [0.16, 0.24, 0.32],
        [0.32, 0.40, 0.48],
        [0.48, 0.54, 0.60],
        [0.60, 0.65, 0.70],
        [0.70, 0.73, 0.76],
        [0.76, 0.78, 0.80],
        [0.80, 0.82, 0.84],
        [0.84, 0.86, 0.88],
        [0.88, 0.90, 0.92],
        [0.92, 0.93, 0.94],
        [0.94, 0.95, 0.96],
        [0.96, 0.97, 0.98],
        [0.98, 0.985, 0.99],
        [0.99, 0.995, 1.0],
    ];

    let num_x = x_array.len();

    (num_x, x_array)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_expected_number_of_intervals() {
        let (num_x, x_array) = set_x_array_simpsons();

        assert_eq!(num_x, 15);
        assert_eq!(x_array.len(), 15);
    }

    #[test]
    fn first_and_last_intervals_match_fortran_effective_values() {
        let (_num_x, x_array) = set_x_array_simpsons();

        assert_eq!(x_array[0], [0.0, 0.08, 0.16]);
        assert_eq!(x_array[14], [0.99, 0.995, 1.0]);
    }

    #[test]
    fn each_interval_is_strictly_increasing() {
        let (_num_x, x_array) = set_x_array_simpsons();

        for [left, mid, right] in x_array {
            assert!(left < mid);
            assert!(mid < right);
        }
    }

    #[test]
    fn adjacent_intervals_share_boundaries() {
        let (_num_x, x_array) = set_x_array_simpsons();

        for window in x_array.windows(2) {
            assert_eq!(window[0][2], window[1][0]);
        }
    }

    #[test]
    fn all_values_are_between_zero_and_one_inclusive() {
        let (_num_x, x_array) = set_x_array_simpsons();

        for interval in x_array {
            for value in interval {
                assert!((0.0..=1.0).contains(&value));
            }
        }
    }
}
