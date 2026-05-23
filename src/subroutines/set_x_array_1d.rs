/// Initializes the values of the x_array, which will be used to integrate
/// Eq. (A24). Because we will be using Simpson's rule, these are organized
/// into sets of three points with repeated boundary points.
///
/// # Returns
///
/// * `(num_x, x_array)` - number of x values and the corresponding array
pub fn set_x_array_1d() -> (usize, Vec<f64>) {
    // Set x_array values
    let x_array = vec![
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.70, 0.75, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90,
        0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.993, 0.996, 0.999,
    ];

    let num_x = x_array.len();

    (num_x, x_array)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_expected_number_of_points() {
        let (num_x, x_array) = set_x_array_1d();

        assert_eq!(num_x, 28);
        assert_eq!(x_array.len(), 28);
    }

    #[test]
    fn first_and_last_values_match_fortran() {
        let (_num_x, x_array) = set_x_array_1d();

        assert_eq!(x_array[0], 0.0);
        assert_eq!(x_array[27], 0.999);
    }

    #[test]
    fn x_array_is_strictly_increasing() {
        let (_num_x, x_array) = set_x_array_1d();

        for window in x_array.windows(2) {
            assert!(window[1] > window[0]);
        }
    }

    #[test]
    fn all_x_values_are_between_zero_and_one() {
        let (_num_x, x_array) = set_x_array_1d();

        assert!(x_array.iter().all(|x| *x >= 0.0 && *x < 1.0));
    }

    #[test]
    fn contains_expected_high_resolution_tail() {
        let (_num_x, x_array) = set_x_array_1d();

        assert_eq!(&x_array[24..28], &[0.99, 0.993, 0.996, 0.999]);
    }
}
