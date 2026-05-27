/// Sets number of steps and their sizes, to be used in subroutine bw_euler_stepper.
///
/// # Returns
///
/// * `(num_steps, stepsize_array)` - number of Euler steps and array listing
///   the size of each step.
pub fn set_em_step_sizes() -> (usize, Vec<f64>) {
    // Initialize num_steps
    let mut num_steps = 0usize;

    // Set list of step sizes
    let sizes_list = [
        1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 2.0e-2, 1.0e-2,
        3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6,
    ];

    // MAV TODO: Adaptive Mesh Refinement?
    // Set number of steps for each step size. First array is lower-res and
    // faster, as it takes 42% fewer steps than the second array does
    let isteps_list = [
        4usize, 2, // Step size 1e-6,  3e-6
        3, 2, // Step size 1e-5,  3e-5
        3, 2, // Step size 1e-4,  3e-4
        3, 2,  // Step size 1e-3,  3e-3
        5,  // Step size 1e-2
        44, // Step size 2e-2
        5,  // Step size 1e-2
        2, 3, // Step size 3e-3,  1e-3
        2, 3, // Step size 3e-4,  1e-4
        2, 3, // Step size 3e-5,  1e-5
        2, 3, // Step size 3e-6,  1e-6
    ];

    // Lower-resolution alternative preserved from Fortran:
    //
    // let isteps_list = [
    //     9usize, 7, // Step size 1e-6,  3e-6
    //     6, 7, // Step size 1e-5,  3e-5
    //     6, 7, // Step size 1e-4,  3e-4
    //     6, 7, // Step size 1e-3,  3e-3
    //     6, // Step size 1e-2
    //     41, // Step size 2e-2
    //     6, // Step size 1e-2
    //     7, 6, // Step size 3e-3,  1e-3
    //     7, 6, // Step size 3e-4,  1e-4
    //     7, 6, // Step size 3e-5,  1e-5
    //     7, 8, // Step size 3e-6,  1e-6
    // ];

    let total_steps: usize = isteps_list.iter().sum();
    let mut stepsize_array = vec![0.0; total_steps];

    // Loop over sizes, then steps at each size to fill the array
    for (i_size, &step_size) in sizes_list.iter().enumerate() {
        for _j_step in 0..isteps_list[i_size] {
            stepsize_array[num_steps] = step_size;
            num_steps += 1;
        }
    }

    // Check output to make sure it's sensible
    println!("Number of steps to use in Euler method: {}", num_steps + 1);
    println!(
        "Range covered by Euler steps: {:14.7e}",
        stepsize_array.iter().take(num_steps).sum::<f64>() + stepsize_array[0]
    );

    (num_steps, stepsize_array)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_expected_number_of_steps() {
        let (num_steps, stepsize_array) = set_em_step_sizes();

        assert_eq!(num_steps, 95);
        assert_eq!(stepsize_array.len(), 95);
    }

    #[test]
    fn first_and_last_steps_are_one_e_minus_six() {
        let (_num_steps, stepsize_array) = set_em_step_sizes();

        assert_eq!(stepsize_array[0], 1.0e-6);
        assert_eq!(stepsize_array[stepsize_array.len() - 1], 1.0e-6);
    }

    #[test]
    fn contains_expected_large_middle_block() {
        let (_num_steps, stepsize_array) = set_em_step_sizes();

        let count_two_e_minus_two = stepsize_array
            .iter()
            .filter(|value| **value == 2.0e-2)
            .count();

        assert_eq!(count_two_e_minus_two, 44);
    }

    #[test]
    fn range_covered_matches_fortran_diagnostic() {
        let (num_steps, stepsize_array) = set_em_step_sizes();

        let range_covered = stepsize_array.iter().take(num_steps).sum::<f64>() + stepsize_array[0];

        assert!((range_covered - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn all_steps_are_positive() {
        let (_num_steps, stepsize_array) = set_em_step_sizes();

        assert!(stepsize_array.iter().all(|value| *value > 0.0));
    }
}
