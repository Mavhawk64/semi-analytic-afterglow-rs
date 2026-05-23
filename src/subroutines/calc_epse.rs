use crate::modules::inputs::Inputs;

/// Determine epsilon_e at the shock based on the provided shock speed.
///
/// # Arguments
///
/// * `gamma` - the Lorentz factor at which we wish to find epsilon_e
/// * `inputs` - runtime input parameters
pub fn calc_epse(_gamma: f64, inputs: &Inputs) -> f64 {
    // Right now, epsilon_e is assumed to be a constant. But since we have this
    // subroutine we can easily change that later!
    inputs.epse0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_epse0_from_inputs() {
        let inputs = Inputs {
            epse0: 0.15,
            ..Default::default()
        };

        let result = calc_epse(100.0, &inputs);

        assert_eq!(result, 0.15);
    }

    #[test]
    fn currently_ignores_gamma() {
        let inputs = Inputs {
            epse0: 0.15,
            ..Default::default()
        };

        let low_gamma = calc_epse(10.0, &inputs);
        let high_gamma = calc_epse(1.0e4, &inputs);

        assert_eq!(low_gamma, high_gamma);
    }
}
