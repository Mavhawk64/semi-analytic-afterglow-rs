use crate::modules::inputs::Inputs;

/// Computes coasting Lorentz factor at a polar angle theta_fl, using the
/// fitting formula of Gottlieb+ 2021 (2021MNRAS.500.3511G)
///
/// Returns Gam(theta), isotropic equivalent energy.
///
/// # Arguments
///
/// * `theta_fl` - angle relative to jet axis
/// * `inputs` - runtime input parameters, replacing Fortran module-level inputs
pub fn gam_of_theta(theta_fl: f64, inputs: &Inputs) -> f64 {
    // Compute coasting (peak) Lorentz factor based on whether we're in jet core
    // or not
    //
    // Commented out version is original from paper, and has discontinuous first
    // derivative at theta_j; this leads to cusp-like behavior. It has been
    // replaced with a version that has continuous first derivative.
    if inputs.do_jet_structure {
        // if theta_fl < theta_j {
        //     Gam_of_theta = eta
        // } else {
        //     Gam_of_theta = eta * (theta_j / theta_fl)^jet_p
        // }

        inputs.eta * (1.0 + (theta_fl / inputs.theta_j).powf(4.0 * inputs.jet_p)).powf(-0.25)
    } else {
        // If there is no jet structure, we're using a spherically-symmetric
        // fireball so the coasting Lorentz factor is the same regardless of
        // angle
        inputs.eta
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_eta_when_no_jet_structure() {
        let inputs = Inputs {
            do_jet_structure: false,
            eta: 300.0,
            ..Default::default()
        };

        let result = gam_of_theta(0.5, &inputs);

        assert_eq!(result, 300.0);
    }

    #[test]
    fn returns_structured_lorentz_factor() {
        let inputs = Inputs {
            do_jet_structure: true,
            eta: 300.0,
            theta_j: 0.1,
            jet_p: 2.0,
            ..Default::default()
        };

        let theta_fl: f64 = 0.2;

        let expected =
            inputs.eta * (1.0 + (theta_fl / inputs.theta_j).powf(4.0 * inputs.jet_p)).powf(-0.25);

        let result = gam_of_theta(theta_fl, &inputs);

        assert!((result - expected).abs() < 1.0e-12 * expected.abs());
    }
}
