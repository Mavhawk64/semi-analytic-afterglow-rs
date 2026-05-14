use crate::modules::inputs::Inputs;

/// Computes isotropic equivalent energy at a polar angle theta_fl, using the
/// fitting formula of Gottlieb+ 2021 (2021MNRAS.500.3511G)
///
/// Returns E(theta), isotropic equivalent energy.
/// # Arguments
///
/// * `theta_fl` - polar angle of fluid element
/// * `inputs` - reference to the Inputs struct containing various input parameters
pub fn e_of_theta(theta_fl: f64, inputs: &Inputs) -> f64 {
    // Compute energy based on whether we're in core, jet-cocoon interface, or
    // cocoon.
    //
    // Commented out version is original from paper, and has discontinuous first
    // derivative at theta_j; this leads to cusp-like behavior of R_perp. It
    // has been replaced with a version that has continuous first derivative.
    if inputs.do_jet_structure {
        // !if( theta_fl .lt. theta_j ) then
        // !  E_of_theta = E_iso
        // !elseif( theta_fl .lt. theta_c ) then
        // !  E_of_theta = E_iso * (theta_fl / theta_j)**(-1.d0*jet_delta)
        // !else
        // !  E_of_theta = E_iso * (theta_j / theta_c)**jet_delta                 &!&
        // !              * exp( -jet_f * (theta_fl - theta_c) )
        // !endif
        if theta_fl < inputs.theta_c {
            inputs.e_iso
                * (1.0 + (theta_fl / inputs.theta_j).powf(4.0 * inputs.jet_delta)).powf(-0.25)
        } else {
            inputs.e_iso
                * (inputs.theta_j / inputs.theta_c).powf(inputs.jet_delta)
                * (-inputs.jet_f * (theta_fl - inputs.theta_c)).exp()
        }
    } else {
        // If there is no jet structure, we're using a spherically-symmetric
        // fireball so the energy is the same regardless of angle
        inputs.e_iso
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_e_iso_when_no_jet_structure() {
        let inputs = Inputs {
            do_jet_structure: false,
            e_iso: 1.0e52,
            ..Default::default()
        };

        let result = e_of_theta(0.5, &inputs);

        assert_eq!(result, 1.0e52);
    }

    #[test]
    fn returns_core_energy_inside_theta_c() {
        let inputs = Inputs {
            do_jet_structure: true,
            e_iso: 1.0e52,
            theta_j: 0.1,
            theta_c: 0.3,
            jet_delta: 2.0,
            ..Default::default()
        };

        let theta_fl = 0.05;

        let expected = inputs.e_iso
            * (1.0 + (theta_fl / inputs.theta_j).powf(4.0 * inputs.jet_delta)).powf(-0.25);

        let result = e_of_theta(theta_fl, &inputs);

        assert!((result - expected).abs() < 1.0e-10 * expected.abs());
    }

    #[test]
    fn returns_cocoon_energy_outside_theta_c() {
        let inputs = Inputs {
            do_jet_structure: true,
            e_iso: 1.0e52,
            theta_j: 0.1,
            theta_c: 0.3,
            jet_delta: 2.0,
            jet_f: 5.0,
            ..Default::default()
        };

        let theta_fl = 0.5;

        let expected = inputs.e_iso
            * (inputs.theta_j / inputs.theta_c).powf(inputs.jet_delta)
            * (-inputs.jet_f * (theta_fl - inputs.theta_c)).exp();

        let result = e_of_theta(theta_fl, &inputs);

        assert!((result - expected).abs() < 1.0e-10 * expected.abs());
    }
}
