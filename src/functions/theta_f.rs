/// Given two sides of a spherical triangle (theta_obs and acos(mu)), and the
/// included angle between them (psi), compute the third angle of the
/// triangle (theta_f).
///
/// NOTE: this is theta_fl, psi is phi_fl in the paper.
///
/// Returns theta_f angle between jet axis and location of fluid parcel.
///
/// # Arguments
///
/// * `theta_obs` - observer angle relative to jet axis
/// * `psi` - angle between theta_obs and mu sides of spherical triangle
/// * `mu` - cosine of angle between observer angle and location of fluid parcel
pub fn theta_f(theta_obs: f64, psi: f64, mu: f64) -> f64 {
    // Compute the output value
    (mu * theta_obs.cos() + (1.0 - mu.powi(2)).sqrt() * theta_obs.sin() * psi.cos()).acos()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_theta_obs_when_mu_is_one() {
        let theta_obs = 0.5;
        let psi = 1.2;
        let mu = 1.0;

        let result = theta_f(theta_obs, psi, mu);

        // When mu = 1:
        //
        // sqrt(1 - mu^2) = 0
        //
        // theta_f = acos(cos(theta_obs)) = theta_obs
        assert!((result - theta_obs).abs() < 1.0e-12);
    }
}
