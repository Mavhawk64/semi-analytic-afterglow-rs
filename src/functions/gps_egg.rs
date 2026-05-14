use crate::modules::inputs::Inputs;

/// Start with an equation for chi and Eq. (A16) of GranotSari2002:
///
/// chi = 1 + 2*(4-k)*Gam_emis^2*(1 - r/y)            [1]
/// mu  = 1 - [1-chi*y^(4-k)]/[2(4-k)*Gam_axis^2*y]   [2]
///
/// r is the radius of the emission site in units of R_LOS, y is the shock
/// radius at time of emission in units of R_LOS, and mu is the fraction of
/// y lying along the line of sight.
///
/// Note that Eq. [1] comes from BM76's Eq. (27), with ct replaced by their
/// Eq. (26), and truncated at order Gam^-2.
///
/// Equation [2] can be rewritten
///
/// 0 = 1 - mu - [1 - chi*y^(4-k)]/[2*(4-k)*Gam_axis^2*y]   [3]
///
/// Multiply by (2*(4-k) Gam_axis^2 y) to eliminate a division to arrive at
///
/// 0 = (2*(4-k) Gam_axis^2 y)*(1 - mu) + chi*y^(4-k) - 1   [4]
///
/// Now multiply Eq. [1] by y (again to eliminate a division) and we get
///
/// chi*y = y + 2*(4-k)*Gam_emis^2*(y - r)                  [5]
///
/// Assuming that we know Gam_axis and mu, Eqs. [4] and [5] are a system of
/// two equations for the two unknowns chi and y.
///
/// The quantities r, y, and mu are assumed to lie between -1 and +1 (with r
/// and y being further restricted to positive values only).
///
/// Output argument:
///
/// - `gps_egg`: the difference between the right and left sides of
///   Equation [4] above -- obviously we want that to be 0.
///
/// # Arguments
///
/// * `y_input` - current value being tried for y
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used in evaluating the function
/// * `r_params` - array (length 7) of up to 7 floating point parameters to be used in evaluating the function
/// * `inputs` - runtime input parameters, replacing the Fortran module-level `k_cbm`
pub fn gps_egg(y_input: f64, i_params: &[i32; 7], r_params: &[f64; 7], inputs: &Inputs) -> f64 {
    let k_cbm = inputs.k_cbm as f64;

    // `i_params` elements:
    // 1. if 1, inside the egg, if 2, along the shell given mu, if 3 along the shell given `R_perp`
    // 2. value of `isign` (for when along shell given `R_perp`)
    // 5. if 1, assume on-axis. else off-axis
    //
    // `r_params` elements:
    // 1. `r` or `mu` based on `i_params(1)`
    // 2. `mu` or `Gam_axis` based on `i_params(1)`
    // 3. `Gam_axis`
    // 4. `E_ratio`

    let offaxis = i_params[4] != 1;

    let mu: f64;
    let gam_axis: f64;
    let chiy: f64;
    let mut e_ratio: f64 = 1.0;

    match i_params[0] {
        // The typical use case: inside the egg
        1 => {
            let r = r_params[0];
            mu = r_params[1];
            gam_axis = r_params[2];
            e_ratio = r_params[3];

            let mut gam_emis = (gam_axis.powi(2) * y_input.powf(k_cbm - 3.0)).sqrt();

            if offaxis {
                gam_emis /= e_ratio.sqrt();
            }

            chiy = y_input + 2.0 * (4.0 - k_cbm) * gam_emis.powi(2) * (y_input - r);
        }

        // An alternate use case: along the shell given mu
        2 => {
            mu = r_params[0];
            gam_axis = r_params[1];
            chiy = y_input;
        }

        // Another use case: along shell given R_perp
        3 => {
            let isign = i_params[1] as f64;
            let r = r_params[0]; // R_perp as used here

            // Step to avert roundoff-triggered NaNs
            let mut mu_work = y_input.powi(2) - r.powi(2);

            if mu_work < -1.0e-17 {
                panic!(
                    "ERROR in gps_egg: roundoff inacceptable\ny_input**2 - r**2 = {mu_work:e}\nStopping program now"
                );
            } else {
                mu_work = mu_work.max(0.0);
            }

            mu = isign * mu_work.sqrt() / y_input;
            gam_axis = r_params[1];
            chiy = y_input;
        }

        _ => {
            panic!(
                "ERROR in gps_egg: i_params(1) not a valid choice: {}\n        1: solve for y inside shell\n        2,3: solve for y at shell\nStopping program now",
                i_params[0]
            );
        }
    }

    let mut chi_term = chiy * y_input.powf(3.0 - k_cbm);

    if offaxis {
        chi_term *= e_ratio;
    }

    2.0 * (4.0 - k_cbm) * gam_axis.powi(2) * y_input * (1.0 - mu) + chi_term - 1.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_inputs(k_cbm: i32) -> Inputs {
        Inputs {
            k_cbm,
            ..Default::default()
        }
    }

    #[test]
    fn case_1_inside_egg_on_axis_matches_manual_equation() {
        let inputs = base_inputs(0);

        let mut i_params: [i32; 7] = [0; 7];
        i_params[0] = 1; // inside the egg
        i_params[4] = 1; // on-axis

        let mut r_params: [f64; 7] = [0.0; 7];
        r_params[0] = 0.5; // r
        r_params[1] = 0.8; // mu
        r_params[2] = 2.0; // Gam_axis
        r_params[3] = 1.0; // E_ratio

        let y_input: f64 = 0.75;
        let k_cbm: f64 = inputs.k_cbm as f64;

        let gam_emis: f64 = (r_params[2].powi(2) * y_input.powf(k_cbm - 3.0)).sqrt();
        let chiy: f64 = y_input + 2.0 * (4.0 - k_cbm) * gam_emis.powi(2) * (y_input - r_params[0]);

        let chi_term: f64 = chiy * y_input.powf(3.0 - k_cbm);

        let expected: f64 =
            2.0 * (4.0 - k_cbm) * r_params[2].powi(2) * y_input * (1.0 - r_params[1]) + chi_term
                - 1.0;

        let result: f64 = gps_egg(y_input, &i_params, &r_params, &inputs);

        assert!((result - expected).abs() < 1.0e-12);
    }

    #[test]
    fn case_2_along_shell_given_mu_matches_manual_equation() {
        let inputs = base_inputs(2);

        let mut i_params: [i32; 7] = [0; 7];
        i_params[0] = 2; // along shell given mu
        i_params[4] = 1; // on-axis

        let mut r_params: [f64; 7] = [0.0; 7];
        r_params[0] = 0.7; // mu
        r_params[1] = 3.0; // Gam_axis

        let y_input: f64 = 0.8;
        let k_cbm: f64 = inputs.k_cbm as f64;

        let mu: f64 = r_params[0];
        let gam_axis: f64 = r_params[1];
        let chiy: f64 = y_input;
        let chi_term: f64 = chiy * y_input.powf(3.0 - k_cbm);

        let expected: f64 =
            2.0 * (4.0 - k_cbm) * gam_axis.powi(2) * y_input * (1.0 - mu) + chi_term - 1.0;

        let result: f64 = gps_egg(y_input, &i_params, &r_params, &inputs);

        assert!((result - expected).abs() < 1.0e-12);
    }

    #[test]
    fn case_3_along_shell_given_rperp_matches_manual_equation() {
        let inputs = base_inputs(0);

        let mut i_params: [i32; 7] = [0; 7];
        i_params[0] = 3; // along shell given R_perp
        i_params[1] = 1; // isign
        i_params[4] = 1; // on-axis

        let mut r_params: [f64; 7] = [0.0; 7];
        r_params[0] = 0.6; // R_perp
        r_params[1] = 2.0; // Gam_axis

        let y_input: f64 = 0.8;
        let k_cbm: f64 = inputs.k_cbm as f64;

        let mu: f64 = (y_input.powi(2) - r_params[0].powi(2)).sqrt() / y_input;
        let gam_axis: f64 = r_params[1];
        let chiy: f64 = y_input;
        let chi_term: f64 = chiy * y_input.powf(3.0 - k_cbm);

        let expected: f64 =
            2.0 * (4.0 - k_cbm) * gam_axis.powi(2) * y_input * (1.0 - mu) + chi_term - 1.0;

        let result = gps_egg(y_input, &i_params, &r_params, &inputs);

        assert!((result - expected).abs() < 1.0e-12);
    }

    #[test]
    fn case_1_off_axis_applies_e_ratio_scaling() {
        let inputs = base_inputs(0);

        let mut i_params: [i32; 7] = [0; 7];
        i_params[0] = 1; // inside the egg
        i_params[4] = 0; // off-axis

        let mut r_params: [f64; 7] = [0.0; 7];
        r_params[0] = 0.5; // r
        r_params[1] = 0.8; // mu
        r_params[2] = 2.0; // Gam_axis
        r_params[3] = 4.0; // E_ratio

        let y_input: f64 = 0.75;
        let k_cbm: f64 = inputs.k_cbm as f64;

        let mut gam_emis: f64 = (r_params[2].powi(2) * y_input.powf(k_cbm - 3.0)).sqrt();
        gam_emis /= r_params[3].sqrt();

        let chiy = y_input + 2.0 * (4.0 - k_cbm) * gam_emis.powi(2) * (y_input - r_params[0]);

        let mut chi_term = chiy * y_input.powf(3.0 - k_cbm);
        chi_term *= r_params[3];

        let expected = 2.0 * (4.0 - k_cbm) * r_params[2].powi(2) * y_input * (1.0 - r_params[1])
            + chi_term
            - 1.0;

        let result = gps_egg(y_input, &i_params, &r_params, &inputs);

        assert!((result - expected).abs() < 1.0e-12);
    }

    #[test]
    #[should_panic(expected = "i_params(1) not a valid choice")]
    fn invalid_case_panics() {
        let inputs = base_inputs(0);

        let mut i_params = [0; 7];
        i_params[0] = 99;

        let r_params = [0.0; 7];

        gps_egg(0.8, &i_params, &r_params, &inputs);
    }

    #[test]
    #[should_panic(expected = "roundoff inacceptable")]
    fn case_3_panics_when_roundoff_is_too_negative() {
        let inputs = base_inputs(0);

        let mut i_params = [0; 7];
        i_params[0] = 3;
        i_params[1] = 1;
        i_params[4] = 1;

        let mut r_params = [0.0; 7];
        r_params[0] = 0.9; // R_perp > y_input
        r_params[1] = 2.0; // Gam_axis

        gps_egg(0.8, &i_params, &r_params, &inputs);
    }
}
