use crate::modules::inputs::Inputs;

const H0: f64 = 67.66;
const OMEGA_R: f64 = 0.4165 / (H0 * H0);
const OMEGA_VAC: f64 = 0.6889 - 0.5 * OMEGA_R;
const OMEGA_M: f64 = 0.3111 - 0.5 * OMEGA_R;
const C: f64 = 2.997_924_58e5;
const D_H: f64 = C / H0;
const T_H: f64 = 9.778e11 / H0;
const NUM_STEPS: usize = 1000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CosmoCalcResult {
    pub z: f64,
    pub d_l: f64,
    pub t_look: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CosmoCalcError {
    InvalidInputs { z: f64, d_l: f64 },
    RedshiftSearchExceededBounds { z_new: f64 },
}

/// Calculator to shift between redshift, lookback time, and luminosity distance.
/// Subroutine adapted from Hogg (1999). For more detail see
/// [http://ui.adsabs.harvard.edu/abs/1999astro.ph..5116H]
///
/// Note one major difference: Equation (14) in Hogg uses Omega_r where this
/// program uses Omega_k, and does not include a term for what this program
/// calls Omega_r. For more information, see Wright (2006):
/// [http://ui.adsabs.harvard.edu/abs/2006PASP..118.1711W]
///
/// The program calculates lookback time (in years) as well, but since this
/// result is not necessary for GRB redshift calculation it is not passed
/// back to the calling subroutine.
///
/// # Arguments
///
/// * `z` - redshift. Exactly one of `z` and `d_l` must be positive
/// * `d_l` - luminosity distance in Mpc. Exactly one of `z` and `d_l` must be positive
/// * `inputs` - runtime inputs, used for `debug_mode`
pub fn cosmo_calc(z: f64, d_l: f64, inputs: &Inputs) -> Result<CosmoCalcResult, CosmoCalcError> {
    // Local variables

    let z_save = z;
    let d_l_save = d_l;

    let mut z = z;
    let mut d_l = d_l;
    let mut t_look = 0.0;

    // Status of code
    //  -1: initialization value
    //   1: value of z for provided d_L located. Proceed to integration for
    //      lookback time
    //   2: d_L less than critical value. Skip stepping through z integral and
    //      integration to find t_look
    let mut code_stat = -1;

    //--------------------------------------------------------------------------
    // Set either redshift or comoving distance to a nonzero value. Other
    // should be set to zero
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Quick error check that at least one input is physically reasonable and
    // that only one is positive
    //--------------------------------------------------------------------------

    if ((z <= 0.0) && (d_l <= 0.0)) || ((z > 0.0) && (d_l > 0.0)) {
        return Err(CosmoCalcError::InvalidInputs {
            z: z_save,
            d_l: d_l_save,
        });
    }

    //--------------------------------------------------------------------------
    // If d_L is given, step through integral in z (Equation (14) in Hogg
    // 1999) until next step exceeds d_L
    //--------------------------------------------------------------------------

    if z == 0.0 && d_l > 0.0 && d_l < 0.443 {
        // First, check to make sure d_L is large enough to get a reasonable
        // redshift -- otherwise return zero
        z = 0.0;
        t_look = d_l * 3.2616e6; // convert Mpc to years

        // d_L less than critical value. Skip stepping through
        // z integral and integration to find t_look
        code_stat = 2;
    } else if z == 0.0 && d_l > 0.0 {
        // d_L was a reasonable value
        let mut d_cm_old = 0.0;
        let mut d_cm_new = 0.0;
        let mut z_old = 0.0;

        // E(z) = sqrt[Omega_r*(1+z)^4 + Omega_m*(1+z)^3
        //           + Omega_k*(1+z)^2 + Omega_vac]
        // Equation (13) in Hogg (1999)
        let mut e_old = 1.0;

        for i in 1..=100 {
            let z_new = 1.0e-4 * i as f64;

            // Equation (14) in Hogg (1999)
            let e_new =
                (OMEGA_R * (1.0 + z_new).powi(4) + OMEGA_M * (1.0 + z_new).powi(3) + OMEGA_VAC)
                    .sqrt();

            // Use Equation (15) in Hogg (1999) to update d_new
            d_cm_new += D_H * 1.0e-4 * 0.5 * (1.0 / e_old + 1.0 / e_new);

            // Use Equation (21) in Hogg (1999) to compute d_L from d_CM; remember,
            // assuming a flat Universe so d_M = d_C
            let d_l_old = (1.0 + z_old) * d_cm_old;
            let d_l_new = (1.0 + z_new) * d_cm_new;

            // If current step in z exceeded provided value of d_L, perform linear
            // interpolation between the dL_old and dL_new
            if d_l_old < d_l && d_l_new >= d_l {
                let slope = (d_l_new - d_l_old) / (z_new - z_old);
                z = (d_l - d_l_old) / slope + z_old;

                // value of z for provided d_L located. Proceed to
                // integration for lookback time
                code_stat = 1;
                break;
            }

            // Set variables for next pass through loop
            d_cm_old = d_cm_new;
            e_old = e_new;
            z_old = z_new;
        }

        // If a value of z wasn't already found stepping from z = 0 to 0.01 with
        // steps of size 1.d-4, start stepping up from z = 0.01 with steps of
        // size 0.01
        let mut i = 0usize;

        while code_stat < 1 {
            i += 1;

            let z_new = 1.0e-2 * i as f64;

            // Equation (14) in Hogg (1999)
            let e_new =
                (OMEGA_R * (1.0 + z_new).powi(4) + OMEGA_M * (1.0 + z_new).powi(3) + OMEGA_VAC)
                    .sqrt();

            // Use Equation (15) in Hogg (1999) to update d_new
            d_cm_new += D_H * 1.0e-2 * 0.5 * (1.0 / e_old + 1.0 / e_new);

            // Use Equation (21) in Hogg (1999) to compute d_L from d_CM; remember,
            // assuming a flat Universe so d_M = d_C
            let d_l_old = (1.0 + z_old) * d_cm_old;
            let d_l_new = (1.0 + z_new) * d_cm_new;

            // If current step in z exceeded provided value of d_L, perform linear
            // interpolation between the dL_old and dL_new
            if d_l_old < d_l && d_l_new >= d_l {
                let slope = (d_l_new - d_l_old) / (z_new - z_old);
                z = (d_l - d_l_old) / slope + z_old;

                // value of z for provided d_CM located. Proceed to
                // integration for lookback time
                code_stat = 1;
                break;
            }

            // Set variables for next pass through loop
            d_cm_old = d_cm_new;
            e_old = e_new;
            z_old = z_new;

            // Quick and dirty error check
            if i > 1000 {
                return Err(CosmoCalcError::RedshiftSearchExceededBounds { z_new });
            }
        }
    }

    //--------------------------------------------------------------------------
    // End calculation of z for a provided d_L
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // If z is given, perform integral to find d_L using the trapezoid rule
    // and Equation (15) in Hogg (1999).
    // Independently of the previous, calculate t_look using the trapezoid rule
    // and Equations (30) in Hogg (1999).
    //--------------------------------------------------------------------------

    if (z > 0.0 && d_l == 0.0) || code_stat == 1 {
        let mut z_old = 0.0;

        // E(z) = sqrt[Omega_m*(1+z)^3 + Omega_r*(1+z)^2 + Omega_vac]
        // Equation (14) in Hogg (1999)
        let mut e_old = 1.0;
        let mut d_cm = 0.0;

        t_look = 0.0;

        for i in 1..=NUM_STEPS {
            let z_new = z * i as f64 / NUM_STEPS as f64;

            // Equation (14) in Hogg (1999)
            let e_new =
                (OMEGA_R * (1.0 + z_new).powi(4) + OMEGA_M * (1.0 + z_new).powi(3) + OMEGA_VAC)
                    .sqrt();

            // If d_CM wasn't provided, use Equation (15) in Hogg (1999) to update
            // running total
            if code_stat < 1 {
                d_cm += 0.5 * (1.0 / e_old + 1.0 / e_new);
            }

            // Independently of previous calculation, use Equation (30) in Hogg (1999)
            // to update t_look
            t_look += 0.5 * (1.0 / ((1.0 + z_old) * e_old) + 1.0 / ((1.0 + z_new) * e_new));

            // Set variables for next pass through loop
            z_old = z_new;
            e_old = e_new;
        }

        // Finally, multiply by appropriate scale factors out front and width of
        // the trapezoids used for integration...
        if code_stat < 1 {
            d_cm *= D_H * z / NUM_STEPS as f64;

            // ...and convert comoving distance to luminosity distance
            d_l = (1.0 + z) * d_cm;
        }

        t_look *= T_H * z / NUM_STEPS as f64;
    }

    //--------------------------------------------------------------------------
    // End calculation of t_look, and possibly also d_L
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Debugging output lines
    //--------------------------------------------------------------------------

    if inputs.debug_mode {
        eprintln!("Program finished.");
        eprintln!("Inputted values were z = {z_save:11.4e} and d_L = {d_l_save:11.4e}");
        eprintln!();
        eprintln!("Calculated values: z = {z:11.4e}");
        eprintln!("                d_L  = {d_l:11.4e} Mpc");
        eprintln!("              t_look = {t_look:11.4e} years");
    }

    Ok(CosmoCalcResult { z, d_l, t_look })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn inputs() -> Inputs {
        Inputs {
            debug_mode: false,
            ..Default::default()
        }
    }

    #[test]
    fn rejects_both_inputs_zero() {
        let result = cosmo_calc(0.0, 0.0, &inputs());

        assert!(matches!(
            result,
            Err(CosmoCalcError::InvalidInputs { z: 0.0, d_l: 0.0 })
        ));
    }

    #[test]
    fn rejects_both_inputs_positive() {
        let result = cosmo_calc(1.0, 100.0, &inputs());

        assert!(matches!(
            result,
            Err(CosmoCalcError::InvalidInputs { z: 1.0, d_l: 100.0 })
        ));
    }

    #[test]
    fn computes_luminosity_distance_from_redshift() {
        let result = cosmo_calc(1.0, 0.0, &inputs()).expect("cosmo_calc should succeed");

        assert!(result.z > 0.0);
        assert!(result.d_l > 0.0);
        assert!(result.t_look > 0.0);
    }

    #[test]
    fn computes_redshift_from_luminosity_distance() {
        let result = cosmo_calc(0.0, 1_000.0, &inputs()).expect("cosmo_calc should succeed");

        assert!(result.z > 0.0);
        assert_eq!(result.d_l, 1_000.0);
        assert!(result.t_look > 0.0);
    }

    #[test]
    fn very_small_luminosity_distance_returns_zero_redshift() {
        let result = cosmo_calc(0.0, 0.1, &inputs()).expect("cosmo_calc should succeed");

        assert_eq!(result.z, 0.0);
        assert_eq!(result.d_l, 0.1);
        assert_eq!(result.t_look, 0.1 * 3.2616e6);
    }
}
