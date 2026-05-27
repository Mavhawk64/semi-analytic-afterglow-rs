use std::f64::consts::SQRT_2;

use crate::functions::e_of_theta::e_of_theta;
use crate::functions::gam_of_theta::gam_of_theta;
use crate::functions::gps_egg::gps_egg;
use crate::functions::theta_f::theta_f;
use crate::modules::constants::{CCGS, PII, RM_PROT, SIG_T, XME};
use crate::modules::globals::Globals;
use crate::modules::inputs::Inputs;
use crate::subroutines::dist_rebin::dist_rebin;
use crate::subroutines::newtons_method::{NewtonError, newtons_method};
use crate::subroutines::photon_alph::{PhotonAlphError, PhotonAlphInput, photon_alph};
use crate::subroutines::photon_ssc::{PhotonSscError, PhotonSscInput, photon_ssc};
use crate::subroutines::photon_syn::{PhotonSynError, photon_syn};
use crate::subroutines::set_elec_dist::{SetElecDistError, set_elec_dist};

const SQRT_8: f64 = 2.828_427_124_746_190_3;

pub struct RadTransferInput<'a> {
    pub s: f64,
    pub n_phot: usize,
    pub jnu: &'a mut [f64],
    pub alpha: &'a mut [f64],
    pub r_params: &'a [f64; 7],
    pub i_params: &'a [i32; 7],
    pub inputs: &'a Inputs,
    pub globals: &'a Globals,
}

#[derive(Debug, Clone, PartialEq)]
pub enum RadTransferError {
    JnuArrayWrongLength { len: usize, expected_len: usize },
    AlphaArrayWrongLength { len: usize, expected_len: usize },
    SetElecDist(SetElecDistError),
    PhotonSyn(PhotonSynError),
    PhotonAlph(PhotonAlphError),
    PhotonSsc(PhotonSscError),
    Newton(NewtonError),
}

impl From<SetElecDistError> for RadTransferError {
    fn from(err: SetElecDistError) -> Self {
        RadTransferError::SetElecDist(err)
    }
}

impl From<PhotonSynError> for RadTransferError {
    fn from(err: PhotonSynError) -> Self {
        RadTransferError::PhotonSyn(err)
    }
}

impl From<PhotonAlphError> for RadTransferError {
    fn from(err: PhotonAlphError) -> Self {
        RadTransferError::PhotonAlph(err)
    }
}

impl From<PhotonSscError> for RadTransferError {
    fn from(err: PhotonSscError) -> Self {
        RadTransferError::PhotonSsc(err)
    }
}

impl From<NewtonError> for RadTransferError {
    fn from(err: NewtonError) -> Self {
        RadTransferError::Newton(err)
    }
}

/// Computes two quantities needed for evaluating the right hand side of the
///   standard radiative transfer equation:
///
/// $$
/// \frac{dI_\nu}{ds} = j_\nu - \alpha_\nu I_\nu
/// $$
///
///   at a particular value of s.  (Specifically, finds j_nu and alpha_nu.)
///
/// Does this in ten steps:
///
///   1. Read in parameters from the *_params arrays
///   2. Compute key quantities at the current location in the GPS99 egg
///   3. Determine the shock conditions at the time this fluid element
///      originally crossed the shock
///   4. Create the uncooled electron distribution
///   5. determine local (decayed) magnetic field value
///   6. Cool the electron distribution to the current value of chi
///   7. compute synchrotron production as function of photon frequency
///   8. compute absorption coefficient as function of photon frequency
///   9. compute synchrotron self-Compton production as function of photon
///      frequency
///   10. convert from physical quantities needed for dI/ds into normalized
///       quantities needed for the Euler method
///
/// We will make regular use of the Lorentz-invariance of the quantities
///   j_nu/nu^2 and alpha_nu*nu, as both of these should be evaluated in the
///   rest frame of the emitting/absorbing plasma, but we wish to determine
///   I_nu in the explosion frame.
///
/// # Arguments:
///
/// `input`:
/// * `s` - current location along line of sight. Normalized against R_LOS, so
///   it lies in [0,1]
/// * `n_phot` - length of photon arrays.  Will be compared against parameter
///   n_ph to make sure things match up correctly
/// * `jnu` - value of j_nu at this position
/// * `alpha` - value of alpha_nu at this position
/// * `r_params` - array (length 7) of up to 7 floating point parameters to
///   be used in evaluating the function
/// * `i_params` - array (length 7) of up to 7 integer parameters to be used
///   in evaluating the function
///   Meaning of `i_params` elements
///    1. `i_tm`
///    2. `k_x`
///    3. `n_step`
///    4. `i_endpoint`
///    5. `offaxis`
/// * `inputs` - runtime inputs
/// * `globals` - runtime global variables
pub fn rad_transfer_diffeq(input: RadTransferInput<'_>) -> Result<(), RadTransferError> {
    let RadTransferInput {
        s,
        n_phot,
        jnu,
        alpha,
        r_params,
        i_params,
        inputs,
        globals,
    } = input;
    if jnu.len() != n_phot {
        return Err(RadTransferError::JnuArrayWrongLength {
            len: jnu.len(),
            expected_len: n_phot,
        });
    }

    if alpha.len() != n_phot {
        return Err(RadTransferError::AlphaArrayWrongLength {
            len: alpha.len(),
            expected_len: n_phot,
        });
    }

    //--------------------------------------------------------------------------
    // (1) Read in useful constants from the *_params arrays
    //--------------------------------------------------------------------------

    let i_tm = i_params[0] as usize;
    let k_x = i_params[1];
    let _n_step = i_params[2];
    let i_endpoint = i_params[3];
    let offaxis = i_params[4] != 1;

    let r_perp = r_params[0];
    let _r_perp_max = r_params[1];

    let mut psi_fl = 0.0;
    let mut s_far = 0.0;
    let mut s_near = 0.0;

    // This is flipped from the Fortran code
    // I think it was a bug.
    // TODO: Double check this
    if !offaxis {
        s_far = r_params[2];
        s_near = r_params[3];
    } else {
        psi_fl = r_params[2];
    }

    let del_s = r_params[6];

    //--------------------------------------------------------------------------
    // i_params/r_params handled
    //--------------------------------------------------------------------------

    // Propagate down the on-axis/off-axis switch
    let mut i_params_int = [0i32; 7];
    let mut r_params_int = [0.0f64; 7];

    i_params_int[4] = i_params[4];

    //--------------------------------------------------------------------------
    // (2) Using current shock conditions along line of sight, compute several
    // key quantities at this location
    //--------------------------------------------------------------------------

    // Grab the hydrodynamic state of the forward shock along the jet axis
    let r_los = globals.r_los_array[i_tm];
    let gam_los = globals.gam_los_array[i_tm];

    // Fluid element's location within the jet at the time of emission; for the
    // moment, these remain normalized to the radius of the forward shock
    // along the line of sight
    let r = (s.powi(2) + r_perp.powi(2)).sqrt();
    let mu = s / r;

    let mut theta_fl = 0.0;
    let mut e_ratio = 1.0;

    if offaxis {
        // Use mu and the various angles to compute the amount of energy released in
        // the direction of this fluid element
        theta_fl = theta_f(inputs.theta_obs, psi_fl, mu);
        e_ratio = e_of_theta(inputs.theta_obs, inputs) / e_of_theta(theta_fl, inputs);
    }

    // y = R_emis / R_LOS using GS2002 Eq (A16), given the value of mu we found.
    // Note that s can equal s_far or s_near identically, in which case we can
    // set y = r directly without calling Newton_method
    // When dealing with a structured jet, it is possible for both r and y to be
    // greater than unity. This means the usual bounds for Newton's method
    // (between r and 1) will fail to contain the solution for y. Fortunately
    // GPS_egg is a strictly increasing function of y, so we can easily check
    // to see if we have the proper bounds
    let bound_check = if offaxis {
        s != s_far && s != s_near
    } else {
        i_endpoint == 0 || k_x == 0
    };

    let y = if bound_check {
        i_params_int[0] = 1; // To compute y inside the egg
        r_params_int[0] = r;
        r_params_int[1] = mu;
        r_params_int[2] = gam_los;

        if offaxis {
            r_params_int[3] = e_ratio;

            if gps_egg(1.0, &i_params_int, &r_params_int, inputs) > 0.0 {
                // We can use the normal bounds because solution lies between r and 1
                newtons_method(
                    |y_input, i_params, r_params| gps_egg(y_input, i_params, r_params, inputs),
                    &i_params_int,
                    &r_params_int,
                    0.5 * (r + 1.0),
                    r,
                    1.0,
                )?
            } else {
                // We need to let Newton's method search values larger than 1
                newtons_method(
                    |y_input, i_params, r_params| gps_egg(y_input, i_params, r_params, inputs),
                    &i_params_int,
                    &r_params_int,
                    r,
                    0.99 * r,
                    10.0,
                )?
            }
        } else {
            newtons_method(
                |y_input, i_params, r_params| gps_egg(y_input, i_params, r_params, inputs),
                &i_params_int,
                &r_params_int,
                0.5 * (s_far + 1.0),
                s_far,
                1.0,
            )?
        }
    } else {
        r
    };

    // Use y to find Gam_emis
    let k_cbm = inputs.k_cbm as f64;

    let mut gam_emis = (gam_los.powi(2) * y.powf(k_cbm - 3.0)).sqrt();

    if offaxis {
        gam_emis /= e_ratio.sqrt();
    }

    // Set an upper limit of Gam_emis since GRB shocks can't have arbitrarily
    // large Lorentz factors?
    if inputs.do_gammax {
        if offaxis {
            let gam_coast = gam_of_theta(theta_fl, inputs);
            gam_emis = gam_emis.min(gam_coast);
        } else {
            gam_emis = gam_emis.min(inputs.eta);
        }
    }

    // Find chi_emis using the BM76 definition of chi, since we have both r and
    // R in units of R_LOS
    let mut chi_emis = 1.0 + 2.0 * (4.0 - k_cbm) * gam_emis.powi(2) * (y - r) / y;

    // Roundoff errors can rarely cause issues near the forward shock. Prevent
    // those by placing a floor on chi_emis. Write statement simply alerts
    // that such an event occurred
    if chi_emis < 1.0 {
        eprintln!("chi_emis < 1: {chi_emis:24.14e}");

        if (1.0 - chi_emis) * gam_emis.powi(2) > 1.0e-1 {
            eprintln!("    WARNING: chi_emis disparity comparable to 1/Gam^2");
        }
    }

    chi_emis = chi_emis.max(1.0);

    // Compute little_gam from GS2002 Eq (A9)
    let little_gam = gam_emis / (2.0 * chi_emis).sqrt();

    // If we are far enough from the shock, the BM76 solution breaks down and
    // allows for little_gam < 1. We choose to ignore these regions and
    // assume they contribute negligible emission/absorption to the overall
    // afterglow
    if little_gam < 1.05 {
        jnu.fill(1.0e-99);
        alpha.fill(1.0e-99);
        return Ok(());
    }

    // Compute the Doppler factor between the local plasma frame and the GRB
    // rest frame -- that is, do *not* include cosmological redshift
    let beta = (1.0 - little_gam.powi(-2)).sqrt();
    let o_o_dopp_fac = little_gam * (1.0 - beta * mu);
    let dopp_fac = 1.0 / o_o_dopp_fac;

    //--------------------------------------------------------------------------
    // Current conditions set
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (3) Determine shock conditions at the time this fluid element originally
    // crossed the shock
    // Note that R_zero must be expressed in physical units, not normalized
    // units, so that density equation makes sense.
    // Note also that we are *not* using GS2002 Eq (A15) to find Gam_zero,
    // since t_zero is in the explosion frame and is not the observer time.
    // Instead we're going back to Blandford & McKee (1976) Eq (69).
    //--------------------------------------------------------------------------

    // y = R/R_LOS
    // chi = (R/R_0)^(4-k)    <== GS2002 Eq (A8)
    // R_0 = c*t_0
    let r_zero = r_los * y * chi_emis.powf(-1.0 / (4.0 - k_cbm));
    let t_zero = r_zero / CCGS;

    let n_ext = if inputs.k_cbm == 0 {
        inputs.n_cbm
    } else {
        inputs.a_star * r_zero.powf(-k_cbm)
    };

    let energy = e_of_theta(theta_fl, inputs);

    // Equation (1) of Warren et al. (2022) has R^(3-k). It is just 3 here
    // because we are including the R^k in n_ext above
    let mut gam_zero =
        (17.0 - 4.0 * k_cbm) * energy / (8.0 * PII * n_ext * RM_PROT * r_zero.powi(3));

    gam_zero = gam_zero.sqrt();

    // Set an upper limit of Gam_zero since GRB shocks can't have arbitrarily
    // large Lorentz factors?
    if inputs.do_gammax {
        if offaxis {
            let gam_coast = gam_of_theta(theta_fl, inputs);
            gam_zero = gam_zero.min(gam_coast);
        } else {
            gam_zero = gam_zero.min(inputs.eta);
        }
    }

    let n_zero = SQRT_8 * gam_zero * n_ext;
    let e_zero = 2.0 * gam_zero.powi(2) * n_ext * RM_PROT;

    //--------------------------------------------------------------------------
    // Initial shock conditions found
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (4) Create the un-cooled electron spectrum. Do this in two steps.
    // (a) Find the distribution for a density of 1 cm^-3
    // (b) Rescale the distribution to the density at the emission location
    // Note that elec_dist is given as N(gamma)*dgamma, since that is the
    // quantity that will be conserved as the distribution cools.
    //--------------------------------------------------------------------------

    let mut elec_dist = set_elec_dist(
        gam_zero,
        n_ext,
        &globals.p_pf_cgs,
        &globals.elec_gam_array,
        globals.num_dp,
        inputs,
    )?;

    // Granot & Sari (2002) Equation (A9)
    let n_emis = n_zero * chi_emis.powf(-(13.0 - 2.0 * k_cbm) / (8.0 - 2.0 * k_cbm));

    // Rescale the spectrum to the calculated value of density
    for value in &mut elec_dist {
        if *value > 1.0e-55 {
            *value *= n_emis;
        }
    }

    //--------------------------------------------------------------------------
    // Initial electron spectrum set
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (5) Find the local magnetic field value. This is quite simple when
    // microturbulent decay is not considered, as it is the combination of
    // GS2002 Eq (A9) and Eq (A2)
    //--------------------------------------------------------------------------

    // Local magnetic field value
    let b_sq = 8.0
        * PII
        * inputs.eps_b0
        * e_zero
        * chi_emis.powf(-(26.0 - 4.0 * k_cbm) / (12.0 - 3.0 * k_cbm));

    //--------------------------------------------------------------------------
    // Local magnetic field determined
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (6) Cool the electron spectrum to its current value using GS2002's Eq.
    // (A11). This equation is accurate to order 1/gamma and allows for
    // electron Lorentz factors smaller than 1. However, given that we don't
    // care about the contributions of *extremely* low-energy electrons
    // (defined as below 3 MeV, or gamma < 6), this error is acceptable.
    //--------------------------------------------------------------------------

    let mut elec_gam_cooled = vec![0.0; globals.num_dp + 1];

    // The bottom of the electron distribution never cools, since it already has
    // zero energy
    elec_gam_cooled[0] = globals.elec_gam_array[0];

    // Find maximum possible electron Lorentz factor (corresponding to an
    // initially infinite energy) at this value of chi, using GS2002 Eq. (A12)
    let mut gam_max_cool = if chi_emis > 1.0 {
        SQRT_2
            * (19.0 - 2.0 * k_cbm)
            * PII
            * XME
            * CCGS
            * gam_zero
            * chi_emis.powf((25.0 - 2.0 * k_cbm) / (24.0 - 6.0 * k_cbm))
            / (SIG_T
                * 8.0
                * PII
                * inputs.eps_b0
                * e_zero
                * t_zero
                * (chi_emis.powf((19.0 - 2.0 * k_cbm) / (12.0 - 3.0 * k_cbm)) - 1.0))
    } else {
        globals.elec_gam_array[globals.num_dp]
    };

    // Since our initial electron array didn't actually extend to infinity, it
    // is possible for the above equation to return a value larger than the
    // max for sufficiently small values of chi. Correct for that here to
    // ensure our electrons don't *gain* energy while they "cool"
    gam_max_cool = gam_max_cool.min(globals.elec_gam_array[globals.num_dp]);

    // Apply GS2002's Eq. (A11)
    if chi_emis == 1.0 {
        elec_gam_cooled[1..=globals.num_dp]
            .copy_from_slice(&globals.elec_gam_array[1..=globals.num_dp]);
    } else {
        for (l_mom, gamma_cooled) in elec_gam_cooled
            .iter_mut()
            .enumerate()
            .take(globals.num_dp + 1)
            .skip(1)
        {
            *gamma_cooled = globals.elec_gam_array[l_mom]
                / (chi_emis.powf((13.0 - 2.0 * k_cbm) / (24.0 - 6.0 * k_cbm))
                    + globals.elec_gam_array[l_mom] / gam_max_cool);
        }
    }

    // Rebin the electron distribution if there has been substantial cooling
    let use_rebinned =
        elec_gam_cooled[globals.num_dp] < 0.1 * globals.elec_gam_array[globals.num_dp];

    let elec_dist_rebin = if use_rebinned {
        dist_rebin(
            globals.num_dp,
            &elec_gam_cooled,
            &elec_dist,
            &globals.elec_gam_array,
        )
    } else {
        vec![1.0e-99; globals.num_dp]
    };

    //--------------------------------------------------------------------------
    // Electron spectrum cooled
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (7) Compute synchrotron production in the plasma frame using the
    // Lorentz-invariant quantity j_nu / nu^2
    //--------------------------------------------------------------------------

    // Convert observer-frame photon energies into local plasma frame, now
    // accounting for the cosmological redshift
    let phot_en_pf: Vec<f64> = globals
        .phot_en_cgs
        .iter()
        .take(n_phot)
        .map(|phot_en| phot_en * o_o_dopp_fac * (1.0 + inputs.redshift))
        .collect();

    // Compute j_nu in the plasma frame. Use rebinned array if there has been
    // substantial cooling
    let mut syn_jnu = if use_rebinned {
        photon_syn(
            globals.num_dp,
            &globals.elec_gam_array,
            &elec_dist_rebin,
            b_sq,
            &phot_en_pf,
            i_params,
            globals,
        )?
    } else {
        photon_syn(
            globals.num_dp,
            &elec_gam_cooled,
            &elec_dist,
            b_sq,
            &phot_en_pf,
            i_params,
            globals,
        )?
    };

    //--------------------------------------------------------------------------
    // Synchrotron production calculated
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (8) Compute synchrotron self-absorption in the explosion frame using the
    // Lorentz-invariant quantity alpha_nu * nu
    //--------------------------------------------------------------------------

    // Explosion-frame photon spectrum has already been converted into the local
    // plasma frame, so we don't need to redo that here
    // However, we do need to convert the 1/nu^2 array, which will greatly
    // speed up the computation for alpha_nu
    let o_o_nu2_pf: Vec<f64> = globals
        .o_o_nu_sq
        .iter()
        .take(n_phot)
        .map(|value| value * dopp_fac.powi(2) / (1.0 + inputs.redshift).powi(2))
        .collect();

    // Compute alpha_nu in the plasma frame. Use rebinned array if there has
    // been substantial cooling
    let mut alpha_nu = if inputs.do_ssa {
        if use_rebinned {
            photon_alph(PhotonAlphInput {
                num_dp: globals.num_dp,
                elec_gam_array: &globals.elec_gam_array,
                dn_pf: &elec_dist_rebin,
                b_sq,
                phot_en_array: &phot_en_pf,
                o_o_nu_sq: &o_o_nu2_pf,
                i_params,
                globals,
            })?
        } else {
            photon_alph(PhotonAlphInput {
                num_dp: globals.num_dp,
                elec_gam_array: &elec_gam_cooled,
                dn_pf: &elec_dist,
                b_sq,
                phot_en_array: &phot_en_pf,
                o_o_nu_sq: &o_o_nu2_pf,
                i_params,
                globals,
            })?
        }
    } else {
        vec![1.0e-99; n_phot]
    };

    // Convert alpha_nu in the plasma frame to alpha_nu in the explosion frame
    for value in &mut alpha_nu {
        if *value > 1.0e-55 {
            *value /= dopp_fac;
        }
    }

    //--------------------------------------------------------------------------
    // Synchrotron self-absorption calculated
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (9) Compute synchrotron self-Compton production
    //--------------------------------------------------------------------------

    // Convert 1/alpha array into plasma frame (this is not the same alpha as
    // computed in step (8) above)
    let o_o_alpha_pf: Vec<f64> = globals
        .o_o_alpha_array
        .iter()
        .take(n_phot)
        .map(|value| value * dopp_fac / (1.0 + inputs.redshift))
        .collect();

    // Set the parameter that must be passed to photon_ssc: the plasma-frame
    // distance of the step we're currently taking in the Euler method
    r_params_int[0] = del_s * r_los / little_gam;

    // Compute j_nu in the plasma frame
    let mut ssc_jnu = if inputs.do_ssc {
        photon_ssc(PhotonSscInput {
            num_dp: globals.num_dp,
            elec_gam_array: &elec_gam_cooled,
            dn_pf: &elec_dist,
            phot_en_array: &phot_en_pf,
            o_o_alpha_array: &o_o_alpha_pf,
            syn_jnu: &syn_jnu,
            i_params: &i_params_int,
            r_params: &r_params_int,
            globals,
        })?
    } else {
        vec![1.0e-99; n_phot]
    };

    //--------------------------------------------------------------------------
    // Synchrotron self-Compton calculated
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (10) What we found above were the factors needed to compute dI/ds, the
    // rate of change in physical units. However, our Euler method uses
    // position scaled against the line-of-sight radius of the shock. Since
    // ds = R_LOS * ds', dI/ds' = dI/ds * R_LOS, and so our values of jnu and
    // alpha must be re-scaled accordingly.
    //--------------------------------------------------------------------------

    // Convert j_nu in the plasma frame to j_nu in the explosion frame
    for value in &mut syn_jnu {
        if *value > 1.0e-55 {
            *value *= dopp_fac.powi(2);
        }
    }

    for value in &mut ssc_jnu {
        if *value > 1.0e-55 {
            *value *= dopp_fac.powi(2);
        }
    }

    jnu.fill(1.0e-99);
    alpha.fill(1.0e-99);

    for m_ph in 0..n_phot {
        if syn_jnu[m_ph] > 1.0e-55 {
            jnu[m_ph] += syn_jnu[m_ph] * r_los;
        }

        if ssc_jnu[m_ph] > 1.0e-55 {
            jnu[m_ph] += ssc_jnu[m_ph] * r_los;
        }

        if alpha_nu[m_ph] > 1.0e-55 {
            alpha[m_ph] = alpha_nu[m_ph] * r_los;
        }
    }

    //--------------------------------------------------------------------------
    // dI/ds' found
    //--------------------------------------------------------------------------

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_inputs() -> Inputs {
        Inputs {
            k_cbm: 0,
            n_cbm: 1.0,
            a_star: 1.0,
            eta: 100.0,
            eps_b0: 0.01,
            do_ssa: false,
            do_ssc: false,
            do_gammax: false,
            redshift: 0.0,
            theta_obs: 0.0,
            do_jet_structure: false,
            e_iso: 1.0e52,
            pl_only: true,
            p_spec: 2.2,
            epse0: 0.1,
            do_pmax: false,
            ..Default::default()
        }
    }

    fn base_globals(num_dp: usize, num_ph: usize) -> Globals {
        Globals {
            num_dp,
            num_ph,
            p_pf_cgs: (0..=num_dp)
                .map(|i| 1.0e-20 * 10.0_f64.powf(i as f64))
                .collect(),
            elec_gam_array: (0..=num_dp).map(|i| 1.0 + 10.0 * i as f64).collect(),
            phot_en_cgs: (0..num_ph)
                .map(|i| 1.0e-12 * 2.0_f64.powf(i as f64))
                .collect(),
            o_o_nu_sq: vec![1.0; num_ph],
            o_o_alpha_array: vec![1.0; num_ph],
            r_los_array: vec![1.0e17],
            gam_los_array: vec![20.0],
            syn_x: (0..119).map(|i| i as f64).collect(),
            syn_f: (0..119).map(|i| 1.0 + i as f64).collect(),
            o_o_sone: 1.0,
            ..Default::default()
        }
    }

    fn base_params() -> ([i32; 7], [f64; 7]) {
        let mut i_params = [0; 7];
        i_params[0] = 0; // i_tm
        i_params[1] = 0; // k_x
        i_params[2] = 1; // n_step
        i_params[3] = 1; // i_endpoint
        i_params[4] = 1; // on-axis

        let mut r_params = [0.0; 7];
        r_params[0] = 0.1; // R_perp
        r_params[1] = 1.0; // R_perp_max
        r_params[2] = 0.5; // s_far for on-axis after flipped logic
        r_params[3] = 1.0; // s_near
        r_params[6] = 0.01; // del_s

        (i_params, r_params)
    }

    #[test]
    fn returns_error_when_jnu_length_is_wrong() {
        let inputs = base_inputs();
        let globals = base_globals(4, 4);
        let (i_params, r_params) = base_params();

        let mut jnu = vec![1.0e-99; 3];
        let mut alpha = vec![1.0e-99; 4];

        let result = rad_transfer_diffeq(RadTransferInput {
            s: 0.5,
            n_phot: 4,
            jnu: &mut jnu,
            alpha: &mut alpha,
            r_params: &r_params,
            i_params: &i_params,
            inputs: &inputs,
            globals: &globals,
        });

        assert!(matches!(
            result,
            Err(RadTransferError::JnuArrayWrongLength {
                len: 3,
                expected_len: 4
            })
        ));
    }

    #[test]
    fn returns_error_when_alpha_length_is_wrong() {
        let inputs = base_inputs();
        let globals = base_globals(4, 4);
        let (i_params, r_params) = base_params();

        let mut jnu = vec![1.0e-99; 4];
        let mut alpha = vec![1.0e-99; 3];

        let result = rad_transfer_diffeq(RadTransferInput {
            s: 0.5,
            n_phot: 4,
            jnu: &mut jnu,
            alpha: &mut alpha,
            r_params: &r_params,
            i_params: &i_params,
            inputs: &inputs,
            globals: &globals,
        });

        assert!(matches!(
            result,
            Err(RadTransferError::AlphaArrayWrongLength {
                len: 3,
                expected_len: 4
            })
        ));
    }

    #[test]
    fn low_little_gamma_returns_floor_values() {
        let inputs = base_inputs();
        let mut globals = base_globals(4, 4);
        globals.gam_los_array = vec![1.0];

        let (i_params, r_params) = base_params();

        let mut jnu = vec![0.0; 4];
        let mut alpha = vec![0.0; 4];

        let result = rad_transfer_diffeq(RadTransferInput {
            s: 0.5,
            n_phot: 4,
            jnu: &mut jnu,
            alpha: &mut alpha,
            r_params: &r_params,
            i_params: &i_params,
            inputs: &inputs,
            globals: &globals,
        });

        assert!(result.is_ok());
        assert!(jnu.iter().all(|value| (*value - 1.0e-99).abs() < 1.0e-120));
        assert!(
            alpha
                .iter()
                .all(|value| (*value - 1.0e-99).abs() < 1.0e-120)
        );
    }

    #[test]
    fn typical_on_axis_case_returns_finite_arrays() {
        let inputs = base_inputs();
        let globals = base_globals(4, 4);
        let (i_params, r_params) = base_params();

        let mut jnu = vec![0.0; 4];
        let mut alpha = vec![0.0; 4];

        let result = rad_transfer_diffeq(RadTransferInput {
            s: 0.5,
            n_phot: 4,
            jnu: &mut jnu,
            alpha: &mut alpha,
            r_params: &r_params,
            i_params: &i_params,
            inputs: &inputs,
            globals: &globals,
        });

        assert!(result.is_ok(), "result = {result:?}");
        assert!(jnu.iter().all(|value| value.is_finite()));
        assert!(alpha.iter().all(|value| value.is_finite()));
    }
}
