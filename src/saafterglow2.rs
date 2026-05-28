use std::env;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::time::Instant;

use crate::functions::d_rperp_dmu::d_rperp_dmu;
use crate::functions::e_of_theta::e_of_theta;
use crate::functions::fixed_rperp::fixed_rperp;
use crate::functions::gam_of_theta::gam_of_theta;
use crate::functions::r_of_theta::r_of_theta;
use crate::modules::constants::{
    CCGS, ERG_TO_EV, MEV_TO_ERG, PC_TO_CM, PII, RM_ELEC, RM_PROT, XH, XME, XMP,
};
use crate::modules::globals::{Globals, P_ARRAY_MAX, P_ARRAY_MIN, bin_factor_mom, bin_factor_phot};
use crate::modules::opacity::{OpacityTable, calculate_absorptions};
use crate::modules::parameters::{BINS_PER_DECADE_MOM, BINS_PER_DECADE_PHOT, N_DP, N_PH};
use crate::subroutines::brents_minimum::{BrentError, brents_minimum};
use crate::subroutines::bw_euler_stepper::bw_euler_stepper;
use crate::subroutines::calc_syn_fx::calc_syn_fx;
use crate::subroutines::cosmo_calc::{CosmoCalcError, cosmo_calc};
use crate::subroutines::data_input::{DataInputError, data_input_from_str};
use crate::subroutines::newtons_method::{NewtonError, newtons_method};
use crate::subroutines::rad_transfer_diffeq::{RadTransferInput, rad_transfer_diffeq};
use crate::subroutines::set_em_step_sizes::set_em_step_sizes;
use crate::subroutines::set_x_array_1d::set_x_array_1d;

/// Computes broadband afterglow emission numerically, by combining five parts
///
/// 1. A prescription (the Blandford-McKee solution) for the bulk
///    hydrodynamics at all times in the explosion frame
/// 2. A formula for magnetic field strength as a function of position and
///    time, including microturbulent amplification and decay
/// 3. An analytical description for the electron distribution at the shock
/// 4. A method to calculate downstream evolution of the electron dist
/// 5. Photon production and absorption calculations, including relevant
///    relativistic (e.g. beaming) and geometric (e.g. line of sight)
///    considerations
///
/// Paper notations:
///
/// 1. BM76: Blandford & McKee (1976PhFl...19.1130B)
/// 2. GS2002: Granot & Sari (2002ApJ...568..820G)
/// 3. GPS99: Granot, Piran & Sari (1999ApJ...513..679G, but also 1999ApJ...527..236G)
#[derive(Debug)]
pub enum SaaError {
    Io(std::io::Error),
    DataInput(DataInputError),
    CosmoCalc(CosmoCalcError),
    Newton(NewtonError),
    Brent(BrentError),
    InvalidArgumentCount {
        count: usize,
    },
    InvalidDistanceRedshiftPair,
    InvalidStartPair,
    InvalidEndPair,
    NonContinuousElectronDistribution {
        nt_norm: f64,
        p_spec: f64,
        max_nt_norm: f64,
    },
    TimeRangeInvalid {
        t_start: f64,
        t_end: f64,
    },
    GammaRangeInvalid {
        gam_start: f64,
        gam_end: f64,
    },
    PhotonArrayTooSmall {
        num_ph: usize,
        n_ph: usize,
    },
    InvalidOffAxisGeometry {
        r_far: f64,
        r_perp: f64,
    },
}

impl From<std::io::Error> for SaaError {
    fn from(err: std::io::Error) -> Self {
        SaaError::Io(err)
    }
}

impl From<DataInputError> for SaaError {
    fn from(err: DataInputError) -> Self {
        SaaError::DataInput(err)
    }
}

impl From<CosmoCalcError> for SaaError {
    fn from(err: CosmoCalcError) -> Self {
        SaaError::CosmoCalc(err)
    }
}

impl From<NewtonError> for SaaError {
    fn from(err: NewtonError) -> Self {
        SaaError::Newton(err)
    }
}

impl From<BrentError> for SaaError {
    fn from(err: BrentError) -> Self {
        SaaError::Brent(err)
    }
}

const DEGTRD: f64 = std::f64::consts::PI / 180.0;
const RADTDG: f64 = 180.0 / std::f64::consts::PI;

pub fn run() -> Result<(), SaaError> {
    // Timing call
    let start = Instant::now();

    //--------------------------------------------------------------------------
    // Under normal operation, no runtime arguments will be used, and so just
    // open the output files
    // When used in a bulk submission script on AI.Panther, a single runtime
    // argument will be passed -- the job number within the bulk script --
    // which will be used to adjust the names of the input and output files
    //--------------------------------------------------------------------------

    let args: Vec<String> = env::args().skip(1).collect();

    let (input_text, mut out_writer, mut sed_writer): (String, BufWriter<File>, BufWriter<File>) =
        match args.len() {
            0 => {
                // Normal operation
                let mut input_text = String::new();
                std::io::stdin().read_to_string(&mut input_text)?;

                (
                    input_text,
                    BufWriter::new(File::create("saa_out.dat")?),
                    BufWriter::new(File::create("saa_SEDs.dat")?),
                )
            }
            1 => {
                // This argument is the job number associated with the bulk script
                let argument = &args[0];

                let input_file = format!("saa_input_{argument}.dat");
                let output_file = format!("saa_out_{argument}.dat");
                let sed_file = format!("saa_SEDs_{argument}.dat");

                let mut input_text = String::new();
                File::open(input_file)?.read_to_string(&mut input_text)?;

                (
                    input_text,
                    BufWriter::new(File::create(output_file)?),
                    BufWriter::new(File::create(sed_file)?),
                )
            }
            count => return Err(SaaError::InvalidArgumentCount { count }),
        };

    let mut image_1046 = BufWriter::new(File::create("fort.1046")?);
    let mut image_1070 = BufWriter::new(File::create("fort.1070")?);
    let mut image_1089 = BufWriter::new(File::create("fort.1089")?);
    let mut image_1113 = BufWriter::new(File::create("fort.1113")?);
    let mut image_1161 = BufWriter::new(File::create("fort.1161")?);
    let mut image_1185 = BufWriter::new(File::create("fort.1185")?);

    //--------------------------------------------------------------------------
    // Read in data from the input file, allocate the arrays for shock macro
    // states, set a density variable that depends on k_cbm, and ensure that
    // nt_norm and p_spec allow for a continuous solution--see discussion
    // surrounding Eq. (19) in Ressler & Laskar (2017ApJ...845..150R).
    //--------------------------------------------------------------------------

    let mut inputs = data_input_from_str(&input_text)?;

    let mut globals = Globals {
        t_obs_array: vec![0.0; inputs.n_tm],
        t_los_array: vec![0.0; inputs.n_tm],
        r_los_array: vec![0.0; inputs.n_tm],
        gam_los_array: vec![0.0; inputs.n_tm],
        ..Default::default()
    };

    // Instead of carrying around n_cbm and A_star everywhere, use a term that
    // is agnostic to the value of k_cbm and compute everything generically.
    let den_n = if inputs.k_cbm == 0 {
        inputs.n_cbm
    } else {
        inputs.a_star
    };

    // Ensure electron distribution can be continuous
    let max_nt_norm = 3.0 / (inputs.p_spec + 2.0);
    if inputs.nt_norm > max_nt_norm {
        return Err(SaaError::NonContinuousElectronDistribution {
            nt_norm: inputs.nt_norm,
            p_spec: inputs.p_spec,
            max_nt_norm,
        });
    }

    // Convert degrees to radians for jet structure since it's easier to work
    //   with radians for rest of calculation
    inputs.theta_j *= DEGTRD;
    inputs.theta_c *= DEGTRD;
    inputs.theta_obs *= DEGTRD;

    let energy = e_of_theta(inputs.theta_obs, &inputs);

    //--------------------------------------------------------------------------
    // Data read-in complete
    //--------------------------------------------------------------------------

    // Find whichever of redshift or comoving distance was not specified; must
    // convert dist_lum to megaparsecs since that's what cosmo_calc expects
    if (inputs.dist_lum == 0.0 && inputs.redshift > 0.0)
        || (inputs.dist_lum > 0.0 && inputs.redshift == 0.0)
    {
        let dist_lum_mpc = inputs.dist_lum / (PC_TO_CM * 1.0e6);
        let cosmo = cosmo_calc(inputs.redshift, dist_lum_mpc, &inputs)?;

        inputs.redshift = cosmo.z;

        if inputs.dist_lum == 0.0 {
            inputs.dist_lum = cosmo.d_l * (1.0e6 * PC_TO_CM);
        }
    } else {
        return Err(SaaError::InvalidDistanceRedshiftPair);
    }

    //--------------------------------------------------------------------------
    // Set up the hydro data arrays with values:
    // (1) Ensure that input is in the format desired
    // (2) Fill in other parts of input data
    // (3) Check for physical reasonableness of input data
    // (4) Fill in arrays for time, radius, Lorentz factor
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (1) Check t_start/Gam_start and t_end/Gam_end to make sure both pairs are
    // set as expected: one element positive, the other zero
    //--------------------------------------------------------------------------

    // Check on t_start/Gam_start
    if !((inputs.t_start == 0.0 && inputs.gam_start > 0.0)
        || (inputs.t_start > 0.0 && inputs.gam_start == 0.0))
    {
        return Err(SaaError::InvalidStartPair);
    }

    // Check on t_end/Gam_end
    if !((inputs.t_end == 0.0 && inputs.gam_end > 0.0)
        || (inputs.t_end > 0.0 && inputs.gam_end == 0.0))
    {
        return Err(SaaError::InvalidEndPair);
    }

    //--------------------------------------------------------------------------
    // (2) Compute whichever of t/Gam_start and t/Gam_end wasn't specified. Use
    // Eq. (A15) of GS2002, which itself is built from BM76's Eq. (69) and the
    // relation R = 2*(4-k)*Gam^2*c*t_obs. Note that (A15) uses little gamma
    // (the Lorentz factor of the fluid behind the shock), not capital Gamma
    // (the shock Lorentz factor); Gamma = gamma*sqrt(2) for Gamma >> 1.
    // This formulation has an error term on the order of 1/[2*(4-k)*Gam^2],
    // which is ~ a few percent when Gamma = 2.5. This is an acceptable
    // trade-off for easier manipulation of the equations.
    //--------------------------------------------------------------------------

    let k_cbm = inputs.k_cbm as f64;

    // Handle t/Gam_start
    if inputs.gam_start > 0.0 {
        inputs.t_start = energy * (17.0 - 4.0 * k_cbm)
            / (8.0 * PII * den_n * RM_PROT * inputs.gam_start.powf(8.0 - 2.0 * k_cbm));

        inputs.t_start = inputs.t_start.powf(1.0 / (3.0 - k_cbm))
            / (2.0 * (4.0 - k_cbm) * CCGS)
            / (1.0 + inputs.redshift);
    } else {
        // t_start specified
        inputs.gam_start = energy * (17.0 - 4.0 * k_cbm)
            / (8.0
                * PII
                * den_n
                * RM_PROT
                * (2.0 * (4.0 - k_cbm) * CCGS * inputs.t_start / (1.0 + inputs.redshift))
                    .powf(3.0 - k_cbm));

        inputs.gam_start = inputs.gam_start.powf(0.5 / (4.0 - k_cbm));
    }

    // Handle t/Gam_end
    if inputs.gam_end > 0.0 {
        inputs.t_end = energy * (17.0 - 4.0 * k_cbm)
            / (8.0 * PII * den_n * RM_PROT * inputs.gam_end.powf(8.0 - 2.0 * k_cbm));

        inputs.t_end = inputs.t_end.powf(1.0 / (3.0 - k_cbm))
            / (2.0 * (4.0 - k_cbm) * CCGS)
            / (1.0 + inputs.redshift);
    } else {
        // t_end specified
        inputs.gam_end = energy * (17.0 - 4.0 * k_cbm)
            / (8.0
                * PII
                * den_n
                * RM_PROT
                * (2.0 * (4.0 - k_cbm) * CCGS * inputs.t_end / (1.0 + inputs.redshift))
                    .powf(3.0 - k_cbm));

        inputs.gam_end = inputs.gam_end.powf(0.5 / (4.0 - k_cbm));
    }

    //--------------------------------------------------------------------------
    // (3) Make sure values of t/Gam_start and t/Gam_end are both sensible and
    // physically plausible:
    // (a) t_end > t_start, Gam_end < Gam_start
    // (b) Gam_start <= eta (since eta is the peak Lorentz factor
    //     reached by the shock)
    // (c) Gam_end >= 2.5 (to limit 1/Gam^2 error)
    //--------------------------------------------------------------------------

    // (a) t_end > t_start, Gam_end < Gam_start
    if inputs.t_end <= inputs.t_start {
        return Err(SaaError::TimeRangeInvalid {
            t_start: inputs.t_start,
            t_end: inputs.t_end,
        });
    }

    if inputs.gam_start <= inputs.gam_end {
        return Err(SaaError::GammaRangeInvalid {
            gam_start: inputs.gam_start,
            gam_end: inputs.gam_end,
        });
    }

    // (b) Gam_start <= eta(theta) (since eta is the coasting Lorentz factor
    //     reached by the shock)
    let gam_coast = gam_of_theta(inputs.theta_obs, &inputs);
    if inputs.gam_start > gam_coast {
        if inputs.do_progress_log {
            eprintln!("WARNING: Gam_start > Gam_coast");
            eprintln!(
                "      Gam_start: {:12.3e}; Gam_coast: {:12.3e}",
                inputs.gam_start, gam_coast
            );
            eprintln!("Replacing Gam_start with Gam_coast");
        }

        writeln!(out_writer, "WARNING: Gam_start > Gam_coast")?;
        writeln!(
            out_writer,
            "      Gam_start: {:12.3e}; Gam_coast: {:12.3e}",
            inputs.gam_start, gam_coast
        )?;
        writeln!(out_writer, "Replacing Gam_start with Gam_coast")?;

        inputs.gam_start = gam_coast;

        inputs.t_start = energy * (17.0 - 4.0 * k_cbm)
            / (8.0 * PII * den_n * RM_PROT * gam_coast.powf(8.0 - 2.0 * k_cbm));

        inputs.t_start = inputs.t_start.powf(1.0 / (3.0 - k_cbm)) / (2.0 * (4.0 - k_cbm) * CCGS)
            * (1.0 + inputs.redshift);
    }

    // (c) Gam_end >= 2.5
    if inputs.gam_end < 2.5 {
        if inputs.do_progress_log {
            eprintln!("WARNING: Gam_end < 2.5");
            eprintln!("      Gam_end: {:12.3e}", inputs.gam_end);
            eprintln!("Replacing Gam_end with 2.5");
        }

        writeln!(out_writer, "WARNING: Gam_end < 2.5")?;
        writeln!(out_writer, "      Gam_end: {:12.3e}", inputs.gam_end)?;
        writeln!(out_writer, "Replacing Gam_end with 2.5")?;
        writeln!(out_writer)?;

        inputs.gam_end = 2.5;

        inputs.t_end = energy * (17.0 - 4.0 * k_cbm)
            / (8.0 * PII * den_n * RM_PROT * 2.5_f64.powf(8.0 - 2.0 * k_cbm));

        inputs.t_end = inputs.t_end.powf(1.0 / (3.0 - k_cbm)) / (2.0 * (4.0 - k_cbm) * CCGS)
            * (1.0 + inputs.redshift);
    }

    //--------------------------------------------------------------------------
    // (4) We now have physically consistent values for start/end times and
    // start/end Gammas. This is sufficient to let us initialize all hydro
    // state arrays
    //--------------------------------------------------------------------------

    let t_fac = (inputs.t_end / inputs.t_start).powf(1.0 / (inputs.n_tm as f64 - 1.0));
    let mut t_curr = inputs.t_start;

    for i_tm in 0..inputs.n_tm {
        // Use GS2002 Eq. (A15) to compute radius and Lorentz factor
        let r_los =
            energy * t_curr / (1.0 + inputs.redshift) * (17.0 - 4.0 * k_cbm) * (4.0 - k_cbm)
                / (4.0 * PII * den_n * XMP * CCGS);

        let r_los = r_los.powf(1.0 / (4.0 - k_cbm));

        let gam_los = energy * (17.0 - 4.0 * k_cbm)
            / (8.0
                * PII
                * den_n
                * RM_PROT
                * (2.0 * (4.0 - k_cbm) * CCGS * t_curr / (1.0 + inputs.redshift))
                    .powf(3.0 - k_cbm));

        let gam_los = gam_los.powf(0.5 / (4.0 - k_cbm));

        // Use Eq. (26) of BM76 to compute t in engine rest frame
        let t_los = r_los / CCGS * (1.0 + 1.0 / (2.0 * (4.0 - k_cbm) * gam_los.powi(2)));

        // Place the computed values in the arrays
        globals.t_obs_array[i_tm] = t_curr;
        globals.r_los_array[i_tm] = r_los;
        globals.gam_los_array[i_tm] = gam_los;
        globals.t_los_array[i_tm] = t_los;

        // Increase t_curr for the next time through
        t_curr *= t_fac;
    }

    //--------------------------------------------------------------------------
    // Hydro arrays initialized
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Set up arrays to be used throughout the calculation:
    // (1) Step size for radative transfer ODE
    // (2) List of x values to use for integration in Eq. (A24)
    // (3) Synchrotron x and F(x) (different "x" than in above line)
    // (4) Plasma-frame electron momenta and associated Lorentz factors
    // (5) Photon energies and useful inverses
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // (1) Step sizes for radiative transfer
    //--------------------------------------------------------------------------

    let (em_steps, em_size_array) = set_em_step_sizes();
    globals.em_steps = em_steps;
    globals.em_size_array = em_size_array;

    //--------------------------------------------------------------------------
    // (2) x values at which we'll do the radiative transfer
    //--------------------------------------------------------------------------

    let (num_x, x_array1d) = set_x_array_1d();
    globals.num_x = num_x;

    //--------------------------------------------------------------------------
    // (3) Synchrotron x and F(x), which will save many, many calls to the
    // Bessel function during the running of the program
    //--------------------------------------------------------------------------

    let (syn_f, syn_x) = calc_syn_fx(N_PH);
    globals.syn_f = syn_f;
    globals.syn_x = syn_x;
    globals.o_o_sone = 1.0 / globals.syn_x[1];

    //--------------------------------------------------------------------------
    // (4) Array of electron momenta in the plasma frame, and Lorentz factors
    //   associated with same
    //--------------------------------------------------------------------------

    globals.p_pf_cgs = vec![0.0; N_DP + 1];
    globals.elec_gam_array = vec![0.0; N_DP + 1];

    globals.p_pf_cgs[0] = 1.0e-99;
    globals.elec_gam_array[0] = 1.0;
    globals.num_dp = 1;
    globals.p_pf_cgs[1] = P_ARRAY_MIN * XMP * CCGS;
    globals.elec_gam_array[1] = (1.0 + (globals.p_pf_cgs[1] / (XME * CCGS)).powi(2)).sqrt();

    for l_mom in 2..=N_DP {
        globals.num_dp += 1;

        globals.p_pf_cgs[l_mom] = globals.p_pf_cgs[l_mom - 1] * bin_factor_mom();
        globals.elec_gam_array[l_mom] =
            (1.0 + (globals.p_pf_cgs[l_mom] / (XME * CCGS)).powi(2)).sqrt();

        if globals.p_pf_cgs[l_mom] >= P_ARRAY_MAX * XMP * CCGS {
            break;
        }
    }

    // It is possible that we hit num_dp = n_dp without having hit p_pf = p_max,
    //   if p_array_max/p_array_min is large and bin_factor_mom is small.
    //   Should that happen, throw a warning
    if globals.p_pf_cgs[globals.num_dp] < P_ARRAY_MAX * XMP * CCGS {
        if inputs.do_progress_log {
            eprintln!("WARNING: need to expand momentum/gamma arrays");
            eprintln!(
                "      p_array_min: {:12.3e}; p_array_max: {:12.3e}; bins_per_decade_mom: {}",
                P_ARRAY_MIN, P_ARRAY_MAX, BINS_PER_DECADE_MOM,
            );
        }

        writeln!(out_writer, "WARNING: need to expand momentum/gamma arrays")?;
        writeln!(
            out_writer,
            "      p_array_min: {:12.3e}; p_array_max: {:12.3e}; bins_per_decade_mom: {}",
            P_ARRAY_MIN, P_ARRAY_MAX, BINS_PER_DECADE_MOM,
        )?;
        writeln!(out_writer)?;
    }

    //--------------------------------------------------------------------------
    // (5) Photon energies and useful inverses
    //--------------------------------------------------------------------------

    let phot_en_min_mev = 1.0e-17;
    let phot_en_max_mev = 1.0e9;
    globals.num_ph =
        1 + (f64::log10(phot_en_max_mev / phot_en_min_mev) * BINS_PER_DECADE_PHOT as f64) as usize;

    if globals.num_ph > N_PH {
        return Err(SaaError::PhotonArrayTooSmall {
            num_ph: globals.num_ph,
            n_ph: N_PH,
        });
    }

    // Fill arrays of photon energies to be used in radiative transfer calc
    globals.phot_en_cgs = vec![0.0; globals.num_ph];
    globals.o_o_nu_sq = vec![0.0; globals.num_ph];
    globals.o_o_alpha_array = vec![0.0; globals.num_ph];

    globals.phot_en_cgs[0] = phot_en_min_mev * MEV_TO_ERG;
    globals.o_o_nu_sq[0] = 1.0 / (globals.phot_en_cgs[0] / XH).powi(2);
    globals.o_o_alpha_array[0] = RM_ELEC / globals.phot_en_cgs[0];

    for m_ph in 1..globals.num_ph {
        globals.phot_en_cgs[m_ph] = globals.phot_en_cgs[m_ph - 1] * bin_factor_phot();
        globals.o_o_nu_sq[m_ph] = 1.0 / (globals.phot_en_cgs[m_ph] / XH).powi(2);
        globals.o_o_alpha_array[m_ph] = RM_ELEC / globals.phot_en_cgs[m_ph];
    }

    // AI.Panther change
    writeln!(sed_writer, "List of log10(photon energies), units of eV")?;
    for m_ph in 0..globals.num_ph {
        writeln!(
            sed_writer,
            "{:10.3e}",
            (globals.phot_en_cgs[m_ph] * ERG_TO_EV).log10()
        )?;
    }
    writeln!(sed_writer)?;

    //--------------------------------------------------------------------------
    // Other arrays initialized
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Handle cosmological effects: redshift and absorption by the CMB
    //--------------------------------------------------------------------------

    // Convert photon energies observed at Earth to photon energies in the rest
    //   frame of the central engine
    let phot_en_source: Vec<f64> = globals
        .phot_en_cgs
        .iter()
        .map(|phot_en| phot_en * (1.0 + inputs.redshift))
        .collect();

    // Compute the absorption due to the extragalactic background light
    let mut _transmittance = vec![1.0; globals.num_ph];
    let path = format!("{}/opacity-table.dat", env!("CARGO_MANIFEST_DIR"));
    let opacity_table = OpacityTable::from_file(&path).expect("failed to load opacity-table.dat");
    if inputs.do_absorb {
        calculate_absorptions(
            &opacity_table,
            globals.num_ph,
            inputs.redshift,
            &phot_en_source,
            &mut _transmittance,
        );
    }

    //--------------------------------------------------------------------------
    // Cosmology considered
    //--------------------------------------------------------------------------

    // Timing call
    //--------------------------------------------------------------------------
    let wall_time = start.elapsed().as_secs_f64();
    if inputs.do_progress_log {
        println!("Grid arrays initialized.  time = {:8.2e} sec", wall_time);
    }
    writeln!(
        out_writer,
        "Grid arrays initialized.  time = {:8.2e} sec",
        wall_time
    )?;
    writeln!(out_writer)?;

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //
    // Main computational loop.  At each time step, integrate flux over solid
    //   angle using GS2002's Eq. (A24).  For each line of sight, solve the
    //   radiative transfer equation, Eq. (A23), to find I_nu at that position.
    //
    // Loop proceeds in two parts.  First part finds extent of GPS99 "egg".
    //   This is a separate loop to aid in plotting the egg.  The second part
    //   uses the limits found previously to perform the radiative transfer
    //   calculation.
    //
    // Notation explanation:
    //   y: R_emis/R_axis       Ratio of shock radius at time of emission &
    //                            shock radius on axis (must be positive)
    //   x: R_perp/R_perp_max   Ratio of perpendicular location and max extent
    //                            of GPS99 egg (non-negative)
    //   s: r_par /R_axis       Position along line-of-sight integration
    //                            (may be negative).  Capital R for locations at
    //                            shock; lowercase r for interior points
    //
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    let i_tm = 39usize;

    for k_int in (0..=180).step_by(5) {
        // Set psi_fl
        let psi_fl = k_int as f64 * DEGTRD;

        //-----------------------------------------------------------------------
        // Find mu where R_perp_max occurs.  Do this in two stages for maximum
        //   robustness in the face of a sharply varying function.  First use
        //   Brent's method to get close to the minimum (approx 1 part in 10^6,
        //   which is roughly as close as this method can get).  Then use that
        //   as the initial position for Newton's method to find the root of the
        //   derivative to high precision.
        //-----------------------------------------------------------------------
        let mut i_params = [0i32; 7];
        let mut r_params = [0.0f64; 7];

        i_params[4] = 2; // off-axis
        r_params[0] = inputs.theta_obs;
        r_params[1] = psi_fl;
        r_params[2] = (17.0 - 4.0 * k_cbm) / (8.0 * PII * den_n * RM_PROT);
        r_params[3] = CCGS * globals.t_obs_array[i_tm] / (1.0 + inputs.redshift);
        r_params[4] = globals.r_los_array[i_tm];
        r_params[5] = 0.0;

        let mu_close = brents_minimum(
            |mu_input, i_params, r_params| fixed_rperp(mu_input, i_params, r_params, &inputs),
            &i_params,
            &r_params,
            0.5,
            0.9999,
            0.999999999,
        )?
        .x_min;

        let mu_max = newtons_method(
            |mu_input, i_params, r_params| d_rperp_dmu(mu_input, i_params, r_params, &inputs),
            &i_params,
            &r_params,
            mu_close,
            0.5,
            1.0,
        )?;

        // Convert mu_max into R_perp_max; remove the minus sign added in the
        //   function call
        r_params[5] = 0.0;
        let r_perp_max = -fixed_rperp(mu_max, &i_params, &r_params, &inputs);

        //------------------------------------------------------------------------
        let mut last_j_x = 0usize;

        for j_x in 0..globals.num_x {
            last_j_x = j_x + 1;

            if inputs.debug_mode {
                println!("{} {} {:e}", k_int, j_x + 1, mu_max);
            }

            // Choose a value between 0 and R_perp_max and find the two mu locations
            //   where R_perp = that value.  If R_perp = 0 then we can skip this
            //   step as it is possible to simply write down values for s_far and
            //   s_near
            let (r_perp, mu_near, mu_far) = if j_x == 0 {
                (0.0, 1.0, 0.0)
            } else {
                let r_perp = x_array1d[j_x] * r_perp_max;

                i_params = [0i32; 7];
                r_params[0] = inputs.theta_obs;
                r_params[1] = psi_fl;
                r_params[2] = (17.0 - 4.0 * k_cbm) / (8.0 * PII * den_n * RM_PROT);
                r_params[3] = CCGS * globals.t_obs_array[i_tm] / (1.0 + inputs.redshift);
                r_params[4] = globals.r_los_array[i_tm];
                r_params[5] = r_perp;

                let mu_close = brents_minimum(
                    |mu_input, i_params, r_params| {
                        fixed_rperp(mu_input, i_params, r_params, &inputs)
                    },
                    &i_params,
                    &r_params,
                    mu_max,
                    0.5 * (1.0 + mu_max),
                    1.0 - 1.0e-9,
                )?
                .x_min;

                let mu_near = newtons_method(
                    |mu_input, i_params, r_params| {
                        fixed_rperp(mu_input, i_params, r_params, &inputs)
                    },
                    &i_params,
                    &r_params,
                    mu_close,
                    mu_max,
                    1.0,
                )?;

                let mu_close = brents_minimum(
                    |mu_input, i_params, r_params| {
                        fixed_rperp(mu_input, i_params, r_params, &inputs)
                    },
                    &i_params,
                    &r_params,
                    -1.0 + 1.0e-9,
                    0.5 * (mu_max - 1.0),
                    mu_max,
                )?
                .x_min;

                let mu_far = newtons_method(
                    |mu_input, i_params, r_params| {
                        fixed_rperp(mu_input, i_params, r_params, &inputs)
                    },
                    &i_params,
                    &r_params,
                    mu_close,
                    -1.0,
                    mu_max,
                )?;

                (r_perp, mu_near, mu_far)
            };

            let (mut s_far, mut s_near) = if j_x == 0 {
                (
                    1.0e-4 * globals.r_los_array[i_tm],
                    globals.r_los_array[i_tm],
                )
            } else {
                i_params = [0i32; 7];
                i_params[4] = 2; // off-axis
                r_params[0] = inputs.theta_obs;
                r_params[1] = psi_fl;
                r_params[2] = (17.0 - 4.0 * k_cbm) / (8.0 * PII * den_n * RM_PROT);
                r_params[3] = CCGS * globals.t_obs_array[i_tm] / (1.0 + inputs.redshift);

                let r_far = r_of_theta(mu_far, &i_params, &r_params, &inputs);
                let r_near = r_of_theta(mu_near, &i_params, &r_params, &inputs);

                // DEBUGLINE
                if r_far <= r_perp {
                    println!("R_far < R_perp");
                    return Err(SaaError::InvalidOffAxisGeometry { r_far, r_perp });
                }

                (
                    (r_far.powi(2) - r_perp.powi(2)).sqrt(),
                    (r_near.powi(2) - r_perp.powi(2)).sqrt(),
                )
            };

            // Initialize input and output values for specific intensity
            let i_nu_in = vec![1.0e-99; globals.num_ph];

            // Set the parameters to be passed to rad_transfer_diffeq
            i_params[0] = i_tm as i32;
            i_params[1] = j_x as i32;
            i_params[2] = 0;
            i_params[3] = 0;
            i_params[4] = 2; // off-axis
            i_params[5] = 0;
            i_params[6] = 0;

            r_params[0] = r_perp / globals.r_los_array[i_tm];
            r_params[1] = r_perp_max / globals.r_los_array[i_tm];
            r_params[2] = psi_fl;
            r_params[3] = 0.0;
            r_params[4] = 0.0;
            r_params[5] = 0.0;
            r_params[6] = 0.0;

            // Rescale all values according to R_LOS
            s_far /= globals.r_los_array[i_tm];
            s_near /= globals.r_los_array[i_tm];

            // Evaluate the ODE
            let i_nu_out = bw_euler_stepper(
                |s, n_phot, jnu, alpha, r_params, i_params| {
                    rad_transfer_diffeq(RadTransferInput {
                        s,
                        n_phot,
                        jnu,
                        alpha,
                        r_params,
                        i_params,
                        inputs: &inputs,
                        globals: &globals,
                    })
                    .expect("rad_transfer_diffeq failed");
                },
                s_far,
                s_near,
                &i_nu_in,
                &r_params,
                &i_params,
                &globals,
            );

            let x_coord = mu_near.acos() * RADTDG * psi_fl.sin();
            let y_coord = -mu_near.acos() * RADTDG * psi_fl.cos() + inputs.theta_obs * RADTDG;

            // Write I_nu at certain energy ranges to different files
            writeln!(
                image_1046,
                "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
                k_int,
                j_x + 1,
                x_coord,
                y_coord,
                i_nu_out[45].log10()
            )?;
            writeln!(
                image_1070,
                "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
                k_int,
                j_x + 1,
                x_coord,
                y_coord,
                i_nu_out[69].log10()
            )?;
            writeln!(
                image_1089,
                "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
                k_int,
                j_x + 1,
                x_coord,
                y_coord,
                i_nu_out[88].log10()
            )?;
            writeln!(
                image_1113,
                "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
                k_int,
                j_x + 1,
                x_coord,
                y_coord,
                i_nu_out[112].log10()
            )?;
            writeln!(
                image_1161,
                "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
                k_int,
                j_x + 1,
                x_coord,
                y_coord,
                i_nu_out[160].log10()
            )?;
            writeln!(
                image_1185,
                "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
                k_int,
                j_x + 1,
                x_coord,
                y_coord,
                i_nu_out[184].log10()
            )?;
        }

        // Write I_nu at certain energy ranges to different files
        let x_coord = mu_max.acos() * RADTDG * psi_fl.sin();
        let y_coord = -mu_max.acos() * RADTDG * psi_fl.cos() + inputs.theta_obs * RADTDG;
        let sentinel_j_x = last_j_x + 1;

        writeln!(
            image_1046,
            "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
            k_int, sentinel_j_x, x_coord, y_coord, -99.0
        )?;
        writeln!(
            image_1070,
            "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
            k_int, sentinel_j_x, x_coord, y_coord, -99.0
        )?;
        writeln!(
            image_1089,
            "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
            k_int, sentinel_j_x, x_coord, y_coord, -99.0
        )?;
        writeln!(
            image_1113,
            "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
            k_int, sentinel_j_x, x_coord, y_coord, -99.0
        )?;
        writeln!(
            image_1161,
            "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
            k_int, sentinel_j_x, x_coord, y_coord, -99.0
        )?;
        writeln!(
            image_1185,
            "{:4}{:4}{:12.3e}{:12.3e}{:12.3e}",
            k_int, sentinel_j_x, x_coord, y_coord, -99.0
        )?;
    }

    // Blank lines
    writeln!(image_1046)?;
    writeln!(image_1070)?;
    writeln!(image_1089)?;
    writeln!(image_1113)?;
    writeln!(image_1161)?;
    writeln!(image_1185)?;

    // Energy information
    //                            Photon en. [eV]                 Photon freq. [Hz]
    writeln!(
        image_1046,
        "{:4}{:12.3e}{:12.3e}",
        "",
        (globals.phot_en_cgs[45] * ERG_TO_EV).log10(),
        (globals.phot_en_cgs[45] / XH).log10()
    )?;
    writeln!(
        image_1070,
        "{:4}{:12.3e}{:12.3e}",
        "",
        (globals.phot_en_cgs[69] * ERG_TO_EV).log10(),
        (globals.phot_en_cgs[69] / XH).log10()
    )?;
    writeln!(
        image_1089,
        "{:4}{:12.3e}{:12.3e}",
        "",
        (globals.phot_en_cgs[88] * ERG_TO_EV).log10(),
        (globals.phot_en_cgs[88] / XH).log10()
    )?;
    writeln!(
        image_1113,
        "{:4}{:12.3e}{:12.3e}",
        "",
        (globals.phot_en_cgs[112] * ERG_TO_EV).log10(),
        (globals.phot_en_cgs[112] / XH).log10()
    )?;
    writeln!(
        image_1161,
        "{:4}{:12.3e}{:12.3e}",
        "",
        (globals.phot_en_cgs[160] * ERG_TO_EV).log10(),
        (globals.phot_en_cgs[160] / XH).log10()
    )?;
    writeln!(
        image_1185,
        "{:4}{:12.3e}{:12.3e}",
        "",
        (globals.phot_en_cgs[184] * ERG_TO_EV).log10(),
        (globals.phot_en_cgs[184] / XH).log10()
    )?;

    // Timing call
    //--------------------------------------------------------------------------
    let wall_time = start.elapsed().as_secs_f64();

    if inputs.do_progress_log {
        println!("Emission calculated!  time = {:8.2e} sec", wall_time);
    }

    writeln!(out_writer)?;
    writeln!(
        out_writer,
        "Emission calculated!  time = {:8.2e} sec",
        wall_time
    )?;

    // Deallocate arrays to prevent memory leaks
    //--------------------------------------------------------------------------
    // Rust drops owned vectors automatically here.

    // And close the output files
    //--------------------------------------------------------------------------
    out_writer.flush()?;
    sed_writer.flush()?;
    image_1046.flush()?;
    image_1070.flush()?;
    image_1089.flush()?;
    image_1113.flush()?;
    image_1161.flush()?;
    image_1185.flush()?;

    Ok(())
}
