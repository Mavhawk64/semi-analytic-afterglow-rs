use crate::modules::parameters::{BINS_PER_DECADE_MOM, BINS_PER_DECADE_PHOT, N_DP, N_PH};

/// Min momenta in initial dN/dp arrays, in units of m_p*c
pub const P_ARRAY_MIN: f64 = 1.0e-4;

/// Max momenta in initial dN/dp arrays, in units of m_p*c
pub const P_ARRAY_MAX: f64 = 1.0e12;

/// NOTE: These have to be functions instead of constants because
/// they depend on other functions, and Rust doesn't allow for
/// that in `const` definitions.
///
/// Number of bins the initial dN/dp arrays need
pub fn bin_factor_mom() -> f64 {
    10.0f64.powf(1.0 / BINS_PER_DECADE_MOM as f64)
}

/// Number of photon energies
pub fn bin_factor_phot() -> f64 {
    10.0f64.powf(1.0 / BINS_PER_DECADE_PHOT as f64)
}

/// Runtime global state translated from `globals.f90`.
#[derive(Debug, Clone)]
pub struct Globals {
    /// Arrays holding bulk hydrodynamic quantities
    pub t_obs_array: Vec<f64>,
    pub t_los_array: Vec<f64>,
    pub r_los_array: Vec<f64>,
    pub gam_los_array: Vec<f64>,

    /// Momentum at which to switch between thermal distribution and power-law tail
    pub num_dp: usize,

    /// Arrays to hold momenta and electron distributions for passing to subroutines
    pub p_pf_cgs: Vec<f64>,
    pub elec_gam_array: Vec<f64>,

    /// Number of photon energies
    pub num_ph: usize,

    pub o_o_sone: f64,

    pub syn_f: Vec<f64>,
    pub syn_x: Vec<f64>,

    /// Array of photon energies in ergs
    pub phot_en_cgs: Vec<f64>,

    pub o_o_nu_sq: Vec<f64>,
    pub o_o_alpha_array: Vec<f64>,

    /// Number of steps to take in radiative transfer
    pub em_steps: usize,

    /// Size of the steps to take in radiative transfer
    pub em_size_array: Vec<f64>,

    /// Number of lines of sight to use in angular integration
    pub num_x: usize,

    /// Locations of lines of sight between 0 and 1
    pub x_array: [[f64; 15]; 3],

    pub x_array_1d: [f64; 50],
}

/// The Fortran code kind of assumes things are already initialized or
/// allocated before they're used. In Rust, we have to be explicit, so
/// everything gets allocated and zeroed out here to keep things predictable.
/// Hence, the `Default` implementation.
///
/// Also, any arrays that were `dimension(0:n)` in Fortran are created
/// with length `n + 1` so the indexing lines up the same way.
impl Default for Globals {
    fn default() -> Self {
        Self {
            t_obs_array: Vec::new(),
            t_los_array: Vec::new(),
            r_los_array: Vec::new(),
            gam_los_array: Vec::new(),

            num_dp: 0,
            p_pf_cgs: vec![0.0; N_DP + 1],
            elec_gam_array: vec![0.0; N_DP + 1],

            num_ph: 0,
            o_o_sone: 0.0,

            syn_f: vec![0.0; N_PH + 1],
            syn_x: vec![0.0; N_PH + 1],

            phot_en_cgs: vec![0.0; N_PH],
            o_o_nu_sq: vec![0.0; N_PH],
            o_o_alpha_array: vec![0.0; N_PH],

            em_steps: 0,
            em_size_array: vec![0.0; N_PH],

            num_x: 0,
            x_array: [[0.0; 15]; 3],
            x_array_1d: [0.0; 50],
        }
    }
}
