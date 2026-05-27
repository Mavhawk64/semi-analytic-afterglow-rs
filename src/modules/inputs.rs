/// Runtime input parameters for the semi-analytic afterglow model.
#[derive(Debug, Clone, Default)]
pub struct Inputs {
    /// Number of time steps to use
    pub n_tm: usize,

    /// Is circumburst medium wind-like (k=2) or constant-density (k=0)?
    pub k_cbm: i32,

    /// Number density of CBM when constant
    pub n_cbm: f64,

    /// Density parameter of CBM when windlike
    pub a_star: f64,

    /// Compute EBL absorption between source & Earth?
    pub do_absorb: bool,

    /// Compute synchrotron self-absorption?
    pub do_ssa: bool,

    /// Compute SSC emission? (takes a lot of time)
    pub do_ssc: bool,

    /// Place limit on electron momentum?
    pub do_pmax: bool,

    /// Place limit (=eta) on shock Lorentz factor?
    ///
    /// This flag applies only in rad_transfer, not to Gam_start, which is always <= eta.
    pub do_gammax: bool,

    /// Is electron distribution just a power law, or pl+thermal?
    /// true -> pure power law; false -> pl+thermal
    pub pl_only: bool,

    /// Normalization of non-thermal tail for mixed pl+thermal distributions
    pub nt_norm: f64,

    /// If electron distribution is a power law, this is its spectral index
    pub p_spec: f64,

    /// Isotropic energy equivalent along jet axis
    pub e_iso: f64,

    /// E/m for ejecta, and also peak Lorentz factor of jet
    pub eta: f64,

    /// Pre-decay B field strength
    pub eps_b0: f64,

    /// Default electron energy fraction
    pub epse0: f64,

    /// Time at which observations start, in the observer (including redshift) frame
    pub t_start: f64,

    /// Jet Lorentz factor at t_start
    pub gam_start: f64,

    /// Time at which observations end, also in the observer frame
    pub t_end: f64,

    /// Jet Lorentz factor at t_end
    pub gam_end: f64,

    /// Degrees between jet axis and observer line of sight
    pub observer_los_deg: f64,

    /// Comoving distance (cm) to GRB
    pub dist_lum: f64,

    /// Associated redshift
    pub redshift: f64,

    pub do_progress_log: bool,

    /// Jet structure, or spherically symmetric?
    /// true -> structured jet; false -> spherical symmetry
    pub do_jet_structure: bool,

    /// Degrees between jet axis and observer line of sight
    pub theta_obs: f64,

    /// Outermost angle of jet core
    pub theta_j: f64,

    /// Innermost angle of jet cocoon
    pub theta_c: f64,

    /// Energy decay within jet-cocoon region
    pub jet_delta: f64,

    /// Energy decay within cocoon region
    pub jet_f: f64,

    /// Lorentz factor decay outside jet core
    pub jet_p: f64,

    /// Whether we should turn on debug printing
    pub debug_mode: bool,
}
