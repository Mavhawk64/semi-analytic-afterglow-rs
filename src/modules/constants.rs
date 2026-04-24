/// 1/3 for avoiding divisions
pub const THIRD: f64 = 1.0 / 3.0;

/// π constant
pub const PII: f64 = std::f64::consts::PI;

/// Degrees to radians
pub const DEG_TO_RAD: f64 = PII / 180.0;

/// Radians to degrees
pub const RAD_TO_DEG: f64 = 1.0 / DEG_TO_RAD;

/// Proton mass in grams
/// m_p = 1.6726218e-24 g
pub const XMP: f64 = 1.6726218e-24;

/// Electron mass
/// m_e = 9.10938291e-28 g
pub const XME: f64 = 9.10938291e-28;

/// Speed of light (cm/s)
/// c = 2.99792458e10 cm/s
pub const CCGS: f64 = 2.99792458e10;

/// Proton charge in cgs units (Electrostatic Units)
/// q = 4.80320451e-10 statcoulombs
pub const QCGS: f64 = 4.80320451e-10;

/// Thomson cross section (cm^-2)
/// σ_T = (8π/3) * (q^4 / (m_e^2 * c^4))
pub const SIG_T: f64 = 6.652458e-25;

/// Boltzmann's constant (cgs)
/// k_B = 1.3806488e-16 erg/K
pub const XKB: f64 = 1.3806488e-16;

/// Planck constant (erg·s)
/// h = 6.626070e-27 erg·s
pub const XH: f64 = 6.626070e-27;

/// Reduced Planck constant ħ = h / (2π)
pub const XH_BAR: f64 = XH / (2.0 * PII);

/// Proton rest mass energy [erg]
/// E_p = m_p * c^2
pub const RM_PROT: f64 = XMP * (CCGS * CCGS);

/// Electron rest mass energy [erg]
/// E_e = m_e * c^2
pub const RM_ELEC: f64 = XME * (CCGS * CCGS);

/// 1 / (m_e * c)
pub const O_O_MEC: f64 = 1.0 / (XME * CCGS);

/// Reciprocal electron rest energy [1/erg]
pub const O_O_RME: f64 = 1.0 / RM_ELEC;

/// Conversion: parsec to cm
pub const PC_TO_CM: f64 = 3.084e18;

/// Conversion: ergs to eV
pub const ERG_TO_EV: f64 = 6.242e11;

/// Conversion: ergs to keV
pub const ERG_TO_KEV: f64 = 6.242e8;

/// Conversion: ergs to MeV
pub const ERG_TO_MEV: f64 = 6.242e5;

/// Conversion: eV to ergs
pub const EV_TO_ERG: f64 = 1.602e-12;

/// Conversion: keV to ergs
pub const KEV_TO_ERG: f64 = 1.602e-9;

/// Conversion: MeV to ergs
pub const MEV_TO_ERG: f64 = 1.602e-6;
