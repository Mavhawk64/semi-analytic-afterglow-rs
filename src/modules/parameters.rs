/// Number of threads to use in OpenMP calculations (Fortran reference)
pub const OMP_NUMTH: usize = 4;

/// Max size of dN/dp array in momentum
pub const N_DP: usize = 2500;

/// Max size of photon array in energy
pub const N_PH: usize = 250;

/// Number of bins per decade in log space for SEDs
pub const BINS_PER_DECADE_PHOT: usize = 8;

/// Number of bins per decade in log space for dN/dp
pub const BINS_PER_DECADE_MOM: usize = 16;
