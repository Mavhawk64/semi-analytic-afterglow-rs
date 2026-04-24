use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::modules::constants::MEV_TO_ERG;

/// Here are some constants added that weren't in the original Fortran code,
/// but they are helpful (rather than using magic numbers in the code).
const OPAC_N_PHOT: usize = 50;
const OPAC_N_Z: usize = 201;
const OPAC_TAU_MAX: f64 = 35.0;

#[derive(Debug, Clone)]
pub struct OpacityTable {
    /// Array of photon energies.
    pub opac_phot_table: [f64; OPAC_N_PHOT],

    /// Table for storing opacity data from Franceschini+ 2008 (2008A&A...487..837F).
    /// This should only be accessed by `get_opacity()` and `init_opacity_table()`
    ///
    /// Fortran shape was `dimension(50,0:200)`.
    /// Rust stores this as `[redshift_row][photon_index]`.
    pub opac_table: [[f64; OPAC_N_PHOT]; OPAC_N_Z],
}

impl OpacityTable {
    /// Read opacity data from an external file.
    ///
    /// The file must have:
    /// - first line: comment
    /// - second line: 50 photon energies
    /// - third line: comment
    /// - next 201 lines: opacity rows
    pub fn init_opacity_table(filename: &str) -> std::io::Result<Self> {
        let file = File::open(filename)?;
        let mut lines = BufReader::new(file).lines();

        // The first line holds just comments
        lines.next().transpose()?;

        // Read the photon energies
        let phot_line = lines
            .next()
            .transpose()?
            .ok_or_else(|| std::io::Error::other("missing photon energy row"))?;

        let phot_values = parse_f64_row(&phot_line)?;

        if phot_values.len() != OPAC_N_PHOT {
            return Err(std::io::Error::other(format!(
                "expected {OPAC_N_PHOT} photon energies, got {}",
                phot_values.len()
            )));
        }

        let mut opac_phot_table = [0.0; OPAC_N_PHOT];
        opac_phot_table.copy_from_slice(&phot_values);

        // Next line is also a comment
        lines.next().transpose()?;

        let mut opac_table = [[0.0; OPAC_N_PHOT]; OPAC_N_Z];

        // For the next 201 lines, fill in a row of opac_table
        for (row, opac_row) in opac_table.iter_mut().enumerate() {
            let line = lines
                .next()
                .transpose()?
                .ok_or_else(|| std::io::Error::other(format!("missing opacity row {row}")))?;

            let values = parse_f64_row(&line)?;

            if values.len() != OPAC_N_PHOT {
                return Err(std::io::Error::other(format!(
                    "expected {OPAC_N_PHOT} values in opacity row {row}, got {}",
                    values.len()
                )));
            }

            opac_row.copy_from_slice(&values);
        }

        Ok(Self {
            opac_phot_table,
            opac_table,
        })
    }

    /// Uses the opacity tables from Franceschini+ 2008 (2008A&A...487..837F) to
    /// compute the transmittance of high-energy photons originating from a
    /// source at the specified redshift.
    pub fn get_opacity(&self, redshift: f64) -> ([f64; OPAC_N_PHOT], [f64; OPAC_N_PHOT]) {
        // Output arguments: opac_phot_e and opac_tau_e
        // Array of photon energies
        // ============================================================================
        //  Photon energy table; this is the same for all redshifts
        //  Values are in linear space, with units of MeV
        // ============================================================================
        let opac_phot_e = self.opac_phot_table;
        // Array of optical depths for those photon energies
        let mut opac_tau_e = [0.0; OPAC_N_PHOT];

        if redshift >= 2.0 {
            eprintln!("WARNING: in get_opacity, redshifts > 2 are treated as if z = 2.");

            for (i, tau) in opac_tau_e.iter_mut().enumerate() {
                // Place an upper limit on opacity since we'll be taking
                // e^-tau later
                // NOTE: This is to cut off double-precision / floating-point
                // errors around 0.
                *tau = self.opac_table[200][i].min(OPAC_TAU_MAX);
            }
        } else {
            // Calculate the row/rows needed for finding the opacity
            let z_scaled = redshift * 100.0;
            let row_lo = z_scaled.floor() as usize;
            let row_hi = z_scaled.ceil() as usize;

            // This rust code differs from the current (2026-04-24) Fortran code
            // but it was a suggestion I made a while back to solve a possible
            // divide by zero issue.
            if row_lo == row_hi {
                for (i, tau) in opac_tau_e.iter_mut().enumerate() {
                    // In this case, we just set the opac_tau_e values to the
                    // corresponding row in the table, with an upper limit
                    *tau = self.opac_table[row_lo][i].min(OPAC_TAU_MAX);
                }
            } else {
                // Otherwise, we can proceed as normal
                let z_lo = row_lo as f64 * 1.0e-2;
                let z_hi = row_hi as f64 * 1.0e-2;

                // Linearly interpolate between opacities to find value for redshift
                for (i, tau) in opac_tau_e.iter_mut().enumerate() {
                    let slope =
                        (self.opac_table[row_hi][i] - self.opac_table[row_lo][i]) / (z_hi - z_lo);

                    *tau = self.opac_table[row_hi][i] - slope * (z_hi - redshift);

                    // Place an upper limit on opacity since we'll be taking e^-tau later
                    *tau = tau.min(OPAC_TAU_MAX);
                }
            }
        }

        for tau in &mut opac_tau_e {
            *tau = (-*tau).exp();
        }

        (opac_phot_e, opac_tau_e)
    }
}

/// TODO: add doc comment
pub fn calculate_absorptions(
    opacity_table: &OpacityTable,
    num_ph: usize, // Number of photon energies (Actual size of photon array)
    redshift: f64,
    phot_en_source: &[f64],
    transmittance: &mut [f64],
) {
    // NOTE:
    // In the original Fortran, arrays were defined with indices (0:n), and only
    // indices 1..n were filled by `get_opacity`. The 0th index was then manually
    // set to 0 as a boundary/sentinel value.
    //
    // In this Rust version, we use standard 0-based indexing (0..n-1), so the
    // data is already aligned correctly and no extra sentinel element is needed.
    let (opac_phot_e, opac_xmit) = opacity_table.get_opacity(redshift);

    for m_ph in 0..num_ph {
        // Small optimization from the fortran code
        // Don't find element in array everytime we need it
        let energy = phot_en_source[m_ph];

        if energy < opac_phot_e[0] * MEV_TO_ERG {
            transmittance[m_ph] = 1.0;
        } else if energy > opac_phot_e[OPAC_N_PHOT - 1] * MEV_TO_ERG {
            transmittance[m_ph] = opac_xmit[OPAC_N_PHOT - 1];
        } else {
            let mut m_lo = 0;
            let mut m_hi = 1;

            // Linear search since OPAC_N_PHOT is small
            // O(log n) ~ O(n) for n=50
            // This is also a rusty way of enumerating through the array
            // (same kind of process as python)
            for (m_si, &phot_e) in opac_phot_e.iter().enumerate().skip(1) {
                if energy < phot_e * MEV_TO_ERG {
                    m_lo = m_si - 1;
                    m_hi = m_si;
                    break;
                }
            }

            // Possible binary search implementation if needed:
            // match opac_phot_e
            //     .binary_search_by(|val| (val * MEV_TO_ERG).partial_cmp(&energy).unwrap())
            // {
            //     Ok(idx) => {
            //         m_lo = idx;
            //         m_hi = idx;
            //     }
            //     Err(idx) => {
            //         m_lo = idx - 1;
            //         m_hi = idx;
            //     }
            // }

            let delta_energy = (opac_phot_e[m_hi] - opac_phot_e[m_lo]) * MEV_TO_ERG;
            let slope = (opac_xmit[m_hi] - opac_xmit[m_lo]) / delta_energy;
            let delta = energy - opac_phot_e[m_lo] * MEV_TO_ERG;

            transmittance[m_ph] = opac_xmit[m_lo] + slope * delta;
        }
    }
}

/// Helper function to parse a line of whitespace-separated floats into a Vec<f64>.
/// (not in the original Fortran code)
fn parse_f64_row(line: &str) -> std::io::Result<Vec<f64>> {
    line.split_whitespace()
        .map(|s| {
            s.parse::<f64>()
                .map_err(|err| std::io::Error::other(format!("failed to parse float '{s}': {err}")))
        })
        .collect()
}

/// Unit tests for the opacity module.
#[cfg(test)]
mod tests {
    use super::*;

    fn test_table() -> OpacityTable {
        let mut opac_phot_table = [0.0; OPAC_N_PHOT];
        for (i, value) in opac_phot_table.iter_mut().enumerate() {
            *value = (i + 1) as f64;
        }

        let mut opac_table = [[0.0; OPAC_N_PHOT]; OPAC_N_Z];
        for (z_row, row) in opac_table.iter_mut().enumerate() {
            for (phot_i, value) in row.iter_mut().enumerate() {
                *value = z_row as f64 + phot_i as f64;
            }
        }

        OpacityTable {
            opac_phot_table,
            opac_table,
        }
    }

    #[test]
    fn get_opacity_uses_exact_row_when_redshift_lands_on_grid() {
        let table = test_table();

        let (_, tau) = table.get_opacity(0.50);

        let expected_raw_tau: f64 = 50.0;
        let expected = (-expected_raw_tau.min(OPAC_TAU_MAX)).exp();

        assert!((tau[0] - expected).abs() < 1e-14);
    }

    #[test]
    fn get_opacity_clamps_redshift_above_two() {
        let table = test_table();

        let (_, tau) = table.get_opacity(3.0);

        let expected_raw_tau: f64 = 200.0;
        let expected = (-expected_raw_tau.min(OPAC_TAU_MAX)).exp();

        assert!((tau[0] - expected).abs() < 1e-14);
    }

    #[test]
    fn calculate_absorptions_sets_low_energy_to_one() {
        let table = test_table();

        let phot_en_source = vec![0.5 * MEV_TO_ERG];
        let mut transmittance = vec![0.0];

        calculate_absorptions(&table, 1, 0.0, &phot_en_source, &mut transmittance);

        assert_eq!(transmittance[0], 1.0);
    }

    #[test]
    fn init_opacity_table_loads_repo_opacity_table() {
        let path = format!("{}/opacity-table.dat", env!("CARGO_MANIFEST_DIR"));

        let table =
            OpacityTable::init_opacity_table(&path).expect("failed to load opacity-table.dat");

        assert_eq!(table.opac_phot_table.len(), OPAC_N_PHOT);
        assert_eq!(table.opac_table.len(), OPAC_N_Z);

        assert!(table.opac_phot_table.iter().all(|x| x.is_finite()));
        assert!(table.opac_table.iter().flatten().all(|x| x.is_finite()));
    }
}
