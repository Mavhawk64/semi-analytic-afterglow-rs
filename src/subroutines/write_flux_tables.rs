use std::io::{Result as IoResult, Write};

use crate::modules::constants::{ERG_TO_EV, XH};
use crate::modules::globals::Globals;

/// Write the flux data to a SED file.
///
/// # Arguments
///
/// * `writer` - output writer to write to
/// * `i_tm` - index of current time stamp
/// * `f_nu` - array of photon flux
/// * `nu_fnu` - array of energy flux
/// * `globals` - runtime global variables containing `t_obs_array`,
///   `num_ph`, and `phot_en_cgs`
pub fn write_flux_tables<W: Write>(
    writer: &mut W,
    i_tm: usize,
    f_nu: &[f64],
    nu_fnu: &[f64],
    globals: &Globals,
) -> IoResult<()> {
    writeln!(
        writer,
        "i_tm = {:3};  t_obs = {:12.3e}",
        i_tm, globals.t_obs_array[i_tm]
    )?;

    for m_ph in 0..globals.num_ph {
        writeln!(
            writer,
            "{:4}{:4}{:14.4e}{:14.4e}{:13.4e}{:14.4e}",
            i_tm,
            m_ph + 1,
            (globals.phot_en_cgs[m_ph] * ERG_TO_EV).log10(),
            (globals.phot_en_cgs[m_ph] / XH).log10(),
            nu_fnu[m_ph].log10(),
            f_nu[m_ph].log10(),
        )?;
    }

    // Write a line to explain the columns
    writeln!(
        writer,
        "{:12}log10(eV)     log10(Hz)    log10(nuFnu)  log10(Jy)",
        ""
    )?;

    // Skip a line before next time step
    writeln!(writer)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_globals() -> Globals {
        Globals {
            num_ph: 2,
            t_obs_array: vec![10.0],
            phot_en_cgs: vec![1.0e-12, 2.0e-12],
            ..Default::default()
        }
    }

    #[test]
    fn writes_header_with_time_index_and_observer_time() {
        let globals = base_globals();
        let f_nu = vec![1.0e-3, 2.0e-3];
        let nu_fnu = vec![1.0e-10, 2.0e-10];

        let mut buffer = Vec::new();

        write_flux_tables(&mut buffer, 0, &f_nu, &nu_fnu, &globals)
            .expect("write_flux_tables should succeed");

        let output = String::from_utf8(buffer).expect("output should be valid UTF-8");

        assert!(output.contains("i_tm =   0;  t_obs ="));
    }

    #[test]
    fn writes_one_row_per_photon_bin() {
        let globals = base_globals();
        let f_nu = vec![1.0e-3, 2.0e-3];
        let nu_fnu = vec![1.0e-10, 2.0e-10];

        let mut buffer = Vec::new();

        write_flux_tables(&mut buffer, 0, &f_nu, &nu_fnu, &globals)
            .expect("write_flux_tables should succeed");

        let output = String::from_utf8(buffer).expect("output should be valid UTF-8");

        assert!(output.contains("   0   1"));
        assert!(output.contains("   0   2"));
    }

    #[test]
    fn writes_column_explanation_line() {
        let globals = base_globals();
        let f_nu = vec![1.0e-3, 2.0e-3];
        let nu_fnu = vec![1.0e-10, 2.0e-10];

        let mut buffer = Vec::new();

        write_flux_tables(&mut buffer, 0, &f_nu, &nu_fnu, &globals)
            .expect("write_flux_tables should succeed");

        let output = String::from_utf8(buffer).expect("output should be valid UTF-8");

        assert!(output.contains("log10(eV)"));
        assert!(output.contains("log10(Hz)"));
        assert!(output.contains("log10(nuFnu)"));
        assert!(output.contains("log10(Jy)"));
    }
}
