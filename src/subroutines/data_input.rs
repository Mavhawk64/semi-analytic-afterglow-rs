use crate::modules::inputs::Inputs;

#[derive(Debug, Clone, PartialEq)]
pub enum DataInputError {
    UnknownKeyword(String),
    MissingValue(String),
    InvalidValue { keyword: String, message: String },
    ConflictingInputs { first: String, second: String },
    ParseError { keyword: String, value: String },
}

#[derive(Debug, Clone, Default)]
struct PartialInputs {
    n_tm: Option<usize>,
    k_cbm: Option<i32>,
    do_absorb: Option<bool>,
    do_ssa: Option<bool>,
    do_ssc: Option<bool>,
    pl_only: Option<bool>,
    nt_norm: Option<f64>,
    p_spec: Option<f64>,
    e_iso: Option<f64>,
    n_cbm: Option<f64>,
    a_star: Option<f64>,
    eta: Option<f64>,
    eps_b0: Option<f64>,
    epse0: Option<f64>,
    t_start: Option<f64>,
    gam_start: Option<f64>,
    t_end: Option<f64>,
    gam_end: Option<f64>,
    dist_lum: Option<f64>,
    redshift: Option<f64>,
    do_pmax: Option<bool>,
    do_gammax: Option<bool>,
    theta_obs: Option<f64>,
    theta_j: Option<f64>,
    theta_c: Option<f64>,
    jet_delta: Option<f64>,
    jet_f: Option<f64>,
    jet_p: Option<f64>,
    do_jet_structure: Option<bool>,
    do_progress_log: Option<bool>,
    debug_mode: Option<bool>,
}

pub fn data_input_from_str(text: &str) -> Result<Inputs, DataInputError> {
    let mut partial = PartialInputs::default();

    for raw_line in text.lines() {
        let line = raw_line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut parts = line.split_whitespace();
        let keyword = parts
            .next()
            .ok_or_else(|| DataInputError::MissingValue("keyword".to_string()))?;

        if keyword == "ENDIN" {
            break;
        }

        match keyword {
            "ASTAR" => {
                let value = parse_f64(keyword, parts.next())?;
                require_nonnegative(keyword, value, "must be positive or 0.0 to ignore")?;
                partial.a_star = Some(value);
            }
            "D_LUM" => {
                let value = parse_f64(keyword, parts.next())?;
                require_nonnegative(keyword, value, "must be positive or 0.0 to ignore")?;
                partial.dist_lum = Some(value);
            }
            "DEBUG" => {
                partial.debug_mode = Some(parse_flag(keyword, parts.next())?);
            }
            "DOABS" => {
                partial.do_absorb = Some(parse_flag(keyword, parts.next())?);
            }
            "DOGMX" => {
                partial.do_gammax = Some(parse_flag(keyword, parts.next())?);
            }
            "DOPMX" => {
                partial.do_pmax = Some(parse_flag(keyword, parts.next())?);
            }
            "DOSSA" => {
                partial.do_ssa = Some(parse_flag(keyword, parts.next())?);
            }
            "DOSSC" => {
                partial.do_ssc = Some(parse_flag(keyword, parts.next())?);
            }
            "E_ISO" => {
                let value = parse_f64(keyword, parts.next())?;
                require_positive(keyword, value, "must be positive")?;
                partial.e_iso = Some(value);
            }
            "EPSB0" => {
                let value = parse_f64(keyword, parts.next())?;
                require_positive(keyword, value, "must be positive")?;
                partial.eps_b0 = Some(value);
            }
            "EPSE0" => {
                let value = parse_f64(keyword, parts.next())?;
                require_positive(keyword, value, "must be positive")?;
                partial.epse0 = Some(value);
            }
            "GAMPK" => {
                let value = parse_f64(keyword, parts.next())?;
                if value <= 1.0 {
                    return invalid(keyword, "eta must be > 1");
                }
                partial.eta = Some(value);
            }
            "GAMST" => {
                let value = parse_f64(keyword, parts.next())?;
                if value != 0.0 && value <= 1.0 {
                    return invalid(keyword, "Gam_start must be > 1, or 0.0 to ignore");
                }
                partial.gam_start = Some(value);
            }
            "GMEND" => {
                let value = parse_f64(keyword, parts.next())?;
                if value != 0.0 && value <= 1.0 {
                    return invalid(keyword, "Gam_end must be > 1, or 0.0 to ignore");
                }
                partial.gam_end = Some(value);
            }
            "JTPRM" => {
                let jet_delta = parse_f64(keyword, parts.next())?;
                let jet_f = parse_f64(keyword, parts.next())?;
                let jet_p = parse_f64(keyword, parts.next())?;

                require_positive(keyword, jet_delta, "jet_delta must be positive")?;
                require_positive(keyword, jet_f, "jet_f must be positive")?;
                require_positive(keyword, jet_p, "jet_p must be positive")?;

                partial.jet_delta = Some(jet_delta);
                partial.jet_f = Some(jet_f);
                partial.jet_p = Some(jet_p);
            }
            "JTTHT" => {
                let theta_j = parse_f64(keyword, parts.next())?;
                let theta_c = parse_f64(keyword, parts.next())?;

                require_positive(keyword, theta_j, "theta_j must be positive")?;

                if theta_c <= theta_j {
                    return invalid(keyword, "theta_c must be greater than theta_j");
                }

                partial.theta_j = Some(theta_j);
                partial.theta_c = Some(theta_c);
                partial.do_jet_structure = Some(theta_j > 0.0 && theta_c > 0.0);
            }
            "K_CBM" => {
                partial.k_cbm = Some(parse_i32(keyword, parts.next())?);
            }
            "N_CBM" => {
                partial.n_cbm = Some(parse_f64(keyword, parts.next())?);
            }
            "NTIME" => {
                let value = parse_usize(keyword, parts.next())?;

                if value == 0 {
                    return invalid(keyword, "n_tm must be positive");
                }

                partial.n_tm = Some(value);
            }
            "NTNRM" => {
                let value = parse_f64(keyword, parts.next())?;
                if value <= 0.0 || value > 1.0 {
                    return invalid(keyword, "nt_norm must be in (0, 1]");
                }
                partial.nt_norm = Some(value);
            }
            "NTSIG" => {
                let mut value = parse_f64(keyword, parts.next())?;

                if value <= 0.0 {
                    value = 0.0;
                } else if value <= 2.0 {
                    return invalid(keyword, "p_spec must be > 2, or <= 0 to ignore");
                }

                partial.p_spec = Some(value);
            }
            "PLONL" => {
                partial.pl_only = Some(parse_flag(keyword, parts.next())?);
            }
            "PRGRS" => {
                let value = parse_i32(keyword, parts.next())?;
                partial.do_progress_log = Some(value != 66);
            }
            "RDSHF" => {
                let value = parse_f64(keyword, parts.next())?;
                require_nonnegative(keyword, value, "redshift must be non-negative")?;
                partial.redshift = Some(value);
            }
            "T_END" => {
                let value = parse_f64(keyword, parts.next())?;
                require_nonnegative(keyword, value, "t_end must be non-negative")?;
                partial.t_end = Some(value);
            }
            "THOBS" => {
                let value = parse_f64(keyword, parts.next())?;
                require_nonnegative(keyword, value, "theta_obs may not be negative")?;
                partial.theta_obs = Some(value);
            }
            "TSTRT" => {
                let value = parse_f64(keyword, parts.next())?;
                require_nonnegative(keyword, value, "t_start must be non-negative")?;
                partial.t_start = Some(value);
            }
            _ => return Err(DataInputError::UnknownKeyword(keyword.to_string())),
        }
    }

    finalize_inputs(partial)
}

fn finalize_inputs(partial: PartialInputs) -> Result<Inputs, DataInputError> {
    if partial.redshift.unwrap_or(0.0) > 0.0 && partial.dist_lum.unwrap_or(0.0) > 0.0 {
        return Err(DataInputError::ConflictingInputs {
            first: "RDSHF".to_string(),
            second: "D_LUM".to_string(),
        });
    }

    if partial.gam_start.unwrap_or(0.0) > 1.0 && partial.t_start.unwrap_or(0.0) > 0.0 {
        return Err(DataInputError::ConflictingInputs {
            first: "GAMST".to_string(),
            second: "TSTRT".to_string(),
        });
    }

    if partial.gam_end.unwrap_or(0.0) > 1.0 && partial.t_end.unwrap_or(0.0) > 0.0 {
        return Err(DataInputError::ConflictingInputs {
            first: "GMEND".to_string(),
            second: "T_END".to_string(),
        });
    }

    let theta_obs = partial.theta_obs.unwrap_or(0.0);

    let theta_j = partial
        .theta_j
        .unwrap_or(if theta_obs == 0.0 { 3.0 } else { 0.0 });
    let theta_c = partial
        .theta_c
        .unwrap_or(if theta_obs == 0.0 { 20.0 } else { 0.0 });

    let do_jet_structure = partial
        .do_jet_structure
        .unwrap_or(theta_j > 0.0 && theta_c > theta_j && theta_obs == 0.0);

    Ok(Inputs {
        n_tm: partial.n_tm.unwrap_or(10),
        k_cbm: partial.k_cbm.unwrap_or(0),
        do_absorb: partial.do_absorb.unwrap_or(false),
        do_ssa: partial.do_ssa.unwrap_or(false),
        do_ssc: partial.do_ssc.unwrap_or(false),
        pl_only: partial.pl_only.unwrap_or(true),
        nt_norm: partial.nt_norm.unwrap_or(0.03),
        p_spec: partial.p_spec.unwrap_or(2.23),
        e_iso: partial.e_iso.unwrap_or(1.0e53),
        n_cbm: partial.n_cbm.unwrap_or(1.0),
        a_star: partial.a_star.unwrap_or(1.0e35),
        eta: partial.eta.unwrap_or(1000.0),
        eps_b0: partial.eps_b0.unwrap_or(0.01),
        epse0: partial.epse0.unwrap_or(0.3),
        t_start: partial.t_start.unwrap_or(0.0),
        gam_start: partial.gam_start.unwrap_or(400.0),
        t_end: partial.t_end.unwrap_or(0.0),
        gam_end: partial.gam_end.unwrap_or(2.5),
        dist_lum: partial.dist_lum.unwrap_or(1.0e28),
        redshift: partial.redshift.unwrap_or(0.0),
        do_pmax: partial.do_pmax.unwrap_or(false),
        do_gammax: partial.do_gammax.unwrap_or(false),
        theta_obs,
        theta_j,
        theta_c,
        jet_delta: partial.jet_delta.unwrap_or(1.8),
        jet_f: partial.jet_f.unwrap_or(1.5),
        jet_p: partial.jet_p.unwrap_or(2.4),
        do_jet_structure,
        do_progress_log: partial.do_progress_log.unwrap_or(true),
        debug_mode: partial.debug_mode.unwrap_or(false),
        ..Default::default()
    })
}

fn parse_flag(keyword: &str, value: Option<&str>) -> Result<bool, DataInputError> {
    Ok(parse_i32(keyword, value)? == 66)
}

fn parse_i32(keyword: &str, value: Option<&str>) -> Result<i32, DataInputError> {
    let value = value.ok_or_else(|| DataInputError::MissingValue(keyword.to_string()))?;

    value
        .parse::<i32>()
        .map_err(|_| DataInputError::ParseError {
            keyword: keyword.to_string(),
            value: value.to_string(),
        })
}

fn parse_f64(keyword: &str, value: Option<&str>) -> Result<f64, DataInputError> {
    let value = value.ok_or_else(|| DataInputError::MissingValue(keyword.to_string()))?;

    value
        .parse::<f64>()
        .map_err(|_| DataInputError::ParseError {
            keyword: keyword.to_string(),
            value: value.to_string(),
        })
}

fn parse_usize(keyword: &str, value: Option<&str>) -> Result<usize, DataInputError> {
    let value = value.ok_or_else(|| DataInputError::MissingValue(keyword.to_string()))?;

    value
        .parse::<usize>()
        .map_err(|_| DataInputError::ParseError {
            keyword: keyword.to_string(),
            value: value.to_string(),
        })
}

fn require_positive(keyword: &str, value: f64, message: &str) -> Result<(), DataInputError> {
    if value <= 0.0 {
        invalid(keyword, message)
    } else {
        Ok(())
    }
}

fn require_nonnegative(keyword: &str, value: f64, message: &str) -> Result<(), DataInputError> {
    if value < 0.0 {
        invalid(keyword, message)
    } else {
        Ok(())
    }
}

fn invalid<T>(keyword: &str, message: &str) -> Result<T, DataInputError> {
    Err(DataInputError::InvalidValue {
        keyword: keyword.to_string(),
        message: message.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn applies_defaults_when_only_endin_is_provided() {
        let inputs = data_input_from_str("ENDIN").expect("data_input should succeed");

        assert_eq!(inputs.n_tm, 10);
        assert_eq!(inputs.k_cbm, 0);
        assert_eq!(inputs.n_cbm, 1.0);
        assert_eq!(inputs.nt_norm, 0.03);
        assert_eq!(inputs.p_spec, 2.23);
        assert!(inputs.pl_only);
        assert!(inputs.do_progress_log);
        assert!(!inputs.debug_mode);
    }

    #[test]
    fn parses_basic_keyword_values() {
        let text = "\
NTIME 25
K_CBM 2
N_CBM 3.5
EPSB0 0.02
EPSE0 0.15
DOSSA 66
DOSSC 0
DEBUG 66
ENDIN
";

        let inputs = data_input_from_str(text).expect("data_input should succeed");

        assert_eq!(inputs.n_tm, 25);
        assert_eq!(inputs.k_cbm, 2);
        assert_eq!(inputs.n_cbm, 3.5);
        assert_eq!(inputs.eps_b0, 0.02);
        assert_eq!(inputs.epse0, 0.15);
        assert!(inputs.do_ssa);
        assert!(!inputs.do_ssc);
        assert!(inputs.debug_mode);
    }

    #[test]
    fn parses_multi_value_jet_keywords() {
        let text = "\
JTPRM 1.8 1.5 2.4
JTTHT 3.0 20.0
ENDIN
";

        let inputs = data_input_from_str(text).expect("data_input should succeed");

        assert_eq!(inputs.jet_delta, 1.8);
        assert_eq!(inputs.jet_f, 1.5);
        assert_eq!(inputs.jet_p, 2.4);
        assert_eq!(inputs.theta_j, 3.0);
        assert_eq!(inputs.theta_c, 20.0);
        assert!(inputs.do_jet_structure);
    }

    #[test]
    fn rejects_unknown_keyword() {
        let result = data_input_from_str("BADKY 1\nENDIN");

        assert!(matches!(
            result,
            Err(DataInputError::UnknownKeyword(keyword)) if keyword == "BADKY"
        ));
    }

    #[test]
    fn rejects_missing_value() {
        let result = data_input_from_str("NTIME\nENDIN");

        assert!(matches!(
            result,
            Err(DataInputError::MissingValue(keyword)) if keyword == "NTIME"
        ));
    }

    #[test]
    fn rejects_negative_theta_obs() {
        let result = data_input_from_str("THOBS -1.0\nENDIN");

        assert!(matches!(
            result,
            Err(DataInputError::InvalidValue { keyword, .. }) if keyword == "THOBS"
        ));
    }

    #[test]
    fn rejects_conflicting_redshift_and_luminosity_distance() {
        let text = "\
RDSHF 1.0
D_LUM 1.0e28
ENDIN
";

        let result = data_input_from_str(text);

        assert!(matches!(
            result,
            Err(DataInputError::ConflictingInputs { first, second })
                if first == "RDSHF" && second == "D_LUM"
        ));
    }

    #[test]
    fn rejects_conflicting_start_time_and_start_gamma() {
        let text = "\
GAMST 400.0
TSTRT 10.0
ENDIN
";

        let result = data_input_from_str(text);

        assert!(matches!(
            result,
            Err(DataInputError::ConflictingInputs { first, second })
                if first == "GAMST" && second == "TSTRT"
        ));
    }

    #[test]
    fn rejects_conflicting_end_time_and_end_gamma() {
        let text = "\
GMEND 2.5
T_END 1000.0
ENDIN
";

        let result = data_input_from_str(text);

        assert!(matches!(
            result,
            Err(DataInputError::ConflictingInputs { first, second })
                if first == "GMEND" && second == "T_END"
        ));
    }

    #[test]
    fn ntsig_less_than_or_equal_zero_requests_time_varying_index() {
        let inputs = data_input_from_str("NTSIG -1.0\nENDIN").expect("data_input should succeed");

        assert_eq!(inputs.p_spec, 0.0);
    }

    #[test]
    fn rejects_ntsig_between_zero_and_two() {
        let result = data_input_from_str("NTSIG 1.5\nENDIN");

        assert!(matches!(
            result,
            Err(DataInputError::InvalidValue { keyword, .. }) if keyword == "NTSIG"
        ));
    }
}
