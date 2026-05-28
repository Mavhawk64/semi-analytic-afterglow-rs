fn main() {
    if let Err(error) = semi_analytic_afterglow_rs::saafterglow2::run() {
        eprintln!("{error:?}");
        std::process::exit(1);
    }
}
