fn main() {
    if let Err(error) = semi_analytic_afterglow_rs::saa::run() {
        eprintln!("{error:?}");
        std::process::exit(1);
    }
}
