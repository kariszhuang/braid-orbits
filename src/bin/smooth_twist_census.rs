// Coordination only: mathematical engines are shared with the serial census.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_small_census.rs");
    include!("../smooth_star_census.rs");
    include!("../smooth_twist_census.rs");
}
fn main() {
    if let Err(e) = engine::twist_census_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("twist census: {e}");
        std::process::exit(1);
    }
}
