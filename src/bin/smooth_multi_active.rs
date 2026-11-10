// Multi-active-point unramified lift census (see src/smooth_multi_active.rs).
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_small_census.rs");
    include!("../smooth_star_census.rs");
    include!("../smooth_multi_active.rs");
}
fn main() {
    if let Err(e) = engine::multi_active_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("multi active: {e}");
        std::process::exit(1);
    }
}
