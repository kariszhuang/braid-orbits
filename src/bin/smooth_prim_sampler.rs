// Coordination only: mathematical engines are shared with the serial census.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_small_census.rs");
    include!("../smooth_prim_sampler.rs");
}
fn main() {
    if let Err(e) = engine::prim_sampler_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("prim sampler: {e}");
        std::process::exit(1);
    }
}
