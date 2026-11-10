#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_small_census.rs");
    include!("../smooth_guided.rs");
    include!("../smooth_special.rs");
    include!("../smooth_noncyclic.rs");
    include!("../smooth_discovery_design.rs");
    include!("../smooth_quotient_forest.rs");
    include!("../smooth_mechanism_breaker.rs");
    include!("../smooth_ramified_double.rs");
    include!("../smooth_structural_novelty.rs");
    include!("../smooth_forest_fast.rs");
}
fn main() {
    if let Err(e) = engine::forest_fast_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("fast smooth forest: {e}");
        std::process::exit(1);
    }
}
