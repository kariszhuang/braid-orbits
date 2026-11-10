// Separate catalog objective; historical discovery engines and binaries stay frozen.
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
    include!("../smooth_catalog_growth.rs");
}
fn main() {
    if let Err(e) = engine::catalog_growth_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("CATALOG_GROWTH_ATTENTION: {e}"); std::process::exit(1);
    }
}
