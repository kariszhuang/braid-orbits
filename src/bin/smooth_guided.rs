#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_small_census.rs");
    include!("../smooth_guided.rs");
}
fn main() {
    if let Err(e) = engine::guided_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("guided smooth search: {e}");
        std::process::exit(1);
    }
}
