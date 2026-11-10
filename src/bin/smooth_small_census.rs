// Reuse the exact permutation/closure engine without changing the weak search.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_small_census.rs");
}
fn main() {
    if let Err(e) = engine::small_smooth_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("small smooth census: {e}");
        std::process::exit(1);
    }
}
