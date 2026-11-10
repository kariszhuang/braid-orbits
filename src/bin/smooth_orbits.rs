// Separate target: the frozen production weak-search sources stay unchanged.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
}
fn main() {
    if let Err(e) = engine::smooth_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("smooth: {e}");
        std::process::exit(1);
    }
}
