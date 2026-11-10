// Isolated target: never changes a frozen campaign's executable or cursors.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_special.rs");
    include!("../smooth_noncyclic.rs");
}
fn main() {
    if let Err(error) = engine::noncyclic_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("smooth-noncyclic: {error}");
        std::process::exit(1);
    }
}
