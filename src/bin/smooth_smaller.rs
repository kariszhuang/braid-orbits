// Separate target: existing production campaigns and pinned cursors are untouched.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_special.rs");
    include!("../smooth_noncyclic.rs");
    include!("../smooth_smaller.rs");
}
fn main() {
    if let Err(error)=engine::smaller_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("smooth-smaller: {error}");
        std::process::exit(1);
    }
}
