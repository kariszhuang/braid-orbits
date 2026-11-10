// A separate mathematical target; production and prior frozen missions unchanged.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../smooth_orbits.rs");
    include!("../smooth_special.rs");
}
fn main() {
    if let Err(error) = engine::special_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("smooth-special: {error}");
        std::process::exit(1);
    }
}
