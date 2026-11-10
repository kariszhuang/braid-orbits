//! Exact finite-target experiment; independent of the live producer/dirty main.rs.
#[path = "../genus_three_finite.rs"]
mod finite;
fn main() {
    if let Err(error) = finite::main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("genus-three-finite: {error}");
        std::process::exit(1);
    }
}
