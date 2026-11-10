// A separate target keeps the frozen production search binary unchanged.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../linked_roots.rs");
}

fn main() {
    if let Err(error) = engine::linked_roots_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("linked roots: {error}");
        std::process::exit(1);
    }
}
