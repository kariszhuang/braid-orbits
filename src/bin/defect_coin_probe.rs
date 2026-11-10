// Read-only bounded statistic probe. Shared permutation/braid kernels; no hunt.
#[allow(dead_code)]
mod engine {
    include!("../main.rs");
    include!("../defect_coin_probe.rs");
}
fn main() {
    if let Err(e) = engine::coin_probe_main(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("defect coin probe: {e}");
        std::process::exit(1);
    }
}
