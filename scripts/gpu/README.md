# GPU prototype

An Apple Metal screening prototype (`*.metal`, `*.swift`) with its CPU reference (`*_cpu.rs`) and the CPU/GPU hybrid runner. It was slower than the CPU path, so the searches do not use it.

| file | what it does |
|---|---|
| `benchmark_finite_prefix.py` | Bounded shared-buffer CPU4/Metal screening comparison from a real finite case. |
| `benchmark_gpu.py` | Finite exact GPU-screen benchmark. |
| `check_campaign_hybrid.py` | Bounded paired CPU/Metal controls using an already frozen campaign runtime. |
| `finite_prefix.metal` | — |
| `finite_prefix.swift` | — |
| `finite_prefix_cpu.rs` | — |
| `forest_middle_screen.metal` | — |
| `forest_middle_screen.swift` | — |
| `forest_screen.metal` | — |
| `forest_screen.swift` | — |
| `gpu_screen.metal` | — |
| `gpu_screen.swift` | — |
| `gpu_screen_cpu.rs` | — |
| `hybrid_pilot.py` | Frozen native hybrid runtime and paired experiments for run_campaign.py. |
| `hybrid_pilot.swift` | — |


