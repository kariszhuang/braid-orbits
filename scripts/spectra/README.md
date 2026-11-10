# Braid-graph spectra and defect statistics

Spectra of braid graphs (exact census S8, gap band certificates) and defect distribution studies.

| file | what it does |
|---|---|
| `analyze_gap_band_bipartite.py` | Audit the 0.19 gap band and exact versus near bipartiteness in the saved census. |
| `analyze_spectra.py` | Reproducible matched-null and genus analysis for a Rust spectral census. |
| `analyze_spectral_followups.py` | Analyze structural spectral follow-ups with matched controls and diagnostics. |
| `analyze_weekly_graph_bridge.py` | Recheck the A25 graph bridge and the complete A5/S5 spectral census. |
| `certify_gap_band.py` | Exact characteristic-polynomial proof for the 55 larger and 13 degree-four 0.195262 cases. |
| `compute_spectra.py` | Build, run, and analyze the finite Rust spectral census. |
| `make_all_small_spectra.py` | Figures of the braid-graph spectra for S4 to S8, from saved spectral census output (writes to `data/figures`). |
| `defect_coin_followup.py` | Bounded one-CPU coin-model follow-up; plan-only unless --run is given. |
| `defect_degree_stability.py` | Bounded, ONE auxiliary CPU; plan-only unless --run. |
| `plot_defect_coin_followup.py` | Provisional research figure; exact data and correlated walk estimates separated. |
| `plot_defect_distributions.py` | Plot exact defect histograms from defect-census TSV (no sampling). |
| `run_exact_spectral_certificates.py` | Rerun small recurring spectra and certify their labeled graph classes exactly. |
| `run_s8_exact_shards.py` | Finish the 16 exact S8 signature shards, resuming completed checkpoints. |
| `run_s8_spectra_after_exact.py` | Compute all exact S8 shard spectra after the exact census finishes. |
| `run_spectral_followups.py` | Run the staged structural spectral campaign without overspending every test on every graph. |

Tests (4): `test_defect_coin_followup.py`, `test_defect_coin_model.py`, `test_defect_degree_stability.py`, `test_spectral_theorem_logic.py`
