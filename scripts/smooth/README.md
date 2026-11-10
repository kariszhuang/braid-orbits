# All-smooth searches

Searches for all-smooth braid orbits: forest constructions (`run_forest_fast.py`), exact small-degree censuses, frontier and mechanism searches, and the independent audits that re-check every positive and every literal negative braid word.

| file | what it does |
|---|---|
| `analyze_smaller_actions.py` | Read-only exact group/position inspection of C141's two additional examples. |
| `audit_guided_native.py` | Independent streaming rejection audit; exceptional components use Python. |
| `audit_primitive_family.py` | Independent finite group/frame coverage and literal omitted-seed verification. |
| `benchmark_degree10_audit.py` | Matched prefix audit baseline; full-frame validation is separate. |
| `benchmark_forest_fast.py` | Bounded paired benchmark and independent controls; never launches discovery. |
| `benchmark_forest_hotpaths.py` | Paired fixed-proposal benchmarks: unchanged RNG, certificates and outcomes. |
| `benchmark_forest_middle.py` | Benchmark the changed Metal middle-square sieve, including checked rejection. |
| `benchmark_forest_second.py` | Bounded same-input second-pass timings, including independent verification. |
| `benchmark_forest_third.py` | Matched fixed-proposal third-pass timings, with independent omitted-seed checks. |
| `build_mechanism_inputs.py` | Freeze saved positive banks, mutation controls, and a source-bound profile queue. |
| `build_structural_novelty_input.py` | Freeze resolved-frontier and profile-cursor inputs for the authorized mission. |
| `checkpoint_smooth_handoff.py` | Explicit validated handoff; no deletion or modification of frozen inputs. |
| `compact_smooth_outputs.py` | Losslessly pack repeated smooth certificates and finished profile outputs. |
| `curate_discovery_design.py` | Independent group/block profile and deduplicated catalog for completed pilots. |
| `curate_ramified_double.py` | Full quotient factor maps and distinctness for the new ramified lift pilots. |
| `finalize_mechanism_frontier.py` | Bind ramified-lift discovery/proofs and final six-slot launch inputs. |
| `finalize_mechanism_launch.py` | Check the final launcher-only adjustments and bind preserved pilot leads. |
| `finalize_theory_frontier.py` | Bind final publication guard and mixed-degree partial audits to launch inputs. |
| `investigate_smooth_quotients.py` | Audit quotient genera and test two small assumption-breaking lift strata. |
| `prepare_design_leads.py` | Bind saved, independently checked partial pilot starts for later replay. |
| `profile_mechanism_breaker.py` | Incremental structural/marked-quotient audit of immutable smooth certificates. |
| `profile_smooth_actions.py` | Exact graph-derived passport and pure-braid fibers of audited components. |
| `profile_smooth_small_components.py` | Independent exact group orders and passport descriptions of smooth closures. |
| `profile_smooth_theory.py` | Independent block/quotient invariants for certified smooth components. |
| `record_quotient_framework.py` | Curate source-bound quotient/pilot/cleanup evidence without starting a search. |
| `record_structural_novelty.py` | Small source-bound handoff records; no discovery or external publication. |
| `recover_forest_interruption.py` | Finalize a dead smooth controller, retaining all evidence and pending leads. |
| `refine_smooth_patterns.py` | Exploratory elliptic-quotient pure-braid action tests; no theorem assertion. |
| `replay_saved_forest.py` | Bounded, component-deduplicated replay of VERIFIED forest timeouts. |
| `run_forest_fast.py` | Frozen optimized guided successor. |
| `run_negative_audit.py` | Resumable independent negative-only census audit. |
| `run_smooth_degree9.py` | Frozen degree-nine all-passport census; six Rust threads, then six audit processes. |
| `run_smooth_frontier.py` | 4+2 smooth frontier: exact positive-genus degree10, guided degree11--14. |
| `run_smooth_mission.py` | Frozen, shared-lock exact smooth-orbit mission. |
| `run_smooth_noncyclic.py` | Exact A5 block-quotient pilot; plan-only unless --run is supplied. |
| `run_smooth_small_census.py` | Shared-lock, frozen exhaustive smooth population in degrees 5--8. |
| `run_smooth_smaller.py` | Exact degree12 arbitrary-S4-voltage family; plan-only unless --run. |
| `run_smooth_special.py` | Frozen, shared-lock proof-directed special smoothness mission. |
| `smooth_block_lattice.py` | Complete invariant-partition lattice by joins of rooted orbit closures. |
| `validate_discovery_design.py` | Bounded discovery-design checks; one CPU and a held campaign lock. |
| `validate_final_design.py` | Final input-guard/control rebinding; no repeat of the discovery population. |
| `validate_mechanism_frontier.py` | Bound new constructors/runner/keys to finite regressions and preserved proofs. |
| `validate_native_smooth_components.py` | Validate independent native closures against bound Python-checked catalogs. |
| `validate_negative_audit.py` | Full independent small audit, full degree10 frames and measured paired checks. |
| `validate_structural_novelty.py` | Bounded constructor, omission, cache and lifecycle validation before launch. |
| `validate_theory_frontier.py` | Bounded theory/discovery checks and paired verified-bank workload timings. |
| `validate_triangle_target.py` | Focused prelaunch proof checks for the wide, three-shard TRIANGLE search. |
| `verified_smooth_bank.py` | Immutable per-degree positive snapshots; publish only fully audited closures. |
| `verify_smooth_geometry.py` | Independent exact hypotheses for the fixed-quotient geometric proof. |
| `verify_smooth_outputs.py` | Independent literal closure and SymPy group orders for a frozen smooth run. |
| `verify_smooth_small_run.py` | Freeze an independent re-audit of completed small-census outputs; no rerun. |

Tests (23): `test_compact_smooth_outputs.py`, `test_discovery_design.py`, `test_forest_fast.py`, `test_forest_launchd.py`, `test_forest_recovery.py`, `test_forest_runtime.py`, `test_guided_native.py`, `test_investigate_smooth_quotients.py`, `test_mechanism_breaker.py`, `test_native_components.py`, `test_negative_audit.py`, `test_saved_forest_replay.py`, `test_smooth_block_lattice.py`, `test_smooth_degree9.py`, `test_smooth_frontier.py`, `test_smooth_geometry.py`, `test_smooth_mission.py`, `test_smooth_noncyclic.py`, `test_smooth_small_census.py`, `test_smooth_smaller.py`, `test_smooth_special.py`, `test_structural_novelty.py`, `test_theory_frontier.py`
