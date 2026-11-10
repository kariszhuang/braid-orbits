# Evidence inputs

These are not results to read. They are the few saved files that the code or its tests open by name, copied unchanged
from the original project. Several are bound by SHA-256 hashes that the code checks, so do not edit them.

| folder | used by |
|---|---|
| `certificates/` | The nine certified weak pair-generation failures (A25 and S18). Rust tests compile them in; `scripts/weak` checks them. |
| `spectral_census/` | S18 position spectrum used by a Rust test; gap-band certificates for `scripts/spectra`. |
| `overnight_genus3_20261023/actions/`, `resumed_genus3_20261023/` | Code and hash records of the degree-12 and degree-13 genus-3 exclusions. `run_genus_three.py` checks them before admitting those results. |
| `restriction_admission_20261024/`, `genus_three_twists_20261022/` | Proof admission records read by `run_genus_three.py`. |
| all other dated folders | Earlier catalogues, plans and math notes that a launcher or verifier reads (as a starting point, or to check a hash), and fixtures for their tests. Usually one or two files each. |
