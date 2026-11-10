# Scripts

Python 3 (standard library, except where a script imports SymPy or NumPy). Run everything from the repository
root, for example `python3 scripts/smooth/run_forest_fast.py --help`: paths such as `data/...`, `src/...` and
`target/...` are relative to it.

| folder | contents | files |
|---|---|---|
| [`smooth/`](smooth/README.md) | All-smooth searches | 74 |
| [`genus3/`](genus3/README.md) | Genus-three all-smooth search | 76 |
| [`weak/`](weak/README.md) | Weak pair-generation failures | 27 |
| [`spectra/`](spectra/README.md) | Braid-graph spectra and defect statistics | 19 |
| [`catalog/`](catalog/README.md) | Public catalogue | 22 |
| [`gpu/`](gpu/README.md) | GPU prototype | 15 |

Naming: `run_*` launches a search and writes to `data/`; `verify_*`, `audit_*`, `check_*` and `validate_*` re-check
saved results independently of the code that produced them; `benchmark_*` measures speed; `record_*`, `curate_*` and
`finalize_*` summarize a finished run; `test_X.py` tests `X.py` and sits next to it.

## Imports between folders

Scripts import each other by module name (`from run_forest_fast import ...`) across folders. Each folder holds an
identical `_scripts_path.py`, and a script that needs another one starts with `import _scripts_path`, which puts
all six folders on the import path. `_scripts_path.find(name)` returns the path of a script in any folder.
Some launchers freeze a copy of the scripts into their run folder before starting. The import is written so that
such a copy still runs.
