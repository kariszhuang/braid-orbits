#!/usr/bin/env python3
"""Build, run, and analyze the finite Rust spectral census."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass

import argparse
from datetime import datetime
import os
from pathlib import Path
import subprocess
import sys


PROJECT = Path(__file__).resolve().parents[2]
COMPONENTS = PROJECT / "runs" / "campaigns" / "campaign_20260907_102423" / "exact"
TARGET = PROJECT / "target" / "spectrum"
sys.path.insert(0, str(PROJECT))
from run_campaign import require_results_volume, RESULTS


def component_files(groups: str, low: int, high: int) -> tuple[list[Path], list[str]]:
    found, missing = [], []
    for group in groups.split(","):
        for degree in range(low, high + 1):
            direct = COMPONENTS / f"{group}{degree}" / "generating_components.csv"
            shards = sorted((COMPONENTS / f"{group}{degree}_shards").glob("shard_*/generating_components.csv"))
            if direct.exists():
                found.append(direct)
            elif shards:
                found.extend(shards)
            elif group == 'S' and degree == 8:
                # The completed S8 census is newer than the legacy catalog.
                base = PROJECT / 'runs/exact/S8_shards'
                catalogs = [base / f'shard_{i}/generating_components.csv' for i in range(16)]
                summaries = [base / f'shard_{i}/SUMMARY.txt' for i in range(16)]
                if all(p.is_file() for p in catalogs + summaries):
                    for i, summary in enumerate(summaries):
                        values = dict(l.split(' = ', 1) for l in summary.read_text().splitlines() if ' = ' in l)
                        if values.get('group') != 'S_8' or values.get('signature shard') != f'{i}/16':
                            raise ValueError(f'Wrong completed S8 shard: {summary}')
                    found.extend(catalogs)
                else:
                    missing.append('S8')
            else:
                missing.append(f"{group}{degree}")
    return found, missing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", default="A,S", help="comma-separated ambient groups")
    parser.add_argument("--min-degree", type=int, default=5)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--full-max", type=int, default=256)
    parser.add_argument("--degree5-full-max", type=int, default=1024)
    parser.add_argument("--lanczos-steps", type=int, default=192)
    parser.add_argument("--baseline-reps", type=int, default=2)
    parser.add_argument("--raw-graphs", action="store_true")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    files, missing = component_files(args.groups, args.min_degree, args.max_degree)
    if not files:
        raise SystemExit("no exact component catalogs found")
    out = args.output or RESULTS / "spectra" / (f"degrees-{args.min_degree}-{args.max_degree}-" + datetime.now().strftime("%Y%m%d-%H%M%S"))
    require_results_volume(out)
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CARGO_TARGET_DIR=str(TARGET), RUSTFLAGS="-C target-cpu=native")
    if not args.skip_tests:
        subprocess.run(["cargo", "test", "--release"], cwd=PROJECT, env=env, check=True)
    subprocess.run(["cargo", "build", "--release"], cwd=PROJECT, env=env, check=True)
    binary = TARGET / "release" / "hurwitz_breaker"
    cmd = [
        str(binary), "spectrum",
        "--components", ",".join(map(str, files)),
        "--out", str(out),
        "--threads", str(args.threads),
        "--full-max", str(args.full_max),
        "--degree5-full-max", str(args.degree5_full_max),
        "--lanczos-steps", str(args.lanczos_steps),
        "--baseline-reps", str(args.baseline_reps),
    ]
    if not args.raw_graphs:
        cmd.append("--summary-only")
    (out / "INPUTS.txt").write_text("\n".join(map(str, files)) + "\n", encoding="utf-8")
    (out / "MISSING.txt").write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8")
    (out / "COMMAND.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
    print(f"output={out}", flush=True)
    print(f"catalogs={len(files)} missing={','.join(missing) or 'none'}", flush=True)
    with (out / "RUN.log").open("w") as log:
        subprocess.run(cmd, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run([
        "uv", "run", "--isolated", "--no-project", "--python", "3.13",
        "--with", "numpy==2.3.5", "--with", "matplotlib==3.11.1", "--with", "scipy==1.16.3",
        "python", str(_scripts_path.find("analyze_spectra.py")), str(out),
    ], cwd=PROJECT, check=True)
    print(out)


if __name__ == "__main__":
    main()
