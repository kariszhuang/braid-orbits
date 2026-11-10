#!/usr/bin/env python3
"""Run the staged structural spectral campaign without overspending every test on every graph."""

from __future__ import annotations
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass

import argparse
import csv
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import sys


PROJECT = Path(__file__).resolve().parents[2]
EXACT = PROJECT / "runs" / "campaigns" / "campaign_20260907_102423" / "exact"
SMART = PROJECT / "runs" / "campaigns" / "smart_20260907_200313_392975" / "jobs"
TARGET = PROJECT / "target" / "spectrum"


def component_files() -> list[Path]:
    files = []
    for group, maximum in (("A", 8), ("S", 7)):
        for degree in range(5, maximum + 1):
            direct = EXACT / f"{group}{degree}" / "generating_components.csv"
            shards = sorted((EXACT / f"{group}{degree}_shards").glob("shard_*/generating_components.csv"))
            if direct.exists():
                files.append(direct)
            else:
                files.extend(shards)
    return files


def write_s8_sample(root: Path, count: int) -> Path:
    rows, seen = [], set()
    for path in sorted(SMART.glob("*S8*/attempt_*/RESULTS.tsv")):
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                seed = row.get("seed", "").strip()
                if not seed or seed in seen or row.get("status") not in {"strong_found", "weak_found", "certified_strong_failure", "certified_weak_and_strong_failure"}:
                    continue
                seen.add(seed)
                rows.append(seed)
                if len(rows) >= count:
                    break
        if len(rows) >= count:
            break
    if not rows:
        raise RuntimeError("no completed S8 result seeds found")
    path = root / "inputs" / "S8" / "generating_components.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["component", "seed"])
        writer.writeheader()
        for index, seed in enumerate(rows, 1):
            writer.writerow({"component": index, "seed": seed})
    return path


def load_cases(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stratified_ids(rows: list[dict[str, str]], limit: int, max_vertices: int | None = None, iterative_only: bool = False) -> list[str]:
    eligible = []
    for row in rows:
        if row.get("kind") != "exact_component":
            continue
        vertices = int(row["quotient_states"])
        if max_vertices is not None and vertices > max_vertices:
            continue
        if iterative_only and not row["method"].startswith("lanczos"):
            continue
        eligible.append(row)
    strata: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in eligible:
        strata.setdefault((row["group"], row["degree"], row["genus"]), []).append(row)
    queues = []
    for key in sorted(strata):
        block = sorted(strata[key], key=lambda row: int(row["quotient_states"]))
        positions = sorted(set([0, len(block) // 4, len(block) // 2, 3 * len(block) // 4, len(block) - 1]))
        queues.append([block[position]["case_id"] for position in positions])
    selected = []
    while len(selected) < limit and any(queues):
        next_queues = []
        for queue in queues:
            if queue and len(selected) < limit:
                case_id = queue.pop(0)
                if case_id not in selected:
                    selected.append(case_id)
            if queue:
                next_queues.append(queue)
        queues = next_queues
    return selected


def write_ids(path: Path, ids: list[str]) -> None:
    path.write_text("\n".join(ids) + "\n", encoding="utf-8")


def status(root: Path, stage: str, state: str, **extra: object) -> None:
    payload = {"stage": stage, "state": state, "updated": datetime.now().astimezone().isoformat(), **extra}
    (root / "STATUS.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_logged(root: Path, stage: str, command: list[str], env: dict[str, str]) -> None:
    destination = root / stage
    destination.mkdir(parents=True, exist_ok=False)
    (destination / "COMMAND.txt").write_text(" ".join(command) + "\n", encoding="utf-8")
    status(root, stage, "running", output=str(destination), command=command)
    with (destination / "RUN.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=PROJECT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    status(root, stage, "complete", output=str(destination))


def spectrum_command(binary: Path, files: list[Path], output: Path, threads: int, *, selection: Path | None = None, baseline_reps: int = 2, baseline_kinds: str = "random_6_regular,random_schreier,cycle_matched_schreier", full_max: int = 256, trace_order: int = 0, trace_probes: int = 8, maxcut_restarts: int = 0, switch_attempts: int = 20, lanczos_checks: str = "") -> list[str]:
    command = [
        str(binary), "spectrum", "--components", ",".join(map(str, files)), "--out", str(output),
        "--threads", str(threads), "--full-max", str(full_max), "--degree5-full-max", "1024",
        "--lanczos-steps", "192", "--baseline-reps", str(baseline_reps), "--baseline-kinds", baseline_kinds,
        "--trace-order", str(trace_order), "--trace-probes", str(trace_probes), "--maxcut-restarts", str(maxcut_restarts),
        "--switch-attempts-per-edge", str(switch_attempts), "--summary-only",
    ]
    if selection:
        command += ["--select-case-ids", str(selection)]
    if lanczos_checks:
        command += ["--lanczos-checks", lanczos_checks]
    return command


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-census", type=Path, default=PROJECT / "runs" / "spectra" / "degrees-5-8-full-20260922")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--s8-samples", type=int, default=24)
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    root = (args.output or PROJECT / "runs" / "spectra" / f"structural-followups-{datetime.now():%Y%m%d-%H%M%S}").resolve()
    root.mkdir(parents=True, exist_ok=False)
    try:
        status(root, "setup", "running")
        env = dict(os.environ, CARGO_TARGET_DIR=str(TARGET), RUSTFLAGS="-C target-cpu=native")
        if not args.skip_tests:
            subprocess.run(["cargo", "test", "--release"], cwd=PROJECT, env=env, check=True)
        subprocess.run(["cargo", "build", "--release"], cwd=PROJECT, env=env, check=True)
        binary = TARGET / "release" / "hurwitz_breaker"

        s8 = write_s8_sample(root, args.s8_samples)
        files = component_files() + [s8]
        write_ids(root / "S8_CASES.txt", [f"S8:component-{index}" for index in range(1, args.s8_samples + 1)])
        base_rows = load_cases(args.base_census / "SPECTRUM.tsv")
        high_ids = stratified_ids(base_rows, 140, max_vertices=50_000)
        trace_ids = stratified_ids([row for row in base_rows if row["case_id"] in set(high_ids)], 56, max_vertices=5_000)
        switch_ids = trace_ids[:24]
        convergence_ids = stratified_ids(base_rows, 32, iterative_only=True)
        high_path = root / "HIGH_REP_CASES.txt"; write_ids(high_path, high_ids)
        trace_path = root / "TRACE_CASES.txt"; write_ids(trace_path, trace_ids)
        switch_path = root / "SWITCH_CASES.txt"; write_ids(switch_path, switch_ids)
        convergence_path = root / "CONVERGENCE_CASES.txt"; write_ids(convergence_path, convergence_ids)

        if args.plan_only:
            status(root, "setup", "planned", catalogs=len(files), s8_samples=args.s8_samples, high_replicate_cases=len(high_ids), trace_cases=len(trace_ids), switch_cases=len(switch_ids), convergence_cases=len(convergence_ids))
            print(root)
            return

        analysis_env = dict(env)
        analysis_prefix = ["uv", "run", "--isolated", "--no-project", "--python", "3.13", "--with", "numpy==2.3.5", "--with", "pandas==2.3.3", "--with", "scipy==1.16.3", "--with", "statsmodels==0.14.6", "--with", "matplotlib==3.11.1", "python"]
        preliminary = root / "preliminary-analysis"
        preliminary.mkdir()
        subprocess.run(analysis_prefix + [str(_scripts_path.find("analyze_spectral_followups.py")), "--census", str(args.base_census), "--output", str(preliminary)], cwd=PROJECT, env=analysis_env, check=True)

        run_logged(root, "full-census", spectrum_command(binary, files, root / "full-census", args.threads), env)
        run_logged(root, "high-replicates", spectrum_command(binary, files, root / "high-replicates", args.threads, selection=high_path, baseline_reps=32, full_max=512), env)
        run_logged(root, "trace-and-bipartite", spectrum_command(binary, files, root / "trace-and-bipartite", args.threads, selection=trace_path, baseline_reps=8, full_max=512, trace_order=8, trace_probes=8, maxcut_restarts=4, switch_attempts=100), env)
        run_logged(root, "lanczos-convergence", spectrum_command(binary, files, root / "lanczos-convergence", args.threads, selection=convergence_path, baseline_reps=0, baseline_kinds="", lanczos_checks="96,192,384,768"), env)
        for attempts in (20, 100, 500):
            stage = f"switch-{attempts}"
            run_logged(root, stage, spectrum_command(binary, files, root / stage, args.threads, selection=switch_path, baseline_reps=16, baseline_kinds="random_6_regular", full_max=512, switch_attempts=attempts), env)

        final = root / "analysis"
        final.mkdir()
        command = analysis_prefix + [
            str(_scripts_path.find("analyze_spectral_followups.py")),
            "--census", str(root / "full-census"), "--deep", str(root / "trace-and-bipartite"),
            "--convergence", str(root / "lanczos-convergence"), "--high-reps", str(root / "high-replicates"),
            "--switch-run", f"20={root / 'switch-20'}", "--switch-run", f"100={root / 'switch-100'}",
            "--switch-run", f"500={root / 'switch-500'}", "--output", str(final),
        ]
        subprocess.run(command, cwd=PROJECT, env=analysis_env, check=True)
        status(root, "analysis", "complete", output=str(final))
    except Exception as error:
        status(root, "failed", "failed", error=repr(error))
        raise
    print(root)


if __name__ == "__main__":
    main()
