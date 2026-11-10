#!/usr/bin/env python3
"""Compute all exact S8 shard spectra after the exact census finishes."""

import csv
import shutil
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXACT = ROOT / "runs/exact/S8_shards"
SPECTRA = ROOT / "runs/spectra/S8_exact_shards"
PILOT = ROOT / "runs/exact/s8-pilot-20260927-shard0"
BIN = ROOT / "target/release/hurwitz_breaker"


def main() -> None:
    deadline = time.monotonic() + 24 * 3600
    while not all((EXACT / f"shard_{i}/SUMMARY.txt").is_file() for i in range(1, 16)):
        if time.monotonic() > deadline:
            raise RuntimeError("Exact S8 census did not finish within 24 hours; spectra not started")
        time.sleep(30)

    shard0 = EXACT / "shard_0"
    shard0.mkdir(parents=True, exist_ok=True)
    for name in ("SUMMARY.txt", "generating_components.csv"):
        shutil.copy2(PILOT / name, shard0 / name)

    for i in range(16):
        source = EXACT / f"shard_{i}/generating_components.csv"
        out = SPECTRA / f"shard_{i}"
        if (out / "SUMMARY.txt").is_file():
            print(f"S8 spectral shard {i}/16 already complete", flush=True)
            continue
        out.mkdir(parents=True, exist_ok=True)
        command = [str(BIN), "spectrum", "--components", str(source),
                   "--out", str(out), "--threads", "8", "--full-max", "256",
                   "--lanczos-steps", "192", "--baseline-reps", "0",
                   "--baseline-kinds", "", "--summary-only"]
        print(f"starting S8 spectral shard {i}/16", flush=True)
        with (out / "RUN.log").open("a", encoding="utf-8") as log:
            subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
        with source.open(newline="", encoding="utf-8") as stream:
            expected = sum(1 for _ in csv.DictReader(stream))
        with (out / "SPECTRUM.tsv").open(newline="", encoding="utf-8") as stream:
            actual = sum(1 for _ in csv.DictReader(stream, delimiter="\t"))
        if actual != expected:
            raise RuntimeError(f"S8 spectral shard {i}: expected {expected} cases, got {actual}")
        print(f"completed S8 spectral shard {i}/16: {actual} cases", flush=True)

    subprocess.run(["uv", "run", "--isolated", "--no-project", "--python", "3.13",
                    "--with", "matplotlib==3.11.1", "--with", "numpy", "python",
                    "scripts/spectra/make_all_small_spectra.py"], cwd=ROOT, check=True)
    print("Complete S4–S8 figure generated", flush=True)


if __name__ == "__main__":
    main()
