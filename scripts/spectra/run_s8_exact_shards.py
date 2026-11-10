#!/usr/bin/env python3
"""Finish the 16 exact S8 signature shards, resuming completed checkpoints."""

import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from queue import Empty, Queue


ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "target/release/hurwitz_breaker"
BASE = ROOT / "runs/exact/S8_shards"


def run_shard(index: int) -> None:
    out = BASE / f"shard_{index}"
    if (out / "SUMMARY.txt").is_file():
        print(f"shard {index}/16 already complete", flush=True)
        return
    out.mkdir(parents=True, exist_ok=True)
    command = [str(BIN), "exact", "--group", "S", "--d", "8",
               "--threads", "8", "--out", str(out), "--checkpoint-every", "500",
               "--keep-going", "--shard-count", "16", "--shard-index", str(index)]
    if (out / "checkpoint.bin").is_file():
        command.append("--resume")
    print(f"starting exact S8 shard {index}/16", flush=True)
    with (out / "RUN.log").open("a", encoding="utf-8") as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not (out / "SUMMARY.txt").is_file():
        raise RuntimeError(f"S8 shard {index} returned without a completion summary")
    print(f"completed exact S8 shard {index}/16", flush=True)


def main() -> None:
    if not BIN.is_file():
        raise SystemExit(f"Build the release binary first: {BIN}")
    pending: Queue[int] = Queue()
    for index in range(5, 16):
        pending.put(index)

    def worker() -> None:
        while not all((BASE / f"shard_{i}/SUMMARY.txt").is_file() for i in (2, 3, 4)):
            time.sleep(20)
        while True:
            try:
                index = pending.get_nowait()
            except Empty:
                return
            run_shard(index)

    # Shards 2-4 are already running separately. Let them finish, then keep
    # at most two processes active to avoid heavy memory compression and swap.
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(worker) for _ in range(2)]
        for future in futures:
            future.result()


if __name__ == "__main__":
    try:
        main()
    except (OSError, subprocess.CalledProcessError, RuntimeError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
