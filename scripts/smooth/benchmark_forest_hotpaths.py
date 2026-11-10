"""Paired fixed-proposal benchmarks: unchanged RNG, certificates and outcomes.

Runs serially so it does not create a second search pool. Existing campaign
load is retained and reported as a limitation. Both versions receive exactly
the same seed, limits and inputs; all persisted results are independently
audited. No discovery-rate inference is made from throughput.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
import random
import subprocess
import time

from audit_guided_native import audit_guided_native
from benchmark_forest_second import interval
from run_smooth_degree9 import sha

CASES = {
    "triangle14": (14, "triangle", 1),
    "triangle18": (18, "triangle", 1),
    "returns18": (18, "returns", 1),
    "returns24": (24, "returns", 1),
    "forest13": (13, "forest", 256),
    "forest18": (18, "forest", 256),
    "quotient24": (24, "quotient", 1),
}
STREAMS = ("FRAMES.tsv", "DESIGN.tsv", "POLICY.tsv", "MECHANISM.tsv",
           "LEADS.tsv", "RESOLUTIONS.tsv", "KNOWN.tsv")


def same_outcomes(old, new):
    a = {k: v for k, v in old["construction"].items() if k != "seconds"}
    b = {k: v for k, v in new["construction"].items() if k != "seconds"}
    if a != b or old["components"] != new["components"]:
        raise ValueError("candidate counters or independently audited outcomes changed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("old", "new", "auditor", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--attempts", type=int, default=200000)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--cases", nargs="+", choices=CASES,
                        default=[name for name in CASES if name != "quotient24"])
    args = parser.parse_args()
    if args.attempts < 1 or args.repeats < 3:
        parser.error("positive attempts and at least three paired replicates required")
    root = args.out.resolve()
    root.mkdir()
    old, new, auditor = (p.resolve() for p in (args.old, args.new, args.auditor))
    order_rng = random.Random(20261021)
    stages = []
    for name in args.cases:
        n, lane, reuse = CASES[name]
        timings = []
        for repeat in range(args.repeats):
            seed = 2026102100 + repeat
            order = ["old", "new"]
            order_rng.shuffle(order)
            pair, folders, reports = {}, {}, {}
            for variant in order:
                folder = root / f"{name}_{repeat}_{variant}"
                folders[variant] = folder
                command = [str(old if variant == "old" else new), "--d", str(n),
                           "--lane", lane, "--frame-reuse", str(reuse),
                           "--seconds", "120", "--attempts", str(args.attempts),
                           "--seed", str(seed), "--deep-level", "3",
                           "--compact-frames", "--out", str(folder)]
                began = time.monotonic()
                subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
                producer = time.monotonic() - began
                began = time.monotonic()
                report = audit_guided_native(folder, auditor, native_components=True)
                audit = time.monotonic() - began
                (folder / "CHECKS.json").write_text(json.dumps(report, indent=2) + "\n")
                reports[variant] = report
                count = report["construction"]
                if int(count["constructed"]) != args.attempts:
                    raise ValueError("time cap hit; not an equal-workload comparison")
                if int(count["unresolved"]):
                    raise ValueError("unresolved outcome; compare this workload separately")
                pair[variant] = dict(producer=producer, audit=audit, total=producer + audit)
            same_outcomes(reports["old"], reports["new"])
            for stream in STREAMS:
                if sha(folders["old"] / stream) != sha(folders["new"] / stream):
                    raise ValueError("RNG/proposal/proof stream changed: " + stream)
            timings.append(dict(seed=seed, order=order, **pair))
            print(name, repeat, "producer", pair["old"]["producer"] / pair["new"]["producer"], flush=True)
        stages.append(dict(case=name, degree=n, lane=lane, reuse=reuse, timings=timings,
                           producer=interval([v["old"]["producer"] / v["new"]["producer"] for v in timings]),
                           pipeline=interval([v["old"]["total"] / v["new"]["total"] for v in timings])))
        # Persist each completed stage; a later failed test never erases evidence.
        (root / "STAGES.json").write_text(json.dumps(stages, indent=2) + "\n")
    result = dict(passed=True, stages=stages, attempts=args.attempts, repeats=args.repeats,
                  binary_sha256=sha(new), old_binary_sha256=sha(old), audit_binary_sha256=sha(auditor),
                  identical_rng_proposals_counters_certificates=True,
                  every_persisted_outcome_independently_verified=True,
                  limitations="Busy shared Mac; fixed biased streams, no known-bank startup included. "
                              "Producer includes process startup and T7 logging; pipeline includes independent audit. "
                              "Bootstrap paired-median intervals are descriptive, not guaranteed future throughput or novelty yield.")
    (root / "RESULTS.json").write_text(json.dumps(result, indent=2) + "\n")
    print("RESULTS=" + str(root / "RESULTS.json"), flush=True)


if __name__ == "__main__":
    main()
