#!/usr/bin/env python3
"""Audit the 0.19 gap band and exact versus near bipartiteness in the saved census."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import Counter, deque
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FULL = ROOT / "runs/spectra/structural-followups-20260922-ecf2d7a/full-census"
CERT = ROOT / "evidence/spectral_census/EXACT_LAMBDA2_CERTIFICATES.tsv"
OUT = ROOT / "runs/spectra/gap-band-bipartite-20260927"
TARGET = (2 - np.sqrt(2)) / 3


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def p4_fibers(transitions: np.ndarray) -> list[int] | None:
    """Find a move-label-preserving map to the four exceptional positions."""
    n = len(transitions)
    for root_position in range(4):
        positions = [-1] * n
        positions[0] = root_position
        queue = deque([0])
        valid = True
        while queue and valid:
            vertex = queue.popleft()
            for move, neighbor in enumerate(transitions[vertex]):
                i = move % 3
                expected = (positions[vertex] + 1 if positions[vertex] == i else
                            positions[vertex] - 1 if positions[vertex] == i + 1 else positions[vertex])
                if positions[neighbor] == -1:
                    positions[neighbor] = expected
                    queue.append(int(neighbor))
                elif positions[neighbor] != expected:
                    valid = False
                    break
        if valid and set(positions) == set(range(4)):
            return positions
    return None


def certified_quotients() -> list[dict[str, object]]:
    output = []
    q5_expected = ((4, 2, 0, 0, 0), (2, 0, 4, 0, 0),
                   (0, 2, 2, 2, 0), (0, 0, 4, 0, 2), (0, 0, 0, 2, 4))
    for row in rows(CERT):
        n = int(row["vertices"])
        transitions = np.array([int(x) for x in row["canonical_labeled_code"].split(",")], dtype=int).reshape(n, 6)
        fibers = p4_fibers(transitions)
        if fibers is not None:
            quotient = np.zeros((4, 4), dtype=int)
            for vertex, position in enumerate(fibers):
                counts = np.bincount([fibers[j] for j in transitions[vertex]], minlength=4)
                if not quotient[position].any():
                    quotient[position] = counts
                else:
                    assert np.array_equal(quotient[position], counts)
            assert np.array_equal(quotient, [[4, 2, 0, 0], [2, 2, 2, 0],
                                             [0, 2, 2, 2], [0, 0, 2, 4]])
            kind = "labeled_P4"
        else:
            adjacency = np.zeros((n, n), dtype=int)
            for vertex, targets in enumerate(transitions):
                for neighbor in targets:
                    adjacency[vertex, neighbor] += 1
            eigenvalues, eigenvectors = np.linalg.eigh(adjacency)
            vector = eigenvectors[:, np.argmin(abs(eigenvalues - (2 + 2 * np.sqrt(2))))]
            assert abs(adjacency @ vector - (2 + 2 * np.sqrt(2)) * vector).max() < 1e-10
            levels = np.round(vector / max(abs(vector)), 6)
            assert len(set(levels)) == 5
            classes = {level: np.flatnonzero(levels == level).tolist() for level in sorted(set(levels))}
            quotient = np.zeros((5, 5), dtype=int)
            for i, vertices in enumerate(classes.values()):
                for vertex in vertices:
                    counts = [sum(int(target in classes[level]) for target in transitions[vertex])
                              for level in classes]
                    if vertex == vertices[0]:
                        quotient[i] = counts
                    else:
                        assert np.array_equal(quotient[i], counts)
            assert tuple(map(tuple, quotient)) == q5_expected
            kind = "unlabeled_equitable_5"
        assert np.min(abs(np.linalg.eigvals(quotient) - (2 + 2 * np.sqrt(2)))) < 1e-10
        output.append({"representative": row["representative_case_id"], "vertices": n,
                       "cases": int(row["case_count"]), "quotient": kind})
    assert len(output) == 8
    assert Counter(item["quotient"] for item in output) == {"labeled_P4": 6, "unlabeled_equitable_5": 2}
    return output


def main() -> None:
    spectra = rows(FULL / "SPECTRUM.tsv")
    diagnostics = {row["case_id"]: row for row in rows(FULL / "GRAPH_DIAGNOSTICS.tsv") if row["source"] == "case"}
    cases = []
    for row in spectra:
        if row["group"] == "S" and row["degree"] == "8":
            continue  # S8 is a sample, not a complete degree census.
        gap = float(row["normalized_gap"])
        absolute = 1 - float(row["abs_nontrivial_radius"]) / 6
        bipartite = diagnostics[row["case_id"]]["exact_bipartite"] == "true"
        assert (absolute < 1e-8) == bipartite
        cases.append({"id": row["case_id"], "group": row["group"] + row["degree"],
                      "genus": int(row["genus"]), "vertices": int(row["quotient_states"]),
                      "passport": row["passport"], "method": row["method"],
                      "gap": gap, "absolute": max(0., absolute), "bipartite": bipartite})
    assert len(cases) == 7319
    band = [case for case in cases if .18 <= case["gap"] < .20]
    target = [case for case in band if abs(case["gap"] - TARGET) < 1e-8]
    exceptions = [case for case in band if case not in target]
    near = [case for case in cases if case["absolute"] < .01]
    assert (len(band), len(target), len(exceptions), len(near)) == (2817, 2812, 5, 3731)
    assert sum(case["bipartite"] for case in near) == 3717
    quotients = certified_quotients()
    assert sum(item["cases"] for item in quotients) == 2757

    OUT.mkdir(parents=True, exist_ok=True)
    for name, selected in (("BAND_CASES.tsv", band),
                           ("NEAR_NONBIPARTITE.tsv", [c for c in near if not c["bipartite"]])):
        with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(cases[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(selected)
    large_ids = {case["id"] for case in target if case["vertices"] > 24}
    (OUT / "TARGET_55.ids").write_text("\n".join(sorted(large_ids)) + "\n", encoding="utf-8")
    if "--rerun-55" in sys.argv:
        input_sources = sorted({row["source"] for row in spectra if row["case_id"] in large_ids})
        assert len(input_sources) == 12 and len(large_ids) == 55
        subprocess.run([str(ROOT / "target/release/hurwitz_breaker"), "spectrum",
                        "--components", ",".join(input_sources), "--select-case-ids", str(OUT / "TARGET_55.ids"),
                        "--out", str(OUT / "target55-rerun"), "--threads", "4", "--full-max", "128",
                        "--degree5-full-max", "128", "--lanczos-steps", "192",
                        "--baseline-reps", "0", "--canonical-max", "96", "--summary-only"],
                       cwd=ROOT, check=True)
    report = {"scope": "complete stored A5-A8/S5-S7; sampled S8 excluded",
              "cases": len(cases), "band_0.18_to_0.20": len(band),
              "numerical_target_within_1e-8": len(target), "other_band_values": exceptions,
              "target_graphs_over_24_vertices": sum(case["vertices"] > 24 for case in target),
              "exactly_bipartite": 3717, "absolute_gap_below_0.01": len(near),
              "nonbipartite_below_0.01": len(near) - 3717,
              "large_relaxation_above_0.15_and_small_absolute_below_0.01": sum(
                  case["gap"] > .15 and case["absolute"] < .01 for case in cases),
              "certified_families": quotients}
    rerun = OUT / "target55-rerun/GRAPH_DIAGNOSTICS.tsv"
    if rerun.exists():
        new_p4 = 0
        for row in rows(rerun):
            if row["source"] != "case":
                continue
            flat = [int(x) for x in row["canonical_labeled_code"].split(",")]
            new_p4 += p4_fibers(np.array(flat, dtype=int).reshape(-1, 6)) is not None
        assert new_p4 == 49
        degree4 = ROOT / "runs/spectra/degrees-4-complete-20260927"
        edge_map: dict[str, dict[tuple[int, int], int]] = {}
        for row in rows(degree4 / "EDGES.tsv"):
            edge_map.setdefault(row["case_id"], {})[int(row["from"]), int(row["move"])] = int(row["to"])
        small_p4 = 0
        for row in rows(degree4 / "SPECTRUM.tsv"):
            if not .18 <= float(row["normalized_gap"]) < .20:
                continue
            n = int(row["quotient_states"])
            transitions = np.array([[edge_map[row["case_id"]][vertex, move] for move in range(6)]
                                    for vertex in range(n)], dtype=int)
            small_p4 += p4_fibers(transitions) is not None
        assert small_p4 == 10
        report["verified_four_position_quotient_cases_including_degree_four"] = (
            sum(int(item["cases"]) for item in quotients if item["quotient"] == "labeled_P4")
            + new_p4 + small_p4)
        assert report["verified_four_position_quotient_cases_including_degree_four"] == 2734
    (OUT / "RESULTS.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    with mpl.rc_context({"font.size": 9, "pdf.fonttype": 42}):
        fig, (ax, detail) = plt.subplots(1, 2, figsize=(10.5, 4.2), layout="constrained")
        for group, marker in (("A", "o"), ("S", "s")):
            for bipartite, color in ((False, "#0072B2"), (True, "#D55E00")):
                selected = [c for c in cases if c["group"].startswith(group) and c["bipartite"] == bipartite]
                for panel in (ax, detail):
                    panel.scatter([c["gap"] for c in selected], [c["absolute"] for c in selected],
                                  s=8, marker=marker, color=color, alpha=.18, rasterized=True,
                                  label=f"{group}, {'bipartite' if bipartite else 'non-bipartite'}" if panel is ax else None)
        for panel in (ax, detail):
            panel.plot([0, .42], [0, .42], color="#444444", lw=.8, ls="--", label="equal gaps" if panel is ax else None)
            panel.set(xlabel="Relaxation gap $1-\\lambda_2/6$", ylabel="Absolute gap")
        ax.set(xlim=(0, .42), ylim=(-.008, .42), title="All 7,319 complete-census cases")
        detail.set(xlim=(0, .25), ylim=(-.002, .065), title="Small absolute gaps (zoom)")
        ax.legend(loc="upper left", fontsize=7, frameon=True)
        fig.suptitle("Exact bipartiteness explains the zero-gap shelf; not every small gap is periodic")
        for suffix in ("png", "pdf"):
            fig.savefig(OUT / f"PAIRED_GAPS.{suffix}", dpi=300, facecolor="white")
        plt.close(fig)
    print(json.dumps({key: value for key, value in report.items() if key not in ("other_band_values", "certified_families")}, indent=2))


if __name__ == "__main__":
    main()
