#!/usr/bin/env python3
"""Plot every stored degree-5-to-8 quotient case, without a spectral filter."""

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs/spectra/structural-followups-20260922-ecf2d7a/full-census/SPECTRUM.tsv"
DEGREE4 = ROOT / "runs/spectra/degrees-4-complete-20260927/SPECTRUM.tsv"
S8_EXACT = ROOT / "runs/spectra/S8_exact_shards"
OUTPUT = ROOT / "data/figures"
GROUPS = (
    ("Complete degree 5: A5, S5", lambda r: r["degree"] == "5"),
    ("Complete A6–A8, S6–S7", lambda r: r["degree"] != "5" and not (r["group"] == "S" and r["degree"] == "8")),
    ("Sampled S8", lambda r: r["group"] == "S" and r["degree"] == "8"),
)


def main() -> None:
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 7343 and len({r["case_id"] for r in rows}) == len(rows)
    assert Counter(r["method"] for r in rows) == {"lanczos-192": 4032, "jacobi-exact": 3311}
    with DEGREE4.open(newline="", encoding="utf-8") as stream:
        degree4 = list(csv.DictReader(stream, delimiter="\t"))
    assert len(degree4) == 19
    exact_paths = [S8_EXACT / f"shard_{i}/SPECTRUM.tsv" for i in range(16)]
    s8_complete = all(path.is_file() for path in exact_paths)
    if s8_complete:
        exact_s8 = []
        for path in exact_paths:
            with path.open(newline="", encoding="utf-8") as stream:
                exact_s8.extend(csv.DictReader(stream, delimiter="\t"))
        assert len({r["case_id"] for r in exact_s8}) == len(exact_s8)
        assert all(r["group"] == "S" and r["degree"] == "8" for r in exact_s8)
        exhaustive_figures(degree4 + [r for r in rows if not (r["group"] == "S" and r["degree"] == "8")] + exact_s8,
                           exact_paths)
    else:
        exhaustive_figures(degree4 + rows, [])
    counts = []
    bins = np.arange(0, 0.3801, 0.02)
    grids = []
    for label, selector in GROUPS:
        selected = [r for r in rows if selector(r)]
        counts.append(len(selected))
        for gap_name in ("Relaxation gap", "Absolute gap"):
            gaps = [float(r["normalized_gap"]) if gap_name == "Relaxation gap"
                    else 1 - float(r["abs_nontrivial_radius"]) / 6 for r in selected]
            assert all(-1e-8 <= gap < bins[-1] for gap in gaps)
            # Values numerically equal to zero belong in the zero bin.
            gaps = np.maximum(gaps, 0)
            hist, _, _ = np.histogram2d([int(r["genus"]) for r in selected], gaps,
                                          bins=[np.arange(-.5, 7.5, 1), bins])
            assert int(hist.sum()) == len(selected)
            grids.append(hist.T)
    assert counts == [118, 7201, 24]

    cmap = mpl.colormaps["YlGnBu"].copy()
    cmap.set_bad("#F1F1F1")
    with mpl.rc_context({"font.size": 9, "pdf.fonttype": 42}):
        fig, axes = plt.subplots(2, 3, figsize=(11.3, 6.0), layout="constrained",
                                 sharex=True, sharey=True)
        for col, ((label, _), count) in enumerate(zip(GROUPS, counts)):
            axes[0, col].set_title(f"{label}\n{count:,} cases", fontsize=10)
            for row, gap_name in enumerate(("Relaxation gap", "Absolute gap")):
                ax = axes[row, col]
                grid = grids[2 * col + row]
                mesh = ax.pcolormesh(np.arange(-.5, 7.5), bins,
                                     np.ma.masked_where(grid == 0, grid),
                                     cmap=cmap, norm=LogNorm(vmin=1, vmax=1000),
                                     shading="flat", edgecolors="white", linewidth=.2)
                ax.set(xlim=(-.5, 6.5), ylim=(0, .38), xticks=range(7))
                ax.set_yticks(np.arange(0, .381, .04))
                if col == 0:
                    ax.set_ylabel(f"{gap_name}\n(normalized)")
                if row == 1:
                    ax.set_xlabel("Cover genus")
                if col != 1:
                    for y, x in zip(*np.nonzero(grid)):
                        value = int(grid[y, x])
                        if value:
                            ax.text(x, bins[y] + .01, str(value), ha="center", va="center",
                                    fontsize=5.6, color="white" if value >= 30 else "black")
        fig.colorbar(mesh, ax=axes, ticks=[1, 10, 100, 1000], label="Cases per genus × 0.02-gap bin",
                     shrink=.83, pad=.015)
        fig.suptitle("All stored small-degree braid graphs: genus versus spectral gap", fontsize=12)
        OUTPUT.mkdir(parents=True, exist_ok=True)
        for suffix in ("pdf", "png"):
            fig.savefig(OUTPUT / f"all_small_genus_gaps.{suffix}", dpi=300, facecolor="white")
        plt.close(fig)

    provenance = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "cases": dict(zip((label for label, _ in GROUPS), counts)),
        "gap_definitions": {"relaxation": "1-lambda2/6 (normalized_gap field)",
                            "absolute": "1-abs_nontrivial_radius/6"},
        "bins": "[0,0.02), [0.02,0.04), ..., [0.36,0.38); negative floating roundoff clipped to zero",
        "spectral_methods": {"jacobi-exact": 3311, "lanczos-192": 4032},
        "scope": "All 7343 source rows; no eigenvalue, graph-size, or genus filter. A5–A8 and S5–S7 are exhaustively enumerated; S8 is sampled.",
    }
    (OUTPUT / "all_small_genus_gaps.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (OUTPUT / "all_small_genus_gaps.alt.txt").write_text(
        "Six heatmaps show counts of all 7,343 stored small-degree cases by genus and normalized spectral gap, "
        "in 0.02-wide gap bins. Columns separate complete A5/S5 (118 cases), complete A6–A8/S6–S7 "
        "(7,201 cases), and sampled S8 (24 cases). Rows show relaxation and absolute gaps. "
        "Gray cells contain no cases. Many absolute gaps are zero, indicating bipartite graphs. "
        "The plotted numerical gaps include 4,032 Lanczos estimates; only S8 is sampled.\n"
    )


def exhaustive_figures(rows: list[dict[str, str]], s8_sources: list[Path]) -> None:
    """All complete small-degree cases, with sampled S8 clearly distinguished."""
    largest_gap = max(max(float(r["normalized_gap"]), 1 - float(r["abs_nontrivial_radius"]) / 6)
                      for r in rows)
    upper = max(.44, .02 * (math.ceil(largest_gap / .02) + 1))
    bins = np.arange(0, upper + .0001, .02)
    cmap = mpl.colormaps["YlGnBu"].copy()
    cmap.set_bad("#F1F1F1")
    expected = {"A": [6, 43, 272, 1061, 3834],
                "S": [13, 75, 350, 1684, len([r for r in rows if r["group"] == "S" and r["degree"] == "8"])]}
    for kind in ("A", "S"):
        degrees = tuple(range(4, 9))
        selected = [[r for r in rows if r["group"] == kind and r["degree"] == str(d)]
                    for d in degrees]
        assert [len(group) for group in selected] == expected[kind]
        with mpl.rc_context({"font.size": 9, "pdf.fonttype": 42}):
            fig, axes = plt.subplots(2, len(degrees), figsize=(3.0 * len(degrees) + .5, 5.8), layout="constrained",
                                     sharex=True, sharey=True)
            for col, (degree, group) in enumerate(zip(degrees, selected)):
                scope = " · sample" if kind == "S" and degree == 8 and not s8_sources else ""
                axes[0, col].set_title(f"{kind}{degree} · {len(group):,} cases{scope}", fontsize=11)
                for row, label in enumerate(("Relaxation gap", "Absolute gap")):
                    gap = [float(r["normalized_gap"]) if row == 0 else
                           1 - float(r["abs_nontrivial_radius"]) / 6 for r in group]
                    assert all(-1e-8 <= value < bins[-1] for value in gap)
                    hist, _, _ = np.histogram2d([int(r["genus"]) for r in group],
                                                  np.maximum(gap, 0),
                                                  bins=[np.arange(-.5, 8.5), bins])
                    assert int(hist.sum()) == len(group)
                    ax = axes[row, col]
                    mesh = ax.pcolormesh(np.arange(-.5, 8.5), bins,
                                         np.ma.masked_where(hist.T == 0, hist.T),
                                         cmap=cmap, norm=LogNorm(vmin=1, vmax=1000),
                                         shading="flat", edgecolors="white", linewidth=.2)
                    ax.set(xlim=(-.5, 7.5), ylim=(0, upper), xticks=range(8))
                    ax.set_yticks(np.arange(0, upper + .001, .04))
                    if col == 0:
                        ax.set_ylabel(f"{label}\n(normalized)")
                    if row == 1:
                        ax.set_xlabel("Cover genus")
                    if col == 0:
                        for y, x in zip(*np.nonzero(hist.T)):
                            value = int(hist.T[y, x])
                            ax.text(x, bins[y] + .01, str(value), ha="center", va="center",
                                    fontsize=5.7, color="white" if value >= 30 else "black")
            fig.colorbar(mesh, ax=axes, ticks=[1, 10, 100, 1000],
                         label="Cases per genus × 0.02-gap bin", shrink=.84, pad=.015)
            title = ("Complete A4–A8 quotient census: genus versus spectral gap" if kind == "A" else
                     "Complete S4–S8 quotient census: genus versus spectral gap" if s8_sources else
                     "S4–S7 complete; S8 sampled: genus versus spectral gap")
            fig.suptitle(title,
                         fontsize=12)
            OUTPUT.mkdir(parents=True, exist_ok=True)
            stem = ("exhaustive_A4_A8_genus_gaps" if kind == "A" else
                    "exhaustive_S4_S8_genus_gaps" if s8_sources else "S4_S8_sampled_S8_genus_gaps")
            for suffix in ("png", "pdf"):
                fig.savefig(OUTPUT / f"{stem}.{suffix}", dpi=300, facecolor="white")
            plt.close(fig)
            (OUTPUT / f"{stem}.alt.txt").write_text(
                f"{2 * len(degrees)} heatmaps show {kind}4–{kind}{degrees[-1]} generating components up to simultaneous "
                f"S_d relabeling: case counts {expected[kind]}. Columns are degrees {degrees}. "
                + ("S8 has only 24 sampled cases; S4–S7 are complete. " if kind == "S" and not s8_sources else
                   "All degrees are complete. ") +
                "top row relaxation gap, bottom row absolute gap. Horizontal axis is cover genus "
                "0–7; vertical axis is normalized gap in bins of width 0.02. Gray means no cases. "
                "Gaps for graphs larger than the full-spectrum threshold are numerical Lanczos estimates.\n"
            )
            (OUTPUT / f"{stem}.provenance.json").write_text(json.dumps({
                "sources": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in ((DEGREE4, SOURCE) + tuple(s8_sources if kind == "S" else []))},
                "selection": ("Every generating component of A4–A8, modulo simultaneous S_d relabeling" if kind == "A" else
                              "Every generating component of S4–S8, modulo simultaneous S_d relabeling" if s8_sources else
                              "Every generating component of S4–S7; 24 selected S8 cases only; modulo simultaneous S_d relabeling"),
                "cases_by_degree": dict(zip(degrees, expected[kind])),
                "bin_width": .02,
                "spectral_methods": dict(Counter(r["method"] for group in selected for r in group)),
            }, indent=2) + "\n")


if __name__ == "__main__":
    main()
