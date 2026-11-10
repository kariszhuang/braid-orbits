#!/usr/bin/env python3
"""Reproducible matched-null and genus analysis for a Rust spectral census."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SEED = 20260922
COLORS = {5: "#0072B2", 6: "#E69F00", 7: "#009E73", 8: "#CC79A7"}
MARKERS = {"A": "o", "S": "^"}


def rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def floats(values):
    return np.asarray([float(value) for value in values if value not in ("", "nan")], dtype=float)


def bootstrap_ci(values, rng, reps=4000):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return np.nan, np.nan
    estimates = np.empty(reps)
    for i in range(reps):
        estimates[i] = np.mean(rng.choice(values, len(values), replace=True))
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def rankdata(values):
    values = np.asarray(values)
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2 + 1
        start = end
    return ranks


def stratified_genus_association(cases, rng, permutations=10000):
    strata = defaultdict(list)
    for case in cases:
        strata[(case["group"], case["degree"])].append(case)
    blocks = []
    for key, block in strata.items():
        genus = np.asarray([row["genus"] for row in block])
        gap = np.asarray([row["gap"] for row in block])
        if len(block) < 3 or len(np.unique(genus)) < 2:
            continue
        x = rankdata(genus); y = rankdata(gap)
        blocks.append((key, x - x.mean(), y - y.mean()))
    x = np.concatenate([block[1] for block in blocks])
    y = np.concatenate([block[2] for block in blocks])
    statistic = float(x @ y / np.sqrt((x @ x) * (y @ y)))
    extreme = 0
    for _ in range(permutations):
        yp = np.concatenate([rng.permutation(block[2]) for block in blocks])
        value = float(x @ yp / np.sqrt((x @ x) * (yp @ yp)))
        extreme += abs(value) >= abs(statistic)
    return statistic, (extreme + 1) / (permutations + 1), len(x), len(blocks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("census", type=Path)
    args = parser.parse_args()
    root = args.census.resolve()
    rng = np.random.default_rng(SEED)

    cases = []
    case_meta = {}
    for row in rows(root / "SPECTRUM.tsv"):
        if row["kind"] != "exact_component":
            continue
        item = {
            "id": row["case_id"], "group": row["group"], "degree": int(row["degree"]),
            "genus": int(row["genus"]), "vertices": int(row["quotient_states"]),
            "gap": float(row["normalized_gap"]), "lambda2": float(row["lambda2"]),
        }
        cases.append(item); case_meta[item["id"]] = item

    matched = []
    for row in rows(root / "MATCHED_BASELINES.tsv"):
        if row["case_id"] not in case_meta:
            continue
        item = dict(case_meta[row["case_id"]])
        item.update({
            "baseline": row["baseline"], "baseline_gap": float(row["mean_baseline_normalized_gap"]),
            "difference": float(row["actual_minus_baseline_normalized_gap"]),
            "replicates": int(row["replicates"]),
        })
        matched.append(item)

    summary_rows = []
    for baseline in ("random_6_regular", "random_schreier"):
        block = [row for row in matched if row["baseline"] == baseline]
        diff = np.asarray([row["difference"] for row in block])
        lo, hi = bootstrap_ci(diff, rng)
        summary_rows.append({
            "analysis": "matched_gap", "stratum": baseline, "n": len(diff),
            "estimate": np.mean(diff), "ci_low": lo, "ci_high": hi,
            "median": np.median(diff), "fraction_negative": np.mean(diff < 0),
        })

    rho, p_value, n_assoc, strata = stratified_genus_association(cases, rng)
    with (root / "ANALYSIS.tsv").open("w", newline="", encoding="utf-8") as handle:
        names = ["analysis", "stratum", "n", "estimate", "ci_low", "ci_high", "median", "fraction_negative"]
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t")
        writer.writeheader(); writer.writerows(summary_rows)
    (root / "GENUS_ASSOCIATION.tsv").write_text(
        "method\tn\tstrata\trho\tpermutation_p\tpermutations\tseed\n"
        f"within-group-degree Spearman\t{n_assoc}\t{strata}\t{rho:.10f}\t{p_value:.10g}\t10000\t{SEED}\n",
        encoding="utf-8",
    )

    mpl.rcParams.update({
        "font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 9,
        "legend.fontsize": 7.5, "pdf.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2), layout="constrained")

    ax = axes[0, 0]
    regular = [row for row in matched if row["baseline"] == "random_6_regular"]
    for group in ("A", "S"):
        for degree in sorted({row["degree"] for row in regular}):
            block = [row for row in regular if row["group"] == group and row["degree"] == degree]
            if block:
                ax.scatter([row["baseline_gap"] for row in block], [row["gap"] for row in block],
                           s=10, alpha=.35, color=COLORS.get(degree, "#666666"), marker=MARKERS[group],
                           label=f"{group}{degree}")
    limits = [0, max([row["baseline_gap"] for row in regular] + [row["gap"] for row in regular]) * 1.03]
    ax.plot(limits, limits, color="#333333", lw=.8, ls="--")
    ax.set(xlim=limits, ylim=limits, xlabel="Matched random 6-regular gap", ylabel="Hurwitz gap",
           title="A  Matched spectral gaps")
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles)); ax.legend(unique.values(), unique.keys(), ncol=2, frameon=False)

    ax = axes[0, 1]
    bins = np.linspace(-1, 1, 61)
    hist = {key: np.zeros(len(bins) - 1) for key in ("A5", "S5", "random_6_regular", "random_schreier")}
    for row in rows(root / "EIGENVALUES.tsv"):
        meta = case_meta.get(row["case_id"])
        if not meta or meta["degree"] != 5 or int(row["eigen_index"]) == 0:
            continue
        source = row["source"]
        key = meta["group"] + "5" if source == "case" else source.rsplit("-", 1)[0]
        if key in hist:
            value = float(row["eigenvalue"]) / 6
            index = min(max(np.searchsorted(bins, value, side="right") - 1, 0), len(bins) - 2)
            hist[key][index] += 1
    centers = (bins[:-1] + bins[1:]) / 2
    styles = {"A5": ("#0072B2", "-"), "S5": ("#D55E00", "-"),
              "random_6_regular": ("#000000", "--"), "random_schreier": ("#777777", ":")}
    for key, counts in hist.items():
        if counts.sum():
            ax.plot(centers, counts / counts.sum() / np.diff(bins), color=styles[key][0], ls=styles[key][1], label=key.replace("_", " "))
    ax.set(xlabel=r"Nontrivial normalized eigenvalue $\lambda/6$", ylabel="Pooled density",
           title="B  A5/S5 eigenvalue distributions")
    ax.legend(frameon=False)

    for ax, group, letter in ((axes[1, 0], "A", "C"), (axes[1, 1], "S", "D")):
        group_cases = [row for row in cases if row["group"] == group]
        for degree in sorted({row["degree"] for row in group_cases}):
            block = [row for row in group_cases if row["degree"] == degree]
            jitter = rng.uniform(-.10, .10, len(block))
            ax.scatter(np.asarray([row["genus"] for row in block]) + jitter, [row["gap"] for row in block],
                       s=8, alpha=.18, color=COLORS.get(degree, "#666666"), marker=MARKERS[group])
            medians = defaultdict(list)
            for row in block: medians[row["genus"]].append(row["gap"])
            x = sorted(medians); y = [np.median(medians[value]) for value in x]
            ax.plot(x, y, color=COLORS.get(degree, "#666666"), marker="o", ms=2.5, lw=1, label=f"{group}{degree}")
        ax.set(xlabel="Genus", ylabel="Normalized spectral gap", title=f"{letter}  Genus and gap in {group}d")
        handles, labels = ax.get_legend_handles_labels()
        if handles: ax.legend(handles, labels, frameon=False, ncol=2)

    fig.suptitle("Hurwitz Cayley–Schreier spectra (six labeled braid moves)", fontsize=12)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(root / f"spectral-overview.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)

    regular_summary = next(row for row in summary_rows if row["stratum"] == "random_6_regular")
    schreier_summary = next(row for row in summary_rows if row["stratum"] == "random_schreier")
    degrees = sorted({(row["group"], row["degree"]) for row in cases})
    text = f"""# Spectral analysis

Analyzed {len(cases):,} exact quotient components across {', '.join(f'{g}{d}' for g, d in degrees)}. The normalized gap is `(6 - lambda2) / 6`; smaller values mix more slowly.

Against same-size random 6-regular graphs, the mean paired gap difference (Hurwitz minus null) was {regular_summary['estimate']:.4f} (component bootstrap 95% interval {regular_summary['ci_low']:.4f} to {regular_summary['ci_high']:.4f}); {regular_summary['fraction_negative']:.1%} of components had the smaller gap. Against random Schreier graphs the corresponding difference was {schreier_summary['estimate']:.4f} ({schreier_summary['ci_low']:.4f} to {schreier_summary['ci_high']:.4f}), with {schreier_summary['fraction_negative']:.1%} negative.

The exploratory within-group-and-degree Spearman association between genus and normalized gap was rho={rho:.4f}, stratified permutation p={p_value:.4g} (n={n_assoc:,} components, {strata} strata, 10,000 permutations). This controls group and degree but not passport or component size, so it is evidence of association rather than a causal genus effect.

Full eigenvalue-distribution curves use nontrivial eigenvalues from exact A5/S5 spectra and same-size null graphs. Larger components contribute only extremal Lanczos estimates and therefore enter gap/genus panels, not the full-distribution panel. Random seed: {SEED}.
"""
    (root / "ANALYSIS.md").write_text(text, encoding="utf-8")
    (root / "FIGURE_ALT_TEXT.txt").write_text(
        "Four-panel figure comparing Hurwitz graph gaps with same-size random 6-regular graphs, "
        "pooled nontrivial A5/S5 eigenvalue densities with random regular and random Schreier nulls, "
        "and genus-versus-gap point clouds with degree-specific median lines for alternating and symmetric groups.\n",
        encoding="utf-8",
    )
    print(root / "spectral-overview.png")


if __name__ == "__main__":
    main()
