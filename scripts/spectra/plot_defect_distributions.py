"""Plot exact defect histograms from defect-census TSV (no sampling)."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "runs/defect/20260928/DEFECT_DISTRIBUTIONS.tsv"
OUT = ROOT / "data/figures"

rows = list(csv.DictReader(DATA.open(), delimiter="\t"))
summary = list(csv.DictReader((ROOT / "evidence/defect/CENSUS_SUMMARY.tsv").open(), delimiter="\t"))
by_case = defaultdict(list)
for row in rows:
    by_case[row["case_id"]].append(row)

selected = [
    ("S6:component-54", r"$S_6$, $g=2$, identity entry"),
    ("S6:component-128", r"$S_6$, $g=2$, mixed"),
    ("S:CERTIFIED_S18_genus1_weak", r"$S_{18}$, $g=1$, always 1"),
    ("A:CERTIFIED_A25_order6_genus12_weak", r"$A_{25}$, $g=12$, always 6"),
    ("A:CERTIFIED_A25_order12_genus12_weak", r"$A_{25}$, $g=12$, always 8"),
]

with plt.rc_context({"pdf.fonttype": 42, "font.size": 10}):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8), layout="constrained",
                                   gridspec_kw={"width_ratios": [1.35, 1]})
    for i, (case, label) in enumerate(selected):
        rs = by_case[case]
        n = int(rs[0]["states"])
        for row in rs:
            defect, count = int(row["defect"]), int(row["count"])
            ax.barh(i, count / n, left=sum(int(r["count"]) for r in rs if int(r["defect"]) < defect) / n,
                    height=.7, color=plt.cm.Blues(.30 + .65 * defect / 8), edgecolor="white")
            if count / n > .07:
                ax.text((sum(int(r["count"]) for r in rs if int(r["defect"]) < defect) + count / 2) / n,
                        i, f"{defect}: {count}/{n}", ha="center", va="center", fontsize=8,
                        color="white" if defect >= 5 else "black")
    ax.set_yticks(range(len(selected)), [label for _, label in selected])
    ax.invert_yaxis()
    ax.set(xlim=(0, 1), xlabel="Fraction of states in the braid component", title="Exact defect distributions")
    ax.spines[["top", "right"]].set_visible(False)

    for j, group in enumerate(("S4", "S5", "S6")):
        for genus in range(1, 6):
            matches = [r for r in summary if r["group"] == group and int(r["genus"]) == genus]
            if not matches:
                continue
            all_zero, total = int(matches[0]["all_zero_components"]), int(matches[0]["all_components"])
            ax2.bar(genus + (j - 1) * .23, all_zero / total, width=.22,
                    label=group if genus == 1 else None, color=("#0072B2", "#E69F00", "#009E73")[j])
            ax2.text(genus + (j - 1) * .23, all_zero / total + .018,
                     f"{all_zero}/{total}", ha="center", va="bottom", fontsize=7, rotation=90)
    ax2.set(xlim=(.55, 5.45), ylim=(0, .66), xticks=range(1, 6), xlabel="Genus",
            ylabel="Fraction of components with defect always 0",
            title=r"Complete $S_4$--$S_6$ census ($g>0$)")
    ax2.legend(frameon=False, ncol=3, loc="upper right")
    ax2.text(.98, .82, "Every all-zero case shown here\nhas an identity entry.",
             transform=ax2.transAxes, ha="right", va="top", fontsize=9)
    ax2.spines[["top", "right"]].set_visible(False)
    identity_free_total = sum(int(r["identity_free_components"]) for r in summary)
    identity_free_zero = sum(int(r["identity_free_all_zero"]) for r in summary)
    fig.suptitle(f"Exact defect census: identity-free g>0 all-zero = {identity_free_zero}/{identity_free_total} in S4–S7", fontsize=12)
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png"):
        fig.savefig(OUT / f"defect_distributions_20260928.{extension}", dpi=250)
    plt.close(fig)
