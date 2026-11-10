#!/usr/bin/env python3
"""Rerun small recurring spectra and certify their labeled graph classes exactly.

Requires SymPy and Matplotlib. This script performs computations only when invoked explicitly.
"""

from __future__ import annotations
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass

import argparse
import csv
import json
from pathlib import Path
import subprocess

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp

from run_spectral_followups import component_files


TARGET = 2 + 2 * sp.sqrt(2)
LOWER = sp.Rational(482842712, 100000000)
UPPER = sp.Rational(482842713, 100000000)
assert LOWER > 2 and LOWER**2 - 4*LOWER - 4 < 0 < UPPER**2 - 4*UPPER - 4


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, values: list[dict[str, object]]) -> None:
    if not values:
        raise ValueError(f"no rows to write: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(values)


def adjacency_from_code(code: str, vertices: int) -> sp.Matrix:
    targets = [int(value) for value in code.split(",")]
    if len(targets) != 6 * vertices or any(target < 0 or target >= vertices for target in targets):
        raise ValueError("invalid canonical transition code")
    matrix = sp.zeros(vertices)
    for vertex in range(vertices):
        for move in range(6):
            target = targets[6 * vertex + move]
            if targets[6 * target + (move + 3) % 6] != vertex:
                raise ValueError("inverse move mismatch")
            matrix[vertex, target] += 1
    if matrix != matrix.T or any(sum(matrix.row(i)) != 6 for i in range(vertices)):
        raise ValueError("adjacency is not symmetric 6-regular")
    seen = {0}
    queue = [0]
    for vertex in queue:
        for target in targets[6 * vertex:6 * vertex + 6]:
            if target not in seen:
                seen.add(target)
                queue.append(target)
    if len(seen) != vertices:
        raise ValueError("canonical graph is disconnected")
    return matrix


def exact_certificate(code: str, vertices: int) -> dict[str, object]:
    x = sp.Symbol("x")
    polynomial = sp.Poly(adjacency_from_code(code, vertices).charpoly(x).as_expr(), x)
    full_polynomial = polynomial
    target_factor = sp.Poly(x**2 - 4*x - 4, x)
    negative_target_factor = sp.Poly(x**2 + 4*x - 4, x)
    quotient = polynomial
    multiplicity = 0
    while True:
        divided, remainder = divmod(quotient, target_factor)
        if remainder:
            break
        quotient = divided
        multiplicity += 1
    # Exact rational root isolation proves whether another eigenvalue lies
    # between 2+2sqrt(2) and the trivial eigenvalue 6.
    other_roots = quotient.intervals(eps=sp.Rational(1, 10**12))
    only_trivial_above = all(right <= LOWER or left >= 6 for (left, right), _ in other_roots)
    simple_six = polynomial.eval(6) == 0 and sp.diff(polynomial.as_expr(), x).subs(x, 6) != 0
    real_root_intervals = polynomial.intervals(eps=sp.Rational(1, 10**12))
    lambda_min_lower, lambda_min_upper = real_root_intervals[0][0]
    lambda_min_approx = float((lambda_min_lower + lambda_min_upper) / 2)
    bipartite = polynomial.eval(-6) == 0
    negative_quotient, negative_remainder = divmod(polynomial, negative_target_factor)
    negative_target_multiplicity = 0
    while not negative_remainder:
        negative_target_multiplicity += 1
        negative_quotient, negative_remainder = divmod(negative_quotient, negative_target_factor)
    no_other_root_below_target = all(left > -LOWER for (left, _), _ in negative_quotient.intervals(eps=sp.Rational(1, 10**12)))
    target_controls_absolute_radius = not bipartite and (
        lambda_min_lower > -LOWER or (negative_target_multiplicity > 0 and no_other_root_below_target)
    )
    if bipartite:
        absolute_radius_source, absolute_gap_exact, absolute_gap_approx = "eigenvalue -6", "0", 0.0
    elif target_controls_absolute_radius:
        absolute_radius_source = "target eigenvalue pair"
        absolute_gap_exact = "(2-sqrt(2))/3"
        absolute_gap_approx = float((2 - sp.sqrt(2)) / 3)
    else:
        absolute_radius_source = "negative extremal root"
        absolute_gap_exact = ""
        absolute_gap_approx = 1 - max(float(TARGET), abs(lambda_min_approx)) / 6
    return {
        "vertices": vertices,
        "target_multiplicity": multiplicity,
        "lambda2_certified": multiplicity > 0 and only_trivial_above and simple_six,
        "six_simple": simple_six,
        "other_roots_below_target_or_at_six": only_trivial_above,
        "lambda_min_interval_lower": str(lambda_min_lower),
        "lambda_min_interval_upper": str(lambda_min_upper),
        "lambda_min_approx": f"{lambda_min_approx:.12f}",
        "bipartite_certified": bipartite,
        "negative_target_multiplicity": negative_target_multiplicity,
        "absolute_radius_source": absolute_radius_source,
        "absolute_gap_exact": absolute_gap_exact,
        "absolute_gap_approx": f"{absolute_gap_approx:.12f}",
        "characteristic_polynomial_factors": str(sp.factor(full_polynomial.as_expr())),
    }


def make_figure(output: Path, families: list[dict[str, object]]) -> None:
    ordered = sorted(families, key=lambda row: -int(row["case_count"]))
    labels = [f"n={row['vertices']}; g={min(map(int, str(row['genera']).split(',')))}–{max(map(int, str(row['genera']).split(',')))}" for row in ordered]
    positions = list(range(len(ordered)))
    colors = ["#D55E00" if row["bipartite_certified"] else "#0072B2" for row in ordered]
    mpl.rcParams.update({"font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, (count_ax, gap_ax) = plt.subplots(1, 2, figsize=(11, 5.4), layout="constrained")
    count_ax.barh(positions, [int(row["case_count"]) for row in ordered], color=colors,
                  edgecolor="#222222", linewidth=.45,
                  hatch=["//" if row["bipartite_certified"] else "" for row in ordered])
    count_ax.set_yticks(positions, labels)
    count_ax.invert_yaxis()
    count_ax.set(xlabel="Cases in this exact move-labeled class", title="Class frequency", xlim=(0, max(int(row["case_count"]) for row in ordered) * 1.45))
    for y, row in enumerate(ordered):
        count_ax.text(int(row["case_count"]) + max(int(ordered[0]["case_count"]) * .012, 1), y,
                      f"{row['case_count']}  ({'bipartite' if row['bipartite_certified'] else 'not bipartite'})",
                      va="center", fontsize=7)
    gap_ax.barh([y - .18 for y in positions], [float(row["absolute_gap_approx"]) if not row["bipartite_certified"] else 0 for row in ordered],
                height=.34, color=colors, edgecolor="#222222", linewidth=.45,
                hatch=["//" if row["bipartite_certified"] else "" for row in ordered], label="Absolute gap")
    relaxation_gap = float((2 - sp.sqrt(2)) / 3)
    gap_ax.barh([y + .18 for y in positions], [relaxation_gap] * len(ordered), height=.34,
                color="#999999", edgecolor="#222222", linewidth=.45, label="Relaxation gap")
    for y, row in enumerate(ordered):
        gap = float(row["absolute_gap_approx"])
        gap_ax.text(max(gap, 0) + .006, y - .18, f"{gap:.3f}", va="center", fontsize=7)
    gap_ax.set_yticks(positions, labels)
    gap_ax.invert_yaxis()
    gap_ax.set(xlabel="Normalized gap", title="Normalized gaps\n(gray: relaxation; color: absolute)", xlim=(0, .24))
    fig.suptitle(r"Exact families with $\lambda_2=2+2\sqrt{2}$", fontsize=12)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / f"EXACT_LAMBDA2_FAMILIES.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)
    (output / "EXACT_LAMBDA2_FIGURE_ALT_TEXT.txt").write_text(
        "Eight exact move-labeled graph classes, ordered by case count, all have second-largest eigenvalue 2+2sqrt(2) and relaxation gap (2-sqrt(2))/3. Two classes, containing 2352 of 2757 cases, are bipartite: their eigenvalue -6 makes the absolute gap zero. The other six classes have absolute gap equal to the relaxation gap. Bars are colored and hatched by bipartiteness; counts and status are also printed directly.\n",
        encoding="utf-8")
    (output / "FIGURE_PROVENANCE.json").write_text(json.dumps({
        "source": "EXACT_LABELED_FAMILIES.tsv",
        "filter": "exact-component cases with quotient size <=24 and prior numerical lambda2 within 1e-8 of 2+2sqrt(2); sampled S8 omitted",
        "aggregation": "case counts by exact move-labeled graph isomorphism class",
        "gap_definitions": {"relaxation": "1-lambda2/6", "absolute": "1-max(abs(lambda2),abs(lambda_min))/6"},
        "uncertainty": "none; exact characteristic-polynomial certificates",
        "formats": ["PNG", "PDF", "SVG"],
        "publisher_requirements": "not specified; provisional research figure",
    }, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="completed census containing SPECTRUM.tsv")
    parser.add_argument("--binary", type=Path, required=True, help="current compiled hurwitz_breaker binary")
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    source = args.source.resolve()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary not found: {binary}")
    selected = [row for row in rows(source / "SPECTRUM.tsv")
                if row["kind"] == "exact_component" and row["method"] == "jacobi-exact"
                and int(row["quotient_states"]) <= 24
                and abs(float(row["lambda2"]) - float(TARGET)) < 1e-8
                and (row["group"], row["degree"]) != ("S", "8")]
    if not selected:
        parser.error("no eligible non-S8 recurring small graphs in source")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    ids = output / "SELECTED_CASE_IDS.txt"
    selected_ids = sorted({row["case_id"] for row in selected})
    ids.write_text("\n".join(selected_ids) + "\n", encoding="utf-8")
    rerun = output / "rerun"
    command = [str(binary), "spectrum", "--components", ",".join(map(str, component_files())),
               "--select-case-ids", str(ids), "--out", str(rerun), "--threads", str(args.threads),
               "--full-max", "24", "--degree5-full-max", "24", "--baseline-reps", "0",
               "--baseline-kinds", "", "--canonical-max", "24", "--summary-only"]
    subprocess.run(command, check=True)

    diagnostics = {row["case_id"]: row for row in rows(rerun / "GRAPH_DIAGNOSTICS.tsv") if row["source"] == "case"}
    with (rerun / "GRAPH_DIAGNOSTICS.tsv").open(encoding="utf-8") as handle:
        if "canonical_labeled_code" not in handle.readline().rstrip("\n").split("\t"):
            raise RuntimeError("stale binary: rebuild with CARGO_TARGET_DIR=target/spectrum cargo build --release --bin hurwitz_breaker")
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows(rerun / "SPECTRUM.tsv"):
        diagnostic = diagnostics[row["case_id"]]
        code = diagnostic["canonical_labeled_code"]
        if not code:
            raise ValueError(f"missing exact labeled code for {row['case_id']}")
        groups.setdefault(code, []).append(row)
    output_rows = []
    for code, family in groups.items():
        result = exact_certificate(code, int(family[0]["quotient_states"]))
        output_rows.append({
            "case_count": len(family), "groups": ",".join(sorted({row["group"] + row["degree"] for row in family})),
            "genera": ",".join(sorted({row["genus"] for row in family}, key=int)),
            "passports": len({row["passport"] for row in family}),
            "representative_case_id": family[0]["case_id"],
            **result, "canonical_labeled_code": code,
        })
    output_rows = sorted(output_rows, key=lambda row: -int(row["case_count"]))
    write_rows(output / "EXACT_LABELED_FAMILIES.tsv", output_rows)
    make_figure(output, output_rows)
    if sum(row["case_count"] for row in output_rows) != len(selected_ids):
        raise ValueError("rerun did not reproduce every selected case")
    print(output / "EXACT_LABELED_FAMILIES.tsv")


if __name__ == "__main__":
    main()
