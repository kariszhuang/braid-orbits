#!/usr/bin/env python3
"""Exact characteristic-polynomial proof for the 55 larger and 13 degree-four 0.195262 cases."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
TARGET55 = ROOT / "runs/spectra/gap-band-bipartite-20260927/target55-rerun"
DEGREE4 = ROOT / "runs/spectra/degrees-4-complete-20260927"
OUTPUT = ROOT / "evidence/spectral_census/GAP_BAND_ADDITIONAL_CERTIFICATES.tsv"
X = sp.Symbol("x")
TARGET_FACTOR = sp.Poly(X * X - 4 * X - 4, X)
LOWER = sp.Rational(482842712, 100000000)
assert LOWER**2 - 4 * LOWER - 4 < 0


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def certificate(transitions: list[list[int]], cases: list[str], kind: str) -> dict[str, object]:
    n = len(transitions)
    assert all(len(row) == 6 and all(0 <= j < n for j in row) for row in transitions)
    matrix = sp.zeros(n)
    for vertex, row in enumerate(transitions):
        for move, neighbor in enumerate(row):
            assert transitions[neighbor][(move + 3) % 6] == vertex
            matrix[vertex, neighbor] += 1
    assert matrix == matrix.T and all(sum(matrix.row(i)) == 6 for i in range(n))
    factors = sp.factor_list(sp.Poly(matrix.charpoly(X).as_expr(), X))[1]
    six = [(factor, multiplicity) for factor, multiplicity in factors if factor == sp.Poly(X - 6, X)]
    target = [(factor, multiplicity) for factor, multiplicity in factors if factor == TARGET_FACTOR]
    assert len(six) == len(target) == 1 and six[0][1] == 1
    assert all(factor.count_roots(LOWER, sp.Integer(6)) == 0
               for factor, _ in factors if factor not in (six[0][0], TARGET_FACTOR))
    return {"kind": kind, "representative": cases[0], "case_count": len(cases),
            "case_ids": ",".join(cases), "vertices": n,
            "target_multiplicity": target[0][1],
            "factorization": "*".join(f"({sp.sstr(factor.as_expr())})^{multiplicity}"
                                     for factor, multiplicity in factors)}


def main() -> None:
    large = defaultdict(list)
    for row in rows(TARGET55 / "GRAPH_DIAGNOSTICS.tsv"):
        if row["source"] == "case":
            large[row["canonical_labeled_code"]].append(row["case_id"])
    assert len(large) == 21 and sum(map(len, large.values())) == 55
    result = []
    for code, cases in large.items():
        flat = [int(x) for x in code.split(",")]
        assert len(flat) % 6 == 0
        result.append(certificate([flat[i:i + 6] for i in range(0, len(flat), 6)], sorted(cases), "large"))

    spectra = rows(DEGREE4 / "SPECTRUM.tsv")
    edges = defaultdict(dict)
    for row in rows(DEGREE4 / "EDGES.tsv"):
        edges[row["case_id"]][(int(row["from"]), int(row["move"]))] = int(row["to"])
    target_value = (2 - 2**0.5) / 3
    small = [row for row in spectra if .18 <= float(row["normalized_gap"]) < .20]
    assert len(small) == 13 and all(abs(float(row["normalized_gap"]) - target_value) < 1e-8 for row in small)
    for row in small:
        case_id, n = row["case_id"], int(row["quotient_states"])
        assert len(edges[case_id]) == 6 * n
        transitions = [[edges[case_id][vertex, move] for move in range(6)] for vertex in range(n)]
        result.append(certificate(transitions, [case_id], "degree_four"))
    assert len(result) == 34 and sum(int(row["case_count"]) for row in result) == 68
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    prior = rows(ROOT / "evidence/spectral_census/EXACT_LAMBDA2_CERTIFICATES.tsv")
    assert len(prior) == 8 and sum(int(row["case_count"]) for row in prior) == 2757
    characteristic_polynomials = {
        tuple(sp.Poly(sp.sympify(row["characteristic_polynomial_factors"], locals={"x": X}), X).all_coeffs())
        for row in prior
    }
    characteristic_polynomials.update(
        tuple(sp.Poly(sp.sympify(row["factorization"].replace("^", "**"), locals={"x": X}), X).all_coeffs())
        for row in result
    )
    assert len(characteristic_polynomials) == 31
    print("Certified 68 additional cases in 34 exact graph computations (55 large, 13 degree-four).")
    print("Together with 2,757 prior certificates: 2,825 cases with exact gap (2-sqrt(2))/3.")
    print("These cases have 31 distinct exact characteristic polynomials, so their full spectra are not all the same.")


if __name__ == "__main__":
    main()
