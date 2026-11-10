#!/usr/bin/env python3
"""Recheck the A25 graph bridge and the complete A5/S5 spectral census."""

from __future__ import annotations

import csv
import hashlib
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "evidence/spectral_census"
A25 = ROOT / "runs/spectra/a25-current-full-20260925"
SMALL = ROOT / "runs/spectra/a5s5-canonical-20260925"
EXACT = ROOT / "runs/spectra/exact-certificates-20260924-final/rerun"
FULL = ROOT / "runs/spectra/degrees-5-8-full-20260922"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], data: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def cycle_type(perm: list[int]) -> tuple[int, ...]:
    assert set(perm) == set(range(len(perm)))
    unseen = set(range(len(perm)))
    lengths = []
    while unseen:
        point = next(iter(unseen))
        length = 0
        while point in unseen:
            unseen.remove(point)
            point = perm[point]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def pair_orbits(a: list[int], b: list[int]) -> tuple[int, ...]:
    unseen = set(range(len(a)))
    sizes = []
    while unseen:
        found = {next(iter(unseen))}
        frontier = list(found)
        while frontier:
            point = frontier.pop()
            for neighbor in (a[point], b[point]):
                if neighbor not in found:
                    found.add(neighbor)
                    frontier.append(neighbor)
        unseen.difference_update(found)
        sizes.append(len(found))
    return tuple(sorted(sizes, reverse=True))


def decode_a25(row: dict[str, str]) -> list[list[int]]:
    first = [[(int(row[f"key{j}"], 16) >> (5 * i)) & 31 for i in range(25)] for j in range(3)]
    for perm in first:
        assert set(perm) == set(range(25))
    product = [first[0][first[1][first[2][i]]] for i in range(25)]
    fourth = [0] * 25
    for point, image in enumerate(product):
        fourth[image] = point
    assert [product[fourth[i]] for i in range(25)] == list(range(25))
    return first + [fourth]


def a25_bridge() -> str:
    states = rows(A25 / "STATES.tsv")
    edges = rows(A25 / "EDGES.tsv")
    spectrum = rows(A25 / "SPECTRUM.tsv")
    diagnostics = rows(A25 / "GRAPH_DIAGNOSTICS.tsv")
    assert len(states) == 4 and len(spectrum) == len(diagnostics) == 1
    assert int(spectrum[0]["genus"]) == 0
    exceptional = (4, 4, 4, 4, 3, 3, 3)
    involution = (2,) * 10 + (1,) * 5
    by_position = {}
    output = []
    for row in states:
        perms = decode_a25(row)
        types = [cycle_type(perm) for perm in perms]
        assert sorted(types) == sorted([involution] * 3 + [exceptional])
        assert (2 - 2 * 25 + sum(25 - len(t) for t in types)) // 2 == 0
        position = types.index(exceptional) + 1
        state = int(row["state_index"])
        by_position[position] = state
        partition = pair_orbits(perms[0], perms[1])
        assert partition == ((15, 10) if position <= 2 else (12, 3, 3, 3, 2, 2))
        output.append({"state_index": state, "exceptional_position": position,
                       "first_pair_orbits": "+".join(map(str, partition)),
                       "first_pair_transitive": "false"})
    assert sorted(by_position) == [1, 2, 3, 4]
    transitions = {(int(row["from"]), int(row["move"])): int(row["to"]) for row in edges}
    assert len(transitions) == 24
    for row in output:
        state, position = row["state_index"], row["exceptional_position"]
        for move in range(6):
            i = move % 3 + 1
            destination = position + (1 if position == i else -1 if position == i + 1 else 0)
            assert transitions[state, move] == by_position[destination]
        row["move_targets_0_to_5"] = ",".join(str(transitions[state, move]) for move in range(6))
    write_tsv(EVIDENCE / "A25_FOUR_STATES.tsv", list(output[0]), sorted(output, key=lambda row: row["exceptional_position"]))

    code = diagnostics[0]["canonical_labeled_code"]
    exact_spectra = {row["case_id"]: row for row in rows(EXACT / "SPECTRUM.tsv")}
    comparators = []
    for row in rows(EXACT / "GRAPH_DIAGNOSTICS.tsv"):
        if row["canonical_labeled_code"] != code:
            continue
        case = exact_spectra[row["case_id"]]
        component = row["case_id"].rsplit("component-", 1)[1]
        source = Path(case["source"])
        with source.open(newline="", encoding="utf-8") as stream:
            raw = next(item for item in csv.DictReader(stream) if item["component"] == component)
        assert raw["weak"] == "true"
        comparators.append({"case_id": row["case_id"], "degree": case["degree"],
                            "genus": case["genus"], "weak": raw["weak"],
                            "passport": case["passport"]})
    assert {row["case_id"] for row in comparators} == {
        "A5:component-9", "A7:shard_1:component-20", "A7:shard_2:component-44"
    }
    assert sorted(int(row["genus"]) for row in comparators) == [0, 0, 2]
    comparators.append({"case_id": "A:CERTIFIED_A25_weak", "degree": "25", "genus": "0",
                        "weak": "false", "passport": spectrum[0]["passport"]})
    write_tsv(EVIDENCE / "A25_GRAPH_COMPARATORS.tsv", list(comparators[0]), comparators)
    return code


def a5_s5_spectra() -> None:
    diagnostics = {row["case_id"]: row for row in rows(SMALL / "GRAPH_DIAGNOSTICS.tsv") if row["source"] == "case"}
    cases = {row["case_id"]: row for row in rows(FULL / "SPECTRUM.tsv")
             if row["kind"] == "exact_component" and row["degree"] == "5"}
    assert len(diagnostics) == len(cases) == 118 and set(diagnostics) == set(cases)
    eigenvalues = defaultdict(list)
    with (FULL / "EIGENVALUES.tsv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["source"] == "case" and row["case_id"] in cases:
                eigenvalues[row["case_id"]].append((int(row["eigen_index"]), float(row["eigenvalue"])))
    by_code = defaultdict(list)
    for case_id, row in cases.items():
        values = eigenvalues[case_id]
        assert [i for i, _ in values] == list(range(int(row["quotient_states"])))
        by_code[diagnostics[case_id]["canonical_labeled_code"]].append(case_id)
    assert len(by_code) == 82

    fingerprint_to_code = {}
    class_rows = []
    for code, ids in by_code.items():
        fingerprints = {tuple(round(value, 6) for _, value in eigenvalues[case_id]) for case_id in ids}
        assert len(fingerprints) == 1, "an exact labeled graph class has inconsistent spectra"
        fingerprint = next(iter(fingerprints))
        assert fingerprint not in fingerprint_to_code, "distinct labeled classes are numerically cospectral"
        fingerprint_to_code[fingerprint] = code
        class_rows.append({"code_sha256_16": hashlib.sha256(code.encode()).hexdigest()[:16],
                           "vertices": len(fingerprint), "cases": len(ids),
                           "groups": ",".join(sorted({cases[case_id]["group"] + "5" for case_id in ids})),
                           "genera": ",".join(map(str, sorted({int(cases[case_id]["genus"]) for case_id in ids}))),
                           "case_ids": ",".join(sorted(ids))})
    cross = sorted((row for row in class_rows if "," in row["genera"]),
                   key=lambda row: (-row["cases"], row["code_sha256_16"]))
    assert len(cross) == 9 and sum(row["cases"] for row in cross) == 37
    assert sum("," in row["groups"] for row in class_rows) == 7
    write_tsv(EVIDENCE / "A5_S5_CROSS_GENUS_CLASSES.tsv", list(cross[0]), cross)

    summary = []
    for group in ("A", "S"):
        for genus in sorted({int(row["genus"]) for row in cases.values() if row["group"] == group}):
            case_ids = [case_id for case_id, row in cases.items()
                        if row["group"] == group and int(row["genus"]) == genus]
            representative = {diagnostics[case_id]["canonical_labeled_code"]: case_id for case_id in case_ids}
            gaps = [float(cases[case_id]["normalized_gap"]) for case_id in representative.values()]
            abs_gaps = [1 - float(cases[case_id]["abs_nontrivial_radius"]) / 6
                        for case_id in representative.values()]
            summary.append({"group": group + "5", "genus": genus, "raw_cases": len(case_ids),
                            "distinct_labeled_graphs": len(representative),
                            "median_relaxation_gap": f"{statistics.median(gaps):.8f}",
                            "mean_relaxation_gap": f"{statistics.mean(gaps):.8f}",
                            "mean_absolute_gap": f"{statistics.mean(abs_gaps):.8f}"})
    assert sum(row["raw_cases"] for row in summary) == 118
    write_tsv(EVIDENCE / "A5_S5_GENUS_SPECTRA.tsv", list(summary[0]), summary)
    print(f"A5/S5: {len(cases)} cases, {len(by_code)} exact labeled graph classes, "
          f"{len(cross)} cross-genus classes; all classes have distinct numerical full spectra at 1e-6")


if __name__ == "__main__":
    a25_bridge()
    a5_s5_spectra()
    print("A25: four states are the four exceptional-entry positions; all first pairs intransitive")
