"""Fast exact comparison of saved component histograms with the proposed coins.

Read-only: prints JSON; does not launch searches or modify saved certificates.
Primary population: identity-free complete S4--S7 components, positive genus.
No iid p-values: states form finite enumerated components, not random samples.
"""
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def coins(genus):
    pmf = [Fraction(1)]
    for j in range(1, genus + 1):
        nxt = [Fraction(0)] * (len(pmf) + 1)
        for k, value in enumerate(pmf):
            nxt[k] += value * Fraction(j, j + 1)
            nxt[k + 1] += value * Fraction(1, j + 1)
        pmf = nxt
    assert sum(pmf) == 1
    assert pmf[0] == Fraction(1, genus + 1)
    assert pmf[-1] == Fraction(1, math.factorial(genus + 1))
    return pmf


def metrics(prob, model):
    mean = sum(k * p for k, p in enumerate(prob))
    return {
        "observed": [float(p) for p in prob],
        "model": [float(p) for p in model],
        "tv": float(sum(abs(a-b) for a, b in zip(prob, model))/2),
        "mean": float(mean),
        "model_mean": float(sum(k*p for k, p in enumerate(model))),
    }


def main():
    files = [ROOT / "runs/defect/20260928/DEFECT_DISTRIBUTIONS.tsv",
             ROOT / "runs/defect/20260928_S7_no_identity/DEFECT_DISTRIBUTIONS.tsv"]
    cases = {}
    hashes = {}
    for path in files:
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                case = cases.setdefault(row["case_id"], {"row": row, "hist": Counter()})
                case["hist"][int(row["defect"])] += int(row["count"])
    sources = {}
    groups = defaultdict(list)
    individual = []
    for name, case in cases.items():
        row, hist = case["row"], case["hist"]
        genus, states = int(row["genus"]), int(row["states"])
        assert sum(hist.values()) == states
        assert all(0 <= k <= genus for k in hist)
        if genus == 0 or not name.startswith(("S4:", "S5:", "S6:", "S7:")):
            continue
        source = ROOT / row["source"]
        if source not in sources:
            with source.open() as handle:
                sources[source] = {r["component"]: r for r in csv.DictReader(handle)}
            hashes[str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
        original = sources[source][name.rsplit("-", 1)[1]]
        assert int(original["quotient_states"]) == states
        identity = "()" in original["seed"]
        degree = int(row["degree"])
        prob = [Fraction(hist[k], states) for k in range(genus+1)]
        detail = {"case": name, "degree": degree, "genus": genus,
                  "states": states, "identity": identity, **metrics(prob, coins(genus))}
        individual.append(detail)
        groups[(degree, genus, "available_saved_including_identity")].append(detail)
        if not identity:
            groups[(degree, genus, "identity_free")].append(detail)
    summaries = []
    for (degree, genus, scope), items in sorted(groups.items()):
        total = sum(i["states"] for i in items)
        for weighting in ("state_weighted", "equal_component"):
            weights = [Fraction(i["states"], total) if weighting == "state_weighted"
                       else Fraction(1, len(items)) for i in items]
            # Reconstruct probabilities from integer histograms, not floats.
            prob = [sum(w * Fraction(cases[i["case"]]["hist"][k], i["states"])
                        for w, i in zip(weights, items)) for k in range(genus+1)]
            summaries.append({"degree": degree, "genus": genus, "scope": scope,
                              "weighting": weighting, "components": len(items),
                              "states": total, **metrics(prob, coins(genus))})
    primary = [i for i in individual if not i["identity"]]
    tvs = sorted(i["tv"] for i in primary)
    result = {"scope": "Exact saved S4--S7 positive-genus quotient components; primary identity-free",
              "source_sha256": hashes, "components": len(primary),
              "states": sum(i["states"] for i in primary),
              "component_tv_median": (tvs[(len(tvs)-1)//2]+tvs[len(tvs)//2])/2,
              "caveat": "S7 source omits identity entries; inclusion sensitivity is complete only for S4--S6. No iid p-values or all-A-group claim.",
              "components_tv_le_005": sum(t <= .05 for t in tvs),
              "components_tv_le_010": sum(t <= .10 for t in tvs),
              "summaries": summaries,
              "best_components": sorted(primary, key=lambda i:i["tv"])[:10],
              "worst_components": sorted(primary, key=lambda i:i["tv"], reverse=True)[:5]}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
