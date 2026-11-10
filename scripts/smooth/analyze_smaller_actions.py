#!/usr/bin/env python3
"""Read-only exact group/position inspection of C141's two additional examples."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
import sys
from collections import Counter

from run_smooth_mission import rows, sha
from run_smooth_smaller import analyze
from verify_smooth_geometry import compose, cycles, parse_tuple, require
from sympy.combinatorics import Permutation, PermutationGroup


def inspect(root):
    catalog = json.loads((root/'CATALOG.json').read_text())
    out = []
    v4 = [(0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)]
    rotation = (1, 2, 0, 3)
    cosets, power = {}, tuple(range(4))
    for exponent in range(3):
        for v in v4:
            cosets[compose(v, power)] = exponent
        power = compose(rotation, power)
    require(len(cosets) == 12, 'A4/V4 coset table')
    for item in catalog['closed_components']:
        if item['known_component']:
            continue
        folder = root/item['folder']
        require(sha(folder/'STATES.tsv') == item['table_sha256'], 'changed state table')
        t = parse_tuple(item['summary']['seed'], 12)
        analyze(t)
        group = PermutationGroup([Permutation(list(p)) for p in t])
        # Given seed coordinates, quotient by independently acting V4 on each
        # of three four-point blocks. Check its full containment, not order alone.
        for block in range(3):
            for local in v4[1:3]:
                p = list(range(12))
                p[4*block:4*block+4] = [4*block+x for x in local]
                require(group.contains(Permutation(p)), 'V4^3 not contained')
        quotient_generators, vectors = [], []
        for p in t:
            image, vector = [], []
            for block in range(3):
                target = p[4*block]//4
                require(all(p[4*block+u]//4 == target for u in range(4)), 'block invariant')
                local = tuple(p[4*block+u] % 4 for u in range(4))
                require(local in cosets, 'not an A4 local action')
                e = cosets[local]
                vector.append(e)
                image.extend(3*target+(z+e) % 3 for z in range(3))
            require(sum(vector) % 3 == 0, 'cubic sum constraint failed')
            vectors.append(vector)
            quotient_generators.append(Permutation(image))
        quotient = PermutationGroup(quotient_generators)
        quotient_order = int(quotient.order())
        require(int(group.order()) == 64*quotient_order, 'kernel bigger than V4^3')
        element_orders = Counter(int(p.order()) for p in quotient.generate_schreier_sims())
        fibers = Counter()
        for row in rows(folder/'STATES.tsv'):
            u = parse_tuple(row['tuple'], 12)
            small = analyze(u)['small_entry']
            fibers[small] += 1
        require(len(fibers) == 4 and len(set(fibers.values())) == 1, 'unequal position fibers')
        out.append(dict(folder=item['folder'], group_order=int(group.order()),
                        normal_subgroup='V4^3', normal_subgroup_order=64,
                        quotient_order=quotient_order, quotient_abelian=bool(quotient.is_abelian),
                        quotient_element_orders=dict(element_orders),
                        local_cubic_vectors=vectors, position_state_counts=dict(fibers),
                        geometry='isotrivial_elliptic; human review pending'))
    require(len(out) == 2, 'expected exactly two additional C141 components')
    return out


if __name__ == '__main__':
    print(json.dumps(inspect(Path(sys.argv[1])), indent=2, sort_keys=True))
