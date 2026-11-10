#!/usr/bin/env python3
"""Distinct genus>=2 all-zero orbits with NON-SOLVABLE group: independent closure (verify_cex),
group data, and the two block quotients (do they have identity entries?)."""
import json, sys, collections
sys.path.insert(0, '.')
from verify_cex import parse, verify, canon, move
from classify_allzero import block_systems, genus, cycles
from sympy.combinatorics import Permutation, PermutationGroup
recs = []
for f in sys.argv[1:]:
    for l in open(f):
        r = json.loads(l)
        if 'tuple' in r and r.get('genus', 0) >= 2: recs.append(r)
out = {}
for r in recs:
    n, t = parse(r['tuple'])
    G = PermutationGroup([Permutation(list(g)) for g in t])
    if G.is_solvable: continue
    # orbit id: least canonical state over the closed orbit
    seen = {canon(tuple(t), n): tuple(t)}; q = [tuple(t)]
    while q:
        u = q.pop()
        for m in range(6):
            v = move(u, m); k = canon(v, n)
            if k not in seen: seen[k] = v; q.append(v)
    key = (n, min(seen))
    if key in out: continue
    v = verify(r['tuple'])
    quots = []
    for bs in block_systems(t, n):
        idx = {x: k for k, B in enumerate(bs) for x in B}
        qt = [[idx[g[B[0]]] for B in bs] for g in t]
        ids = sum(1 for g in qt if all(g[k] == k for k in range(len(bs))))
        quots.append(dict(blocks=len(bs), genus=genus(qt, len(bs)), identity_entries=ids, order=int(PermutationGroup([Permutation(g) for g in qt]).order()),
                          types=[sorted((len(c) for c in cycles(g) if len(c) > 1), reverse=True) for g in qt]))
    out[key] = dict(d=n, genus=v['genus'], states=v['states'], defects=v['defect_histogram'], order=v['order'], primitive=v['primitive'],
                    types=[sorted((len(c) for c in cycles(g) if len(c) > 1), reverse=True) for g in t], quotients=quots, tuple=r['tuple'])
for o in out.values(): print(json.dumps(o))
print(len(out), "distinct non-solvable genus>=2 all-zero orbits", file=sys.stderr)
