#!/usr/bin/env python3
"""Classify every all-zero orbit of genus >= 2 found so far by its STRUCTURE, not just primitivity.

For each orbit (one representative tuple): generated group G (order, abelian, cyclic, solvable,
composition-factor orders), all block systems, and for each block system the quotient tuple's
genus, whether the quotient is all-zero... (by the descent theorem it always is) and its group.
A record is 'C_p-type' if G is cyclic, or if G is solvable with every composition factor cyclic
of prime order AND some block quotient of positive genus is cyclic (a tower over a
Deligne-Mostow cyclic piece). Everything else is flagged NONCYCLIC-TYPE.
usage: classify_allzero.py FILE.jsonl ...
"""
import json, sys, collections
from sympy.combinatorics import Permutation, PermutationGroup

def parse(tup):
    n = max(int(x) for p in tup for x in p.replace('(', ' ').replace(')', ' ').split())
    out = []
    for p in tup:
        a = list(range(n))
        for cyc in p.split(')'):
            v = [int(x) - 1 for x in cyc.replace('(', ' ').split()]
            for i in range(len(v)): a[v[i]] = v[(i + 1) % len(v)]
        out.append(a)
    return n, out

def cycles(a):
    seen, out = set(), []
    for s in range(len(a)):
        if s in seen: continue
        c, x = [], s
        while x not in seen: seen.add(x); c.append(x); x = a[x]
        out.append(c)
    return out
def genus(t, n): return (sum(n - len(cycles(g)) for g in t) - 2 * n + 2) // 2

def block_systems(t, n):
    res = set()
    for j in range(1, n):
        par = list(range(n))
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        q = [(0, j)]
        while q:
            x, y = q.pop(); a, b = f(x), f(y)
            if a == b: continue
            par[a] = b
            for g in t: q.append((g[x], g[y]))
        blocks = collections.defaultdict(list)
        for x in range(n): blocks[f(x)].append(x)
        if 1 < len(blocks) < n:
            res.add(tuple(sorted(tuple(sorted(B)) for B in blocks.values())))
    return res

def group_info(t):
    G = PermutationGroup([Permutation(g) for g in t])
    if G.is_solvable:
        cs = G.composition_series(); cf = sorted(int(H.order() // K.order()) for H, K in zip(cs, cs[1:]))
    else:
        cf = ['nonsolvable']
    return dict(order=int(G.order()), abelian=G.is_abelian, cyclic=G.is_cyclic, solvable=G.is_solvable, comp_factors=cf)

def classify(tup):
    n, t = parse(tup)
    gi = group_info(t)
    quots = []
    for bs in block_systems(t, n):
        idx = {x: k for k, B in enumerate(bs) for x in B}
        q = [[idx[g[B[0]]] for B in bs] for g in t]
        if any(all(q[i][k] == k for k in range(len(bs))) for i in range(4)):
            qg = None  # identity entry in the quotient
        qi = group_info(q) if len(bs) > 1 else None
        quots.append(dict(nblocks=len(bs), block=len(bs[0]), genus=genus(q, len(bs)), cyclic=qi['cyclic'] if qi else True, order=qi['order'] if qi else 1))
    cp_tower = gi['solvable'] and all(isinstance(c, int) and _is_prime(c) for c in gi['comp_factors']) and any(qq['genus'] > 0 and qq['cyclic'] for qq in quots)
    kind = "cyclic" if gi['cyclic'] else ("NONSOLVABLE" if not gi['solvable'] else ("C_p-tower" if cp_tower else "solvable-Cp-factors"))
    return n, genus(t, n), gi, quots, kind

def _is_prime(k): return k > 1 and all(k % p for p in range(2, int(k ** 0.5) + 1))

if __name__ == "__main__":
    seen = set(); summary = collections.Counter()
    for f in sys.argv[1:]:
        for line in open(f):
            r = json.loads(line)
            if 'tuple' not in r or r.get('capped'): continue
            if r.get('genus', 0) < 2: continue
            if 'all_zero' in r and not r['all_zero']: continue      # census files: keep all-zero only
            key = (f.split('/')[-1], tuple(r['tuple']))
            if key in seen: continue
            seen.add(key)
            n, g, gi, quots, kind = classify(r['tuple'])
            summary[kind] += 1
            if kind != "cyclic":
                print(f"{kind:15s} d={n} g={g} src={r.get('group', f.split('/')[-1])} |G|={gi['order']} solvable={gi['solvable']} factors={gi['comp_factors'][:8]} "
                      f"quotients={[(qq['nblocks'], qq['genus'], 'cyc' if qq['cyclic'] else qq['order']) for qq in quots]}", flush=True)
    print("SUMMARY (genus >= 2 all-zero records):", dict(summary))
