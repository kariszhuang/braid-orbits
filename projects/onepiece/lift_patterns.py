#!/usr/bin/env python3
"""Where do the twists sit in genus-3 lifts that are one-piece everywhere?

For each recorded orbit (a tuple of degree m*b), find a block system of block size m whose
quotient has genus 1, read off for every quotient cycle c of entry i the return permutation
(the action of g_i^len(c) on one block of c), and record the excess vector
e = (e1, e2, e3, e4), e_i = sum over cycles of entry i of ind(return), plus the cycle types
of the nontrivial returns. sum e_i = 4 for genus 3 over genus 1.
usage: lift_patterns.py M FILE...
"""
import json, sys, collections

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

def block_systems(t, n, m):
    """all block systems with blocks of size m (minimal block containing 0 and j)"""
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
        if all(len(B) == m for B in blocks.values()):
            res.add(tuple(sorted(tuple(sorted(B)) for B in blocks.values())))
    return res

def genus(t, n):
    s = sum(n - len(cycles(g)) for g in t); return (s - 2 * n + 2) // 2

def analyse(t, n, m):
    for bs in block_systems(t, n, m):
        idx = {x: k for k, B in enumerate(bs) for x in B}; nb = len(bs)
        q = [[idx[g[B[0]]] for B in bs] for g in t]
        if genus(q, nb) != 1: continue
        e = [0, 0, 0, 0]; types = []
        for i, g in enumerate(t):
            for c in cycles(q[i]):
                L = len(c); B = bs[c[0]]; pos = {x: k for k, x in enumerate(B)}
                r = []
                for x in B:
                    y = x
                    for _ in range(L): y = g[y]
                    r.append(pos[y])
                ct = sorted((len(z) for z in cycles(r)), reverse=True)
                ind = m - len(ct)
                if ind: e[i] += ind; types.append((i + 1, L, tuple(z for z in ct if z > 1)))
        return tuple(e), tuple(sorted(types))
    return None, None

def main():
    m = int(sys.argv[1]); files = sys.argv[2:]
    support = collections.Counter(); evec = collections.Counter(); kinds = collections.Counter(); n_or = 0; fails = 0
    for f in files:
        for line in open(f):
            if '"states"' not in line: continue
            r = json.loads(line); n, t = parse(r['tuple'])
            e, types = analyse(t, n, m)
            if e is None: fails += 1; continue
            n_or += 1
            assert sum(e) == 4, (e, r['tuple'])
            support[sum(1 for x in e if x)] += 1
            evec[tuple(sorted(e, reverse=True))] += 1
            kinds[tuple(sorted(tt[2] for tt in types))] += 1
    print(f"m={m}: {n_or} orbits analysed ({fails} without a genus-1 quotient of block size {m})")
    print("  number of entries carrying twist (at the recorded state):", dict(sorted(support.items())))
    print("  excess vectors (sorted):", dict(evec.most_common()))
    print("  twist cycle types:", dict(kinds.most_common(12)))

if __name__ == "__main__":
    main()
