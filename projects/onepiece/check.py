#!/usr/bin/env python3
"""Independent check of onepiece (shares no code with src/main.rs).

Brute force over ALL tuples (g1,g2,g3) in S_d^3 with g4 = (g1 g2 g3)^-1, composition
rightmost-first, identity-free, transitive, of the given genus. Orbits are found by
union-find over the six moves on labelled tuples (no canonical forms), then taken
modulo simultaneous conjugation by counting labelled tuples. For every orbit it
reports clusters (<sigma1,sigma3>-orbits of conjugacy classes of states), bad
clusters, and whether the orbit is one-piece everywhere.

usage: check.py D GENUS [all|onepiece]  -> prints a histogram to compare with
the census output out/all_g2_dD.jsonl or out/op_dD_gG.jsonl.
"""
import itertools, sys, collections

def mul(a, b):  # (a*b)[i] = a[b[i]]
    return tuple(a[x] for x in b)

def inv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)

def cycles(a):
    seen, out = set(), []
    for s in range(len(a)):
        if s in seen: continue
        c, x = [], s
        while x not in seen: seen.add(x); c.append(x); x = a[x]
        out.append(c)
    return out

def orbit_partition(gens, n):
    lab = list(range(n))
    def f(x):
        while lab[x] != x: lab[x] = lab[lab[x]]; x = lab[x]
        return x
    for g in gens:
        for i in range(n): lab[f(i)] = f(g[i])
    groups = collections.defaultdict(list)
    for i in range(n): groups[f(i)].append(i)
    return list(groups.values())

def one_piece(t, n):
    p = mul(t[0], t[1])
    pc = cycles(p)
    sides = [orbit_partition([t[0], t[1]], n), orbit_partition([t[2], t[3]], n)]
    if len(pc) != len(sides[0]) + len(sides[1]) - 1: return False
    pos = 0
    for side, (x, y) in zip(sides, [(t[0], t[1]), (t[2], t[3])]):
        for O in side:
            S = set(O)
            c = sum(1 for C in cycles(x) if C[0] in S) + sum(1 for C in cycles(y) if C[0] in S) + sum(1 for C in pc if C[0] in S)
            if len(O) + 2 - c > 0: pos += 1
    return pos <= 1

def move(t, m):
    j = m % 3; a, b = t[j], t[j + 1]; u = list(t)
    if m < 3: u[j] = mul(a, mul(b, inv(a))); u[j + 1] = a
    else: u[j] = b; u[j + 1] = mul(inv(b), mul(a, b))
    return tuple(u)

def main():
    n, genus = int(sys.argv[1]), int(sys.argv[2])
    perms = list(itertools.permutations(range(n)))
    ind = {p: n - len(cycles(p)) for p in perms}
    target = 2 * n - 2 + 2 * genus
    idt = tuple(range(n))
    tuples = []
    for g1 in perms:
        if g1 == idt: continue
        for g2 in perms:
            if g2 == idt or ind[g1] + ind[g2] > target - 2: continue
            p = mul(g1, g2)
            for g3 in perms:
                if g3 == idt: continue
                q = mul(p, g3)
                if q == idt or ind[g1] + ind[g2] + ind[g3] + ind[q] != target: continue
                t = (g1, g2, g3, inv(q))
                if len(orbit_partition(t[:3], n)) == 1: tuples.append(t)
    index = {t: i for i, t in enumerate(tuples)}
    par = list(range(len(tuples)))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    gensS0 = [tuple([1, 0] + list(range(2, n))), tuple(list(range(1, n)) + [0])]
    for i, t in enumerate(tuples):
        for m in range(3): par[f(i)] = f(index[move(t, m)])
        for c in gensS0:  # orbits are taken modulo simultaneous conjugation
            ci = inv(c); par[f(i)] = f(index[tuple(mul(c, mul(x, ci)) for x in t)])
    orbits = collections.defaultdict(list)
    for i in range(len(tuples)): orbits[f(i)].append(i)
    # conjugacy classes of states = orbits of S_n on labelled tuples; clusters use moves 0 and 2
    hist = collections.Counter(); all_one = 0
    for root, members in orbits.items():
        mem = set(members)
        # union over sigma1, sigma3 and conjugation by all of S_n
        cp = {i: i for i in members}
        def g(x):
            while cp[x] != x: cp[x] = cp[cp[x]]; x = cp[x]
            return x
        for i in members:
            t = tuples[i]
            for m in (0, 2): cp[g(i)] = g(index[move(t, m)])
        # conjugation: it suffices to use generators (0 1) and (0 1 ... n-1)
        gensS = [tuple([1, 0] + list(range(2, n))), tuple(list(range(1, n)) + [0])]
        for i in members:
            t = tuples[i]
            for c in gensS:
                ci = inv(c); u = tuple(mul(c, mul(x, ci)) for x in t)
                if u in index: cp[g(i)] = g(index[u])
        clusters = collections.defaultdict(list)
        for i in members: clusters[g(i)].append(i)
        bad = 0
        for c, ms in clusters.items():
            vals = {one_piece(tuples[i], n) for i in ms}
            assert len(vals) == 1, "one-piece not constant on a cluster"
            if not vals.pop(): bad += 1
        hist[(len(clusters), bad)] += 1
        if bad == 0: all_one += 1
    print(f"d={n} genus={genus} labelled_tuples={len(tuples)} orbits={len(orbits)} all_one_piece_orbits={all_one}")
    print("histogram (clusters, bad_clusters): count", sorted(hist.items()))

if __name__ == "__main__":
    main()
