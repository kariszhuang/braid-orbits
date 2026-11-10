#!/usr/bin/env python3
"""Checks for the descent theorem (DESCENT.md). Independent of the Rust code.

For tuples (g1,g2,g3,g4), product 1 (rightmost-first), preserving the standard system of
n blocks of size m (i.e. in S_m wr S_n), transitive, it checks at the cusp {12|34}:
  T0  every <x,y>-orbit O maps onto its image Obar with fibres of one size k = |O|/|Obar|
  T1  c_O(a) <= k * c_Obar(abar) for a = x, y, p      (cycle-count inequality)
  T2  2 g_O - 2 >= k (2 gbar - 2) for every piece       (piece Riemann-Hurwitz inequality)
  T3  #positive-genus pieces: quotient <= source
  T4  defect (b1 of the cusp graph): quotient <= source
  T5  the induced map H1(Gamma; Q) -> H1(Gammabar; Q) is onto (rank computed exactly)
  T6  one-piece(source) => one-piece(quotient)
Exhaustive over S2 wr S3 and S3 wr S2 (all g1,g2,g3), random for larger wreath products.
usage: descent_check.py [random_trials]
"""
import itertools, random, sys
from fractions import Fraction

def mul(a, b): return tuple(a[i] for i in b)          # (a*b)[i] = a[b[i]]
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
def orbits(gens, pts):
    pts = set(pts); out = []
    while pts:
        s = pts.pop(); O = {s}; st = [s]
        while st:
            x = st.pop()
            for g in gens:
                y = g[x]
                if y not in O: O.add(y); st.append(y); pts.discard(y)
        out.append(frozenset(O))
    return out

def cusp(t, n):
    """Cusp graph at {12|34}: vertices = orbits of <g1,g2> (side A) and <g3,g4> (side B),
    edges = cycles of p = g1 g2. Returns (A, B, edges, genus per vertex)."""
    x, y, z, w = t; p = mul(x, y)
    A = orbits([x, y], range(n)); B = orbits([z, w], range(n))
    edges = [frozenset(c) for c in cycles(p)]
    def gen(O, a, b):
        c = sum(1 for C in cycles(a) if C[0] in O) + sum(1 for C in cycles(b) if C[0] in O) + sum(1 for C in edges if next(iter(C)) in O)
        tg = len(O) + 2 - c; assert tg >= 0 and tg % 2 == 0; return tg // 2
    gA = [gen(O, x, y) for O in A]; gB = [gen(O, z, w) for O in B]
    return A, B, edges, gA, gB, p

def rank_q(rows, ncols):
    """rank over Q by Gaussian elimination"""
    M = [[Fraction(v) for v in r] for r in rows]; rk = 0; col = 0
    while rk < len(M) and col < ncols:
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None: col += 1; continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / M[rk][col]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1; col += 1
    return rk

def cycle_space(V, E):
    """basis of H1 of a graph: V vertex list, E list of (u, v); returns list of edge vectors"""
    # incidence kernel via spanning tree: one fundamental cycle per non-tree edge
    idx = {v: i for i, v in enumerate(V)}; par = list(range(len(V)))
    def f(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    tree = []; non = []
    for e, (u, v) in enumerate(E):
        a, b = f(idx[u]), f(idx[v])
        if a != b: par[a] = b; tree.append(e)
        else: non.append(e)
    adj = {v: [] for v in V}
    for e in tree: u, v = E[e]; adj[u].append((v, e, 1)); adj[v].append((u, e, -1))
    def path(s, t):  # signed edge vector of the tree path s -> t
        prev = {s: None}; st = [s]
        while st:
            a = st.pop()
            for b, e, sg in adj[a]:
                if b not in prev: prev[b] = (a, e, sg); st.append(b)
        vec = [0] * len(E); c = t
        while prev[c] is not None: a, e, sg = prev[c]; vec[e] += sg; c = a
        return vec
    basis = []
    for e in non:
        u, v = E[e]; vec = path(v, u)  # v -> u in the tree, then edge u -> v closes it
        vec[e] += 1; basis.append(vec)
    return basis

def check(t, m, nb, stats):
    n = m * nb
    blk = lambda i: i // m
    tb = tuple(tuple(blk(g[i * m]) for i in range(nb)) for g in t)
    A, B, E, gA, gB, p = cusp(t, n)
    Ab, Bb, Eb, gAb, gBb, pb = cusp(tb, nb)
    # T0 T1 T2 (side A with x,y; side B with z,w)
    for side, Q, gs, gq, (a1, a2) in ((A, Ab, gA, gAb, (0, 1)), (B, Bb, gB, gBb, (2, 3))):
        for O, gO in zip(side, gs):
            img = frozenset(blk(i) for i in O)
            j = Q.index(img)  # T0: image is a quotient orbit
            sizes = {sum(1 for i in O if blk(i) == b) for b in img}; assert len(sizes) == 1; k = sizes.pop()
            assert k * len(img) == len(O)
            for a, ab in ((t[a1], tb[a1]), (t[a2], tb[a2]), (p if side is A else inv(p), pb if side is A else inv(pb))):
                co = sum(1 for C in cycles(a) if C[0] in O); cq = sum(1 for C in cycles(ab) if C[0] in img)
                assert co <= k * cq, "T1"
            assert 2 * gO - 2 >= k * (2 * gq[j] - 2), "T2"
    pos = sum(g > 0 for g in gA + gB); posb = sum(g > 0 for g in gAb + gBb)
    assert posb <= pos, "T3"
    d = len(E) - len(A) - len(B) + 1; db = len(Eb) - len(Ab) - len(Bb) + 1
    assert d >= 0 and db >= 0 and db <= d, "T4"
    # T5: rank of H1(Gamma) -> H1(Gammabar)
    V = [('A', O) for O in A] + [('B', O) for O in B]
    def ends(C, AA, BB): s = next(iter(C)); return (('A', next(O for O in AA if s in O)), ('B', next(O for O in BB if s in O)))
    El = [ends(C, A, B) for C in E]
    Vb = [('A', O) for O in Ab] + [('B', O) for O in Bb]
    Ebl = [ends(C, Ab, Bb) for C in Eb]
    emap = [Eb.index(next(C2 for C2 in Eb if blk(next(iter(C))) in C2)) for C in E]  # edge image (orientation A->B preserved)
    basis = cycle_space(V, El)
    imgs = []
    for vec in basis:
        w = [0] * len(Eb)
        for e, c in enumerate(vec): w[emap[e]] += c
        imgs.append(w)
    if db > 0:
        bq = cycle_space(Vb, Ebl)
        r_all = rank_q(bq, len(Eb)); assert r_all == db
        r_img = rank_q(imgs, len(Eb)) if imgs else 0
        # the images lie in H1(Gammabar) automatically (graph map); onto <=> rank = db
        assert rank_q(imgs + bq, len(Eb)) == db, "images not cycles?"
        assert r_img == db, "T5"
    one = d == 0 and pos <= 1; oneb = db == 0 and posb <= 1
    if one: assert oneb, "T6"
    stats['tuples'] += 1; stats['one'] += one; stats['defect_pos'] += d > 0; stats['quot_defect_pos'] += db > 0
    stats['strict_defect_drop'] += db < d; stats['quot_pos_genus'] += posb > 0

def wreath(m, nb):
    out = []
    for top in itertools.permutations(range(nb)):
        for fib in itertools.product(list(itertools.permutations(range(m))), repeat=nb):
            out.append(tuple(top[i // m] * m + fib[i // m][i % m] for i in range(m * nb)))
    return out

def transitive(t, n): return len(orbits(t, range(n))) == 1

def run_exhaustive(m, nb):
    W = wreath(m, nb); n = m * nb; ident = tuple(range(n))
    stats = dict(tuples=0, one=0, defect_pos=0, quot_defect_pos=0, strict_defect_drop=0, quot_pos_genus=0)
    for g1 in W:
        for g2 in W:
            q = mul(g1, g2)
            for g3 in W:
                g4 = inv(mul(q, g3)); t = (g1, g2, g3, g4)
                if not transitive(t, n): continue
                check(t, m, nb, stats)
    return stats

def run_random(trials, seed=1):
    rnd = random.Random(seed); stats = dict(tuples=0, one=0, defect_pos=0, quot_defect_pos=0, strict_defect_drop=0, quot_pos_genus=0)
    shapes = [(2, 4), (2, 5), (3, 3), (4, 2), (3, 4), (2, 6), (5, 2), (4, 3)]
    def relt(m, nb):
        top = list(range(nb)); rnd.shuffle(top); fib = [rnd.sample(range(m), m) for _ in range(nb)]
        return tuple(top[i // m] * m + fib[i // m][i % m] for i in range(m * nb))
    while stats['tuples'] < trials:
        m, nb = rnd.choice(shapes); n = m * nb
        g = [relt(m, nb) for _ in range(3)]; g.append(inv(mul(g[0], mul(g[1], g[2]))))
        if not transitive(g, n): continue
        check(tuple(g), m, nb, stats)
    return stats

if __name__ == "__main__":
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    for (m, nb) in [(2, 3), (3, 2)]:
        print(f"exhaustive S{m} wr S{nb}:", run_exhaustive(m, nb), flush=True)
    print(f"random ({trials}):", run_random(trials), flush=True)
    print("ALL CHECKS PASSED (T0-T6)")
