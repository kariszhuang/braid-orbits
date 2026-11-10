#!/usr/bin/env python3
"""Second batch of primitive groups for the genus >= 2 attack (writes groups2.txt)."""
import itertools, sys
sys.path.insert(0, '.')
import make_groups as mg  # reuses GF and fmt (importing writes groups.txt again; harmless)
GF, fmt = mg.GF, mg.fmt
mg.IRRED[(7, 2)] = [3, 1, 1]
out = []
# projective lines, larger q
for (p, k) in [(37, 1), (41, 1), (43, 1), (47, 1), (7, 2), (53, 1), (59, 1), (61, 1)]:
    F = GF(p, k); q = F.q; INF = q
    def mob(a, b, c, dd):
        img = []
        for x in range(q + 1):
            if x == INF: num, den = a, c
            else: num, den = F.add[F.mul[a][x]][b], F.add[F.mul[c][x]][dd]
            img.append(INF if den == 0 else F.mul[num][F.inv[den]])
        return img
    t = mob(1, 1, 0, 1); w = F.w
    out.append(fmt(f"PSL(2,{q})", q + 1, [t, mob(F.mul[w][w], 0, 0, 1), mob(0, F.neg[1], 1, 0)]))
    out.append(fmt(f"PGL(2,{q})", q + 1, [t, mob(w, 0, 0, 1), mob(0, 1, 1, 0)]))
    if k > 1:
        fr = [INF if x == INF else F.frob[x] for x in range(q + 1)]
        out.append(fmt(f"PGammaL(2,{q})", q + 1, [t, mob(w, 0, 0, 1), mob(0, 1, 1, 0), fr]))
# affine lines, larger primes (tiny groups; all subgroups covered by AGL(1,p))
for p in [41, 43, 47, 53, 59, 61]:
    w = next(r for r in range(2, p) if len({pow(r, e, p) for e in range(p - 1)}) == p - 1)
    out.append(fmt(f"AGL(1,{p})", p, [[(x + 1) % p for x in range(p)], [(w * x) % p for x in range(p)]]))
# projective spaces PSL(n,q) on points
def psl(n, F, extra_field=False, frob=False):
    q = F.q; pts = []
    for v in itertools.product(range(q), repeat=n):
        if any(v) and next(c for c in v if c) == 1: pts.append(v)
    ix = {v: i for i, v in enumerate(pts)}
    def norm(v): li = F.inv[next(c for c in v if c)]; return tuple(F.mul[li][c] for c in v)
    def act(M):
        img = []
        for v in pts:
            w = []
            for r in range(n):
                s = 0
                for c in range(n): s = F.add[s][F.mul[M[r][c]][v[c]]]
                w.append(s)
            img.append(ix[norm(tuple(w))])
        return img
    gens = []
    for i in range(n):
        for j in range(n):
            if i != j:
                M = [[1 if r == c else 0 for c in range(n)] for r in range(n)]; M[i][j] = 1; gens.append(act(M))
                if extra_field: M[i][j] = F.w; gens.append(act(M))
    diag = None
    if F.q > 2:
        M = [[1 if r == c else 0 for c in range(n)] for r in range(n)]; M[0][0] = F.w; diag = act(M)
    if frob: gens.append([ix[norm(tuple(F.frob[c] for c in v))] for v in pts])
    return pts, gens, diag
F2 = GF(2); pts, g, _ = psl(4, F2); out.append(fmt("PSL(4,2)", len(pts), g))
F5 = GF(5); pts, g, _ = psl(3, F5); out.append(fmt("PSL(3,5)", len(pts), g))
F4 = GF(2, 2); pts, g, dg = psl(3, F4, extra_field=True); out.append(fmt("PGL(3,4)", len(pts), g + [dg]))
pts, g, dg = psl(3, F4, extra_field=True, frob=True); out.append(fmt("PGammaL(3,4)", len(pts), g + [dg]))
# S9, A9 on 2-subsets
for k in [9]:
    subs = list(itertools.combinations(range(k), 2)); ix = {s: i for i, s in enumerate(subs)}
    act = lambda perm: [ix[tuple(sorted(perm[x] for x in s))] for s in subs]
    cyc = [(i + 1) % k for i in range(k)]; tr = [1, 0] + list(range(2, k)); c3 = [1, 2, 0] + list(range(3, k))
    out.append(fmt(f"S{k}_on_2sets", len(subs), [act(cyc), act(tr)]))
    out.append(fmt(f"A{k}_on_2sets", len(subs), [act(cyc), act(c3)]))
# product actions H wr S2 on k^2 for primitive H
def cyc_perm(n, cycles):
    a = list(range(n))
    for c in cycles:
        for i in range(len(c)): a[c[i] - 1] = c[(i + 1) % len(c)] - 1
    return a
F5 = GF(5); INF5 = 5
def mob5(a, b, c, dd, F=F5, q=5):
    img = []
    for x in range(q + 1):
        if x == q: num, den = a, c
        else: num, den = F.add[F.mul[a][x]][b], F.add[F.mul[c][x]][dd]
        img.append(q if den == 0 else F.mul[num][F.inv[den]])
    return img
F7 = GF(7)
def mob7(a, b, c, dd): return mob5(a, b, c, dd, F7, 7)
bases = {
    "S4_on_4": (4, [[1, 2, 3, 0], [1, 0, 2, 3]]),
    "S5_on_6": (6, [mob5(1, 1, 0, 1), mob5(2, 0, 0, 1), mob5(0, 1, 1, 0)]),
    "A5_on_6": (6, [mob5(1, 1, 0, 1), mob5(4, 0, 0, 1), mob5(0, 4, 1, 0)]),
    "L2(7)_on_8": (8, [mob7(1, 1, 0, 1), mob7(2, 0, 0, 1), mob7(0, 6, 1, 0)]),
    "L2(7)_on_7": (7, [cyc_perm(7, [[1, 2, 3, 4, 5, 6, 7]]), cyc_perm(7, [[2, 3], [4, 7]])]),
}
for name, (k, hg) in bases.items():
    pts = list(itertools.product(range(k), repeat=2)); ix = {v: i for i, v in enumerate(pts)}
    L = lambda s: [ix[(s[a], b)] for a, b in pts]; sw = [ix[(b, a)] for a, b in pts]
    out.append(fmt(f"{name}_wrS2_product", k * k, [L(h) for h in hg] + [sw]))
# diagonal type: A5 x A5 acting on A5 by x -> a x b^-1, with inversion and S5 conjugation
A5 = [p for p in itertools.permutations(range(5)) if sum(1 for i in range(5) for j in range(i) if p[j] > p[i]) % 2 == 0]
ixA = {p: i for i, p in enumerate(A5)}
mul = lambda a, b: tuple(a[i] for i in b); inv = lambda a: tuple(sorted(range(5), key=lambda i: a[i]))
a1 = (1, 2, 3, 4, 0); a2 = (1, 2, 0, 3, 4); tr = (1, 0, 2, 3, 4)
left = lambda a: [ixA[mul(a, x)] for x in A5]; right = lambda b: [ixA[mul(x, inv(b))] for x in A5]
invm = [ixA[inv(x)] for x in A5]; conj = [ixA[mul(tr, mul(x, tr))] for x in A5]
out.append(fmt("A5xA5_diag", 60, [left(a1), left(a2), right(a1), right(a2)]))
out.append(fmt("(A5xA5).2_diag_inv", 60, [left(a1), left(a2), right(a1), right(a2), invm]))
out.append(fmt("(A5xA5).2^2_diag", 60, [left(a1), left(a2), right(a1), right(a2), invm, conj]))
# M22 (standard generators; the tool reports the actual order)
out.append(fmt("M22?", 22, [cyc_perm(22, [list(range(1, 12)), list(range(12, 23))]),
    cyc_perm(22, [[1, 4, 5, 9, 3], [2, 8, 10, 7, 6], [12, 15, 16, 20, 14], [13, 19, 21, 18, 17]]),
    cyc_perm(22, [[1, 21], [2, 10, 8, 6], [3, 13, 4, 17], [5, 19, 9, 18], [11, 22], [12, 14, 16, 20]])]))
open("groups2.txt", "w").write("\n".join(out) + "\n"); print(len(out), "groups")
