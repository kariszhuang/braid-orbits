#!/usr/bin/env python3
"""Generators of primitive groups (not containing A_d) for the all-zero disproof attempt.
Writes groups.txt: one line per group: NAME<TAB>DEGREE<TAB>gen;gen;... with each gen given as
space-separated 0-based images. Orders and primitivity are re-checked by the Rust tool.
"""
import itertools, sys

# ---------------- finite fields GF(p^k) as integers 0..q-1 (polynomial basis)
IRRED = {(2, 2): [1, 1, 1], (2, 3): [1, 1, 0, 1], (2, 4): [1, 1, 0, 0, 1], (2, 5): [1, 0, 1, 0, 0, 1],
         (3, 2): [2, 2, 1], (3, 3): [1, 2, 0, 1], (5, 2): [2, 1, 1], (7, 2): [3, 1, 1]}
class GF:
    def __init__(self, p, k=1):
        self.p, self.k, self.q = p, k, p ** k
        def vec(x): return [(x // p ** i) % p for i in range(k)]
        def num(v): return sum(c * p ** i for i, c in enumerate(v))
        self.add = [[num([(a + b) % p for a, b in zip(vec(x), vec(y))]) for y in range(self.q)] for x in range(self.q)]
        if k == 1:
            self.mul = [[x * y % p for y in range(p)] for x in range(p)]
        else:
            f = IRRED[(p, k)]  # monic poly coefficients low..high, degree k
            def pm(a, b):
                r = [0] * (2 * k - 1)
                for i, x in enumerate(a):
                    for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
                for d in range(2 * k - 2, k - 1, -1):
                    c = r[d]
                    if c:
                        for i in range(k + 1): r[d - k + i] = (r[d - k + i] - c * f[i]) % p
                return r[:k]
            self.mul = [[num(pm(vec(x), vec(y))) for y in range(self.q)] for x in range(self.q)]
        self.neg = [next(y for y in range(self.q) if self.add[x][y] == 0) for x in range(self.q)]
        self.inv = [None] + [next(y for y in range(1, self.q) if self.mul[x][y] == 1) for x in range(1, self.q)]
        # primitive element
        for w in range(2, self.q) if self.q > 2 else [1]:
            s, x = set(), 1
            for _ in range(self.q - 1): x = self.mul[x][w]; s.add(x)
            if len(s) == self.q - 1: self.w = w; break
        self.frob = [self.pw(x, p) for x in range(self.q)]
    def pw(self, x, e):
        r = 1
        for _ in range(e): r = self.mul[r][x]
        return r

def fmt(name, d, gens):
    for g in gens: assert sorted(g) == list(range(d)), (name, g)
    return f"{name}\t{d}\t" + ";".join(" ".join(map(str, g)) for g in gens)

out = []

# ---------------- controls
out.append(fmt("C7(control)", 7, [[(i + 1) % 7 for i in range(7)]]))
out.append(fmt("D5(control)", 5, [[(i + 1) % 5 for i in range(5)], [(-i) % 5 for i in range(5)]]))

# ---------------- projective line groups PSL/PGL/PGammaL(2,q) on q+1 points (infinity = q)
for (p, k) in [(5, 1), (7, 1), (2, 3), (3, 2), (11, 1), (13, 1), (2, 4), (17, 1), (19, 1), (23, 1), (5, 2), (3, 3), (29, 1), (31, 1), (2, 5)]:
    F = GF(p, k); q = F.q; INF = q
    def mob(a, b, c, dd):  # x -> (a x + b)/(c x + dd)
        img = []
        for x in range(q + 1):
            if x == INF: num, den = a, c
            else: num, den = F.add[F.mul[a][x]][b], F.add[F.mul[c][x]][dd]
            img.append(INF if den == 0 else F.mul[num][F.inv[den]])
        return img
    t = mob(1, 1, 0, 1); w = F.w
    pgl = [t, mob(w, 0, 0, 1), mob(0, 1, 1, 0)]
    psl = [t, mob(F.mul[w][w], 0, 0, 1), mob(0, F.neg[1], 1, 0)]
    out.append(fmt(f"PSL(2,{q})", q + 1, psl))
    if q % 2 == 1: out.append(fmt(f"PGL(2,{q})", q + 1, pgl))
    if k > 1:
        fr = [INF if x == INF else F.frob[x] for x in range(q + 1)]
        out.append(fmt(f"PGammaL(2,{q})", q + 1, pgl + [fr]))
        out.append(fmt(f"PSigmaL(2,{q})", q + 1, psl + [fr]))

# ---------------- affine groups of dimension 1: C_p^k x| C_m (m | q-1) and AGammaL(1,q)
for (p, k) in [(2, 3), (3, 2), (2, 4), (5, 2), (3, 3), (2, 5), (17, 1), (19, 1), (23, 1), (29, 1), (31, 1), (37, 1)]:
    F = GF(p, k); q = F.q
    tr = [F.add[x][1] for x in range(q)]
    for m in [d for d in range(2, q) if (q - 1) % d == 0]:
        a = F.pw(F.w, (q - 1) // m)
        out.append(fmt(f"AGL1sub(q={q},m={m})", q, [tr, [F.mul[a][x] for x in range(q)]]))
    if k > 1:
        out.append(fmt(f"AGammaL(1,{q})", q, [tr, [F.mul[F.w][x] for x in range(q)], F.frob]))

# ---------------- AGL(n,p) on p^n points, and ASL
def agl(n, p, special=False):
    pts = list(itertools.product(range(p), repeat=n)); ix = {v: i for i, v in enumerate(pts)}
    gens = []
    gens.append([ix[tuple(((v[0] + 1) % p,) + v[1:])] for v in pts])  # translation
    for i in range(n):
        for j in range(n):
            if i != j:  # transvection x_i += x_j
                gens.append([ix[tuple((v[a] + (v[j] if a == i else 0)) % p for a in range(n))] for v in pts])
    if not special and p > 2:
        w = next(r for r in range(2, p) if len({pow(r, e, p) for e in range(p - 1)}) == p - 1)
        gens.append([ix[((v[0] * w) % p,) + tuple(v[1:])] for v in pts])
    return pts, gens
for (n, p) in [(2, 3), (3, 2), (2, 5), (4, 2), (3, 3), (2, 7)]:
    pts, g = agl(n, p); out.append(fmt(f"AGL({n},{p})", p ** n, g))
    if p > 2: pts, g = agl(n, p, True); out.append(fmt(f"ASL({n},{p})", p ** n, g))

# ---------------- PSL(3,q) on the points of PG(2,q)
for (p, k) in [(2, 1), (3, 1), (2, 2)]:
    F = GF(p, k); q = F.q
    pts = []
    for v in itertools.product(range(q), repeat=3):
        if v == (0, 0, 0): continue
        lead = next(c for c in v if c)
        if lead == 1: pts.append(v)
    ix = {v: i for i, v in enumerate(pts)}
    def norm(v):
        lead = next(c for c in v if c); li = F.inv[lead]; return tuple(F.mul[li][c] for c in v)
    def act(M):
        img = []
        for v in pts:
            w = tuple(F.add[F.add[F.mul[M[r][0]][v[0]]][F.mul[M[r][1]][v[1]]]][F.mul[M[r][2]][v[2]]] for r in range(3))
            img.append(ix[norm(w)])
        return img
    gens = []
    for i in range(3):
        for j in range(3):
            if i != j:
                M = [[1 if r == c else 0 for c in range(3)] for r in range(3)]; M[i][j] = 1; gens.append(act(M))
    if q == 4:  # a field-generating transvection too
        M = [[1, F.w, 0], [0, 1, 0], [0, 0, 1]]; gens.append(act(M))
    out.append(fmt(f"PSL(3,{q})", len(pts), gens))

# ---------------- S_k and A_k on 2-subsets and 3-subsets
for k, t in [(5, 2), (6, 2), (7, 2), (8, 2), (7, 3), (6, 3)]:
    subs = list(itertools.combinations(range(k), t)); ix = {s: i for i, s in enumerate(subs)}
    def act(perm): return [ix[tuple(sorted(perm[x] for x in s))] for s in subs]
    cyc = [(i + 1) % k for i in range(k)]; tr = [1, 0] + list(range(2, k))
    out.append(fmt(f"S{k}_on_{t}sets", len(subs), [act(cyc), act(tr)]))
    c3 = [1, 2, 0] + list(range(3, k))
    cyc_even = cyc if k % 2 == 1 else [(i + 1) % (k - 1) for i in range(k - 1)] + [k - 1]
    out.append(fmt(f"A{k}_on_{t}sets", len(subs), [act(cyc_even), act(c3)]))

# ---------------- product action S5 wr S2 and (A5 x A5).2 on 25 points
k = 5; pts = list(itertools.product(range(k), repeat=2)); ix = {v: i for i, v in enumerate(pts)}
cyc = [(i + 1) % k for i in range(k)]; tr = [1, 0, 2, 3, 4]; c3 = [1, 2, 0, 3, 4]
L = lambda s: [ix[(s[a], b)] for a, b in pts]; R = lambda s: [ix[(a, s[b])] for a, b in pts]
sw = [ix[(b, a)] for a, b in pts]
out.append(fmt("S5wrS2_product", 25, [L(cyc), L(tr), R(cyc), sw]))
out.append(fmt("(A5xA5).2_product", 25, [L(cyc), L(c3), R(cyc), R(c3), sw]))

# ---------------- Mathieu groups
def cyc_perm(n, cycles):
    a = list(range(n))
    for c in cycles:
        for i in range(len(c)): a[c[i] - 1] = c[(i + 1) % len(c)] - 1
    return a
m11 = [cyc_perm(11, [list(range(1, 12))]), cyc_perm(11, [[3, 7, 11, 8], [4, 10, 5, 6]])]
out.append(fmt("M11", 11, m11))
m12 = [cyc_perm(12, [list(range(1, 12))]), cyc_perm(12, [[3, 7, 11, 8], [4, 10, 5, 6]]), cyc_perm(12, [[1, 12], [2, 11], [3, 6], [4, 8], [5, 9], [7, 10]])]
out.append(fmt("M12", 12, m12))
# PSL(2,11) on 11 points: the 11-point action inside M11 (normalizer of an 11-cycle is 11:5)
out.append(fmt("L2(11)_on_11", 11, [cyc_perm(11, [list(range(1, 12))]), cyc_perm(11, [[2, 10], [3, 4], [5, 9], [6, 7]])]))

open(sys.argv[1] if len(sys.argv) > 1 else "groups.txt", "w").write("\n".join(out) + "\n")
print(len(out), "groups written")
