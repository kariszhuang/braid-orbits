#!/usr/bin/env python3
"""Deep structure of the PSL(3,q) counterexamples (q = 4, 5, 7).

For each orbit (from the stored tuple, point labelling of make_groups / groups3): reconstruct
the field and the points of PG(2,q), recover each entry's matrix, classify it (homology /
elation / other) with centre and axis, and report the incidence configuration of centres and
axes. Then act with the same matrices on LINES and compute genus, defects and one-piece.
"""
import itertools, json, sys
from fractions import Fraction
sys.path.insert(0, '.')
import make_groups as mg
from verify_cex import parse, mul, inv, cycles, move, canon, orbits, cusp

def plane(q):
    F = mg.GF(*{4: (2, 2), 5: (5, 1), 7: (7, 1), 8: (2, 3), 9: (3, 2)}[q])
    pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
    ix = {v: i for i, v in enumerate(pts)}
    def norm(v): l = next(c for c in v if c); li = F.inv[l]; return tuple(F.mul[li][c] for c in v)
    def dot(a, b):
        s = 0
        for x, y in zip(a, b): s = F.add[s][F.mul[x][y]]
        return s
    lines = [frozenset(i for i, x in enumerate(pts) if dot(v, x) == 0) for v in pts]   # line with coordinates v
    return F, pts, ix, norm, dot, lines

def matvec(F, M, v):
    out = []
    for r in range(3):
        s = 0
        for c in range(3): s = F.add[s][F.mul[M[r][c]][v[c]]]
        out.append(s)
    return tuple(out)

def recover(F, pts, ix, norm, g):
    """3x3 matrix (up to scalar) inducing the point permutation g, or None"""
    E = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    cols = [pts[g[ix[e]]] for e in E]
    q = F.q
    for a, b, c in itertools.product(range(1, q), repeat=3):
        M = [[F.mul[cols[j][i]][s] for j, s in enumerate((a, b, c))] for i in range(3)]
        if all(ix[norm(matvec(F, M, v))] == g[i] for i, v in enumerate(pts)): return M
    return None

def classify(g, lines, n):
    fix = {i for i in range(n) if g[i] == i}
    axes = [L for L in lines if L <= fix]
    order = 1; x = g
    while any(x[i] != i for i in range(n)): x = mul(g, x); order += 1
    if axes:
        A = axes[0]; extra = fix - A
        if extra: return ('homology', order, A, min(extra))
        # elation: centre = the point on the axis through which every fixed line passes
        fixed_lines = [L for L in lines if frozenset(g[i] for i in L) == L]
        through = [p for p in A if all(p in L for L in fixed_lines)]
        return ('elation', order, A, through[0] if through else None)
    return ('other', order, None, None)

def analyse(name, q, tup):
    F, pts, ix, norm, dot, lines = plane(q)
    n, t = parse(tup); t = tuple(t); assert n == len(pts)
    lab = lambda i: i + 1
    print(f"\n=== {name}: PG(2,{q}), {n} points, tuple entries:")
    info = []
    for k, g in enumerate(t):
        kind, order, axis, centre = classify(g, lines, n)
        info.append((kind, order, axis, centre))
        M = recover(F, pts, ix, norm, g)
        print(f"  g{k+1}: {kind}, order {order}, fixed {sum(1 for i in range(n) if g[i]==i)}, "
              f"axis {sorted(lab(i) for i in axis) if axis else '-'}, centre {lab(centre) if centre is not None else '-'}, matrix {M}")
    cent = [c for (_, _, _, c) in info if c is not None]
    axs = [a for (_, _, a, _) in info if a is not None]
    if len(cent) >= 3:
        for trip in itertools.combinations(range(len(cent)), 3):
            col = any({cent[i] for i in trip} <= L for L in lines)
            print(f"  centres of entries {[i+1 for i in trip]} collinear: {col}")
    if len(axs) >= 3:
        for trip in itertools.combinations(range(len(axs)), 3):
            inter = frozenset.intersection(*[axs[i] for i in trip])
            print(f"  axes of entries {[i+1 for i in trip]} concurrent: {bool(inter)}")
    for i, (k1, _, a1, c1) in enumerate(info):
        for j, (k2, _, a2, c2) in enumerate(info):
            if i < j and c1 is not None and a2 is not None:
                print(f"  centre of g{i+1} on axis of g{j+1}: {c1 in a2};  centre of g{j+1} on axis of g{i+1}: {c2 in a1 if (c2 is not None and a1 is not None) else '-'}")
    # same matrices acting on LINES
    lix = {L: i for i, L in enumerate(lines)}
    tl = tuple(tuple(lix[frozenset(g[i] for i in L)] for L in lines) for g in t)
    nl = len(lines)
    gl = (sum(nl - len(cycles(g)) for g in tl) - 2 * nl + 2) // 2
    seen = {canon(tl, nl): tl}; qq = [tl]
    while qq:
        u = qq.pop()
        for m in range(6):
            v = move(u, m); kk = canon(v, nl)
            if kk not in seen: seen[kk] = v; qq.append(v)
    st = list(seen.values()); cs = [cusp(u, nl) for u in st]
    print(f"  ON LINES: genus {gl}, cycle types {[sorted((len(c) for c in cycles(g) if len(c)>1), reverse=True)[:3] for g in tl]}, "
          f"orbit {len(st)} states, defects {sorted({d for d,_ in cs})}, one-piece states {sum(op for _,op in cs)}")
    return t, n

if __name__ == "__main__":
    ex = {}
    for f, q in [('COUNTEREXAMPLES.jsonl', 4), ('GENUS3_COUNTEREXAMPLES.jsonl', 5), ('PSL37_COUNTEREXAMPLES.jsonl', 7)]:
        for k, l in enumerate(open(f)):
            r = json.loads(l)
            if f == 'COUNTEREXAMPLES.jsonl' and r['d'] != 21: continue
            analyse(f"PSL(3,{q}) orbit {k+1} (genus {r.get('genus')})", q, r['tuple'])
