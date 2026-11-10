#!/usr/bin/env python3
"""Independent verification of candidate counterexamples (shares no code with the Rust tools).

For each tuple (cycle notation, 1-based, rightmost-first composition, g1 g2 g3 g4 = 1):
  * product is 1, no entry is the identity, genus from Riemann-Hurwitz
  * generated group: order and primitivity by SymPy
  * braid orbit: closure under the 6 moves, states modulo simultaneous S_d conjugation
    (canonical form = least relabelling by breadth-first search from every base point)
  * at every state: defect at {12|34} and whether it is one-piece
usage: verify_cex.py FILE.jsonl   (lines with "primitive_nonabelian" candidates or any "tuple")
"""
import json, sys
from sympy.combinatorics import Permutation, PermutationGroup

def parse(tup):
    n = max(int(x) for p in tup for x in p.replace('(', ' ').replace(')', ' ').split())
    out = []
    for p in tup:
        a = list(range(n))
        for cyc in p.split(')'):
            v = [int(x) - 1 for x in cyc.replace('(', ' ').split()]
            for i in range(len(v)): a[v[i]] = v[(i + 1) % len(v)]
        out.append(tuple(a))
    return n, out

def mul(a, b): return tuple(a[i] for i in b)
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
def orbits(gens, n):
    lab = [-1] * n; k = 0
    for s in range(n):
        if lab[s] >= 0: continue
        st = [s]; lab[s] = k
        while st:
            x = st.pop()
            for g in gens:
                if lab[g[x]] < 0: lab[g[x]] = k; st.append(g[x])
        k += 1
    return lab, k
def cusp(t, n):
    x, y, z, w = t; p = mul(x, y)
    la, na = orbits([x, y], n); lb, nb = orbits([z, w], n)
    E = len(cycles(p)); defect = E - na - nb + 1
    def genera(lab, k, a, b, q):
        size = [0] * k
        for i in range(n): size[lab[i]] += 1
        c = [0] * k
        for g in (a, b, q):
            for C in cycles(g): c[lab[C[0]]] += 1
        return [(size[v] + 2 - c[v]) // 2 for v in range(k)]
    gs = genera(la, na, x, y, p) + genera(lb, nb, z, w, inv(p))
    pos = sum(1 for g in gs if g > 0)
    return defect, defect == 0 and pos <= 1
def move(t, m):
    j = m % 3; a, b = t[j], t[j + 1]; u = list(t)
    if m < 3: u[j] = mul(a, mul(b, inv(a))); u[j + 1] = a
    else: u[j] = b; u[j + 1] = mul(inv(b), mul(a, b))
    return tuple(u)
def canon(t, n):
    best = None
    for r in range(n):
        lab = {r: 0}; order = [r]; i = 0
        while i < len(order):
            x = order[i]; i += 1
            for g in t[:3]:
                if g[x] not in lab: lab[g[x]] = len(order); order.append(g[x])
        key = tuple(lab[t[j][order[k]]] for j in range(3) for k in range(n))
        if best is None or key < best: best = key
    return best

def verify(tup, cap=200000):
    n, t = parse(tup); t = tuple(t)
    e = tuple(range(n))
    assert mul(t[0], mul(t[1], mul(t[2], t[3]))) == e, "product is not 1"
    assert all(g != e for g in t), "identity entry"
    genus = (sum(n - len(cycles(g)) for g in t) - 2 * n + 2) // 2
    G = PermutationGroup([Permutation(list(g)) for g in t])
    seen = {canon(t, n): t}; q = [t]; defects = {}; onep = 0
    while q:
        u = q.pop()
        d, op = cusp(u, n); defects[d] = defects.get(d, 0) + 1; onep += op
        for m in range(6):
            v = move(u, m); k = canon(v, n)
            if k not in seen:
                seen[k] = v; q.append(v)
                if len(seen) > cap: return dict(d=n, capped=True)
    return dict(d=n, genus=genus, transitive=G.is_transitive(), primitive=G.is_primitive(), order=G.order(),
                states=len(seen), defect_histogram=defects, one_piece_states=onep)

if __name__ == "__main__":
    for line in open(sys.argv[1]):
        r = json.loads(line)
        if 'tuple' not in r or not r.get('primitive') or r.get('abelian'): continue
        v = verify(r['tuple'])
        ok = (not v.get('capped')) and v['primitive'] and v['genus'] >= 1 and set(v['defect_histogram']) == {0}
        print(r['group'], json.dumps(v), "COUNTEREXAMPLE CONFIRMED" if ok else "NOT confirmed", flush=True)
