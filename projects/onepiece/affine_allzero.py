#!/usr/bin/env python3
"""Try to disprove: "a positive-genus primitive braid orbit with defect 0 at every state is C_p".

Exhaustive over the affine groups G = C_p x| C_k <= AGL(1,p) (primitive, p prime), k | p-1.
Elements are pairs (a, b): x -> a x + b, composition (a,b)(c,d) = (ac, ad+b) (apply right first),
matching the project's rightmost-first convention. For every identity-free product-one
4-tuple generating a transitive group, braid orbits are found by union-find over the
three moves sigma_1..sigma_3 (plus inverses implicitly, orbits of a group action), states
modulo nothing (labelled); conjugation by G only identifies states inside an orbit-union so
orbits are unions of conjugacy classes automatically. For each orbit: genus and whether
every state has defect 0 at {12|34} (every cusp of every state is {12|34} of some state).
usage: affine_allzero.py P [P ...]
"""
import sys, collections

def run(p):
    ks = [k for k in range(1, p) if (p - 1) % k == 0]
    # multiplicative generator
    g = next(r for r in range(2, p) if all(pow(r, (p - 1) // q, p) != 1 for q in range(2, p) if (p - 1) % q == 0 and all(q % s for s in range(2, q))))
    out = []
    for k in ks:
        sub = sorted({pow(g, (p - 1) // k * i, p) for i in range(k)})  # C_k in F_p^*
        G = [(a, b) for a in sub for b in range(p)]
        idx = {x: i for i, x in enumerate(G)}
        mul = lambda x, y: (x[0] * y[0] % p, (x[0] * y[1] + x[1]) % p)
        inv = lambda x: (pow(x[0], p - 2, p), (-pow(x[0], p - 2, p) * x[1]) % p)
        e = (1, 0)
        def perm(x): return tuple((x[0] * i + x[1]) % p for i in range(p))
        def ncyc(x):
            a, b = x
            if a == 1: return p if b == 0 else 1
            r = 1; t = a
            while t != 1: t = t * a % p; r += 1
            return 1 + (p - 1) // r
        def norbits(x, y):
            # orbits of <x,y> on F_p
            seen = [False] * p; cnt = 0; px, py = perm(x), perm(y)
            for s in range(p):
                if seen[s]: continue
                cnt += 1; st = [s]; seen[s] = True
                while st:
                    u = st.pop()
                    for v in (px[u], py[u]):
                        if not seen[v]: seen[v] = True; st.append(v)
            return cnt
        def defect(t):
            return ncyc(mul(t[0], t[1])) - norbits(t[0], t[1]) - norbits(t[2], t[3]) + 1
        def ind(x): return p - ncyc(x)
        tuples = []
        for x in G:
            if x == e: continue
            for y in G:
                if y == e: continue
                xy = mul(x, y)
                for z in G:
                    if z == e: continue
                    w = inv(mul(xy, z))
                    if w == e: continue
                    t = (x, y, z, w)
                    # transitive: not all fixing a common point  <=> group contains a translation or fixed points differ
                    fixed = set()
                    trans = False
                    for u in t:
                        if u[0] == 1: trans = True; break
                        fixed.add(u[1] * pow(1 - u[0], p - 2, p) % p)
                    if not trans and len(fixed) == 1: continue
                    tuples.append(t)
        tid = {t: i for i, t in enumerate(tuples)}
        par = list(range(len(tuples)))
        def f(i):
            while par[i] != i: par[i] = par[par[i]]; i = par[i]
            return i
        for i, t in enumerate(tuples):
            for j in range(3):
                a, b = t[j], t[j + 1]
                u = list(t); u[j] = mul(a, mul(b, inv(a))); u[j + 1] = a
                par[f(i)] = f(tid[tuple(u)])
        orbs = collections.defaultdict(list)
        for i in range(len(tuples)): orbs[f(i)].append(i)
        stats = collections.Counter(); hits = []
        for r, mem in orbs.items():
            t0 = tuples[mem[0]]
            gen = (sum(ind(u) for u in t0) - 2 * p + 2) // 2
            allzero = all(defect(tuples[i]) == 0 for i in mem)
            stats[(gen > 0, allzero)] += 1
            if gen > 0 and allzero:
                # the generated group is non-cyclic iff some entry is not a translation
                noncyc = any(u[0] != 1 for u in t0)
                hits.append((gen, len(mem), [(u[0], u[1]) for u in t0], 2 if noncyc else 1))
        out.append((k, len(tuples), len(orbs), dict(stats), hits))
    return out

if __name__ == "__main__":
    for p in map(int, sys.argv[1:]):
        for k, nt, no, stats, hits in run(p):
            pos_allzero = [h for h in hits if h[3] > 1]
            cyc = [h for h in hits if h[3] == 1]
            print(f"p={p} k={k}: tuples {nt}, orbits {no}, (pos_genus, allzero) counts {stats}; positive-genus all-zero: {len(hits)} (generating C_p: {len(cyc)}, generating a NON-CYCLIC group: {len(pos_allzero)})", flush=True)
            for h in pos_allzero[:5]: print("   COUNTEREXAMPLE? genus", h[0], "labelled states", h[1], "tuple (a,b):", h[2])
