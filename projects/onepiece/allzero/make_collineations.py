#!/usr/bin/env python3
"""All central collineations (homologies and elations) of PG(2,q), as permutations of the points,
plus one representative of each PGL(3,q)-conjugacy type (one elation; one homology per
multiplier a != 1). Output: coll_q{q}.txt
  line 1: q n
  line 2: indices of the representatives
  then one permutation per line (0-based images).
"""
import itertools, sys
sys.path.insert(0, '.')
import make_groups as mg

def build(q):
    pk = {4: (2, 2), 5: (5, 1), 7: (7, 1), 8: (2, 3), 9: (3, 2), 11: (11, 1), 13: (13, 1)}[q]
    F = mg.GF(*pk)
    pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
    ix = {v: i for i, v in enumerate(pts)}
    def norm(v): l = next(c for c in v if c); li = F.inv[l]; return tuple(F.mul[li][c] for c in v)
    def add(a, b): return F.add[a][b]
    def mulf(a, b): return F.mul[a][b]
    def act(M):
        out = []
        for v in pts:
            w = []
            for r in range(3):
                s = 0
                for c in range(3): s = add(s, mulf(M[r][c], v[c]))
                w.append(s)
            out.append(ix[norm(tuple(w))])
        return tuple(out)
    # central collineation with centre C (vector) and axis given by linear form f (vector):
    # x -> x + t f(x) C  (elation if f(C) = 0, homology otherwise); as a matrix I + t C f^T.
    perms = {}
    for C in pts:
        for f in pts:
            fC = 0
            for a, b in zip(f, C): fC = add(fC, mulf(a, b))
            for t in range(1, q):
                M = [[add(1 if r == c else 0, mulf(t, mulf(C[r], f[c]))) for c in range(3)] for r in range(3)]
                # skip singular: det(I + t C f^T) = 1 + t f(C)
                if add(1, mulf(t, fC)) == 0: continue
                p = act(M)
                if any(p[i] != i for i in range(len(pts))): perms[p] = True
    perms = list(perms)
    pindex = {p: i for i, p in enumerate(perms)}
    reps = []
    e1, e3 = (1, 0, 0), (0, 0, 1)
    I = lambda: [[1 if r == c else 0 for c in range(3)] for r in range(3)]
    M = I(); M[0][2] = 1; reps.append(pindex[act(M)])                # elation x -> x + x3 e1
    for a in range(2, q):                                            # homologies diag(a,1,1)
        M = I(); M[0][0] = a; p = act(M)
        if p in pindex and pindex[p] not in reps: reps.append(pindex[p])
    return pts, perms, reps

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        pts, perms, reps = build(q)
        with open(f"coll_q{q}.txt", "w") as fh:
            fh.write(f"{q} {len(pts)}\n" + " ".join(map(str, reps)) + "\n")
            for p in perms: fh.write(" ".join(map(str, p)) + "\n")
        print(f"q={q}: {len(pts)} points, {len(perms)} central collineations, {len(reps)} conjugacy types")
