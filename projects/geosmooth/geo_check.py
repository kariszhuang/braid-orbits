"""Independent check of geometric all-smoothness for one four-point tuple.

Reuses only the orbit enumeration of ../common/braid.py (an independent
Python implementation).  The dual-graph test below is written separately
from the Rust tool.

At the {12|34} cusp: vertices = orbits of <g1,g2> and of <g3,g4>, edges =
cycles of g1 g2 (each joins the two orbits containing it).  The stable limit is
smooth iff the graph is a tree with at most one positive-genus vertex.
Usage: python3 geo_check.py D "(1 2 3)" "(...)" "(...)" "(...)"
"""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import parse, compose, inverse, cycles, index, orbit, smooth, genus


def orbit_partition(gens, n):
    lab = [-1] * n
    k = 0
    for s in range(n):
        if lab[s] >= 0:
            continue
        lab[s] = k
        todo = [s]
        while todo:
            x = todo.pop()
            for g in gens:
                y = g[x]
                if lab[y] < 0:
                    lab[y] = k
                    todo.append(y)
        k += 1
    return lab, k


def cycle_list(p):
    seen, out = set(), []
    for s in range(len(p)):
        if s in seen:
            continue
        c, x = [], s
        while x not in seen:
            seen.add(x)
            c.append(x)
            x = p[x]
        out.append(c)
    return out


def local_genus(perms, block):
    """genus of the cover of a 3-pointed sphere restricted to an orbit."""
    s = set(block)
    cyc = [sum(1 for c in cycle_list(p) if c[0] in s) for p in perms]
    two_g = len(block) + 2 - sum(cyc)
    assert two_g >= 0 and two_g % 2 == 0
    return two_g // 2


def geometric_cusp(t):
    n = len(t[0])
    a, b, c, e = t
    p = compose(a, b)
    la, na = orbit_partition([a, b], n)
    lb, nb = orbit_partition([c, e], n)
    edges = [(la[cy[0]], na + lb[cy[0]]) for cy in cycle_list(p)]
    # union-find tree test
    par = list(range(na + nb))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    acyclic = True
    for u, v in edges:
        ru, rv = f(u), f(v)
        if ru == rv:
            acyclic = False
        par[ru] = rv
    connected = len({f(x) for x in range(na + nb)}) == 1
    tree = acyclic and connected
    pinv = inverse(p)
    pos = 0
    for k in range(na):
        blk = [x for x in range(n) if la[x] == k]
        pos += local_genus([a, b, pinv], blk) > 0
    for k in range(nb):
        blk = [x for x in range(n) if lb[x] == k]
        pos += local_genus([c, e, p], blk) > 0
    return tree and pos <= 1


def main():
    n = int(sys.argv[1])
    t = parse(sys.argv[2:6], n)
    assert compose(t[0], compose(t[1], compose(t[2], t[3]))) == tuple(range(n)), 'product is not 1'
    O = orbit(t)
    geo = all(geometric_cusp(u) for u in O)
    comb = all(smooth(u) for u in O)
    print(dict(N=len(O), geometric_all_smooth=geo, combinatorial_all_smooth=comb, genus=genus(t)))


if __name__ == '__main__':
    main()
