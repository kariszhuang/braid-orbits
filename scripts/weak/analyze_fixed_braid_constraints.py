#!/usr/bin/env python3
"""Literal audit of fixed-braid constraints in existing weak certificates.

No search or full-group enumeration. Generation/closure certification is reused
from the supplied certificate catalog. Mathematical proofs are separate from
these finite checks. Maps in JSON use one-based point labels.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import random

from odd_order_spin import compose, cycle_lengths, hurwitz, inverse, parse_seed


def require(condition, message):
    if not condition:
        raise ValueError(message)


def product(*perms):
    require(bool(perms), "empty product")
    result = tuple(range(len(perms[0])))
    for p in perms:
        result = compose(result, p)
    return result


def power(p, exponent):
    result = tuple(range(len(p)))
    for _ in range(exponent):
        result = compose(result, p)
    return result


def conjugate(p, h):
    return product(h, p, inverse(h))


def point_orbits(perms):
    n = len(perms[0])
    unseen = set(range(n))
    result = []
    while unseen:
        root = min(unseen)
        orbit = [root]
        unseen.remove(root)
        for x in orbit:
            for p in perms:
                if p[x] in unseen:
                    unseen.remove(p[x])
                    orbit.append(p[x])
        result.append(sorted(orbit))
    return result


def tuple_conjugator(source, target):
    """Root propagation, O(n^2); reject nonunique maps rather than guess.

    A full A_n/S_n tuple (n >= 4) has trivial simultaneous centralizer.
    Each root image determines the whole equivariant map by transitivity.
    """
    n = len(source[0])
    require(len(point_orbits(source)) == 1, "source is not transitive")
    answers = []
    for root_image in range(n):
        image = [None] * n
        image[0] = root_image
        queue = [0]
        consistent = True
        for x in queue:
            for a, b in zip(source, target):
                xx, yy = a[x], b[image[x]]
                if image[xx] is None:
                    image[xx] = yy
                    queue.append(xx)
                elif image[xx] != yy:
                    consistent = False
                    break
            if not consistent:
                break
        if not consistent or len(queue) != n or len(set(image)) != n:
            continue
        h = tuple(image)
        require(all(conjugate(a, h) == b for a, b in zip(source, target)),
                "propagated map failed literal conjugation")
        answers.append(h)
    require(len(answers) <= 1, "nontrivial simultaneous centralizer: out of scope")
    return answers[0] if answers else None


def incidence(t):
    left, right = point_orbits(t[:2]), point_orbits(t[2:])
    boundary = product(*t[:2])
    edges = []
    for cycle in point_orbits([boundary]):
        li = next(i for i, block in enumerate(left) if cycle[0] in block)
        ri = next(i for i, block in enumerate(right) if cycle[0] in block)
        require(set(cycle) <= set(left[li]) and set(cycle) <= set(right[ri]),
                "boundary cycle straddles a pair orbit")
        edges.append({"left": li, "right": ri, "points": [x + 1 for x in cycle]})
    degrees = [sum(e["left"] == i for e in edges) for i in range(len(left))]
    degrees += [sum(e["right"] == i for e in edges) for i in range(len(right))]
    adjacency = [set() for _ in degrees]
    for e in edges:
        i, j = e["left"], len(left) + e["right"]
        adjacency[i].add(j)
        adjacency[j].add(i)
    reached, queue = {0}, [0]
    for x in queue:
        for y in adjacency[x] - reached:
            reached.add(y)
            queue.append(y)
    require(len(reached) == len(degrees), "incidence graph disconnected")
    return {"left_sizes": sorted(map(len, left), reverse=True),
            "right_sizes": sorted(map(len, right), reverse=True),
            "edges": edges, "vertex_degrees": degrees,
            "defect": len(edges) - len(degrees) + 1}


def maps(**perms):
    return {name: [x + 1 for x in p] for name, p in perms.items()}


def audit_state(t):
    n = len(t[0])
    identity = tuple(range(n))
    require(all(sorted(p) == list(identity) for p in t), "invalid permutation")
    require(product(*t) == identity, "not product one")
    require(all(p != identity for p in t), "identity-entry/triple stratum is out of scope")
    fixed = [tuple_conjugator(t, hurwitz(t, i)) for i in range(3)]
    adjacent, ends = [], None
    for i in range(2):
        h, k = fixed[i:i + 2]
        if h is None or k is None:
            continue
        a, b, c = t[i:i + 3]
        z = product(a, b, c)
        u = product(a, inverse(h), inverse(k))
        require(product(h, k, h) == product(k, h, k), "braid relabeling relation")
        require(z == power(product(h, k), 3), "full twist identity")
        require(power(u, 3) == identity and compose(u, h) == compose(k, u),
                "order-three conjugator identity")
        require(a == product(u, k, h) and compose(a, k) == compose(k, a),
                "adjacent constraint reconstruction")
        require(b == conjugate(a, inverse(h)) and c == conjugate(b, inverse(k)),
                "adjacent entries reconstruction")
        require(all((n - len(cycle_lengths(p))) % 2 == 0 for p in t),
                "adjacent-fixed parity obstruction violated")
        require(compose(u, z) != compose(z, u), "central full twist obstruction")
        blocks = point_orbits([h, k])
        if len(blocks) > 1:
            require(any(set(u[x] for x in block) != set(block) for block in blocks),
                    "intransitive construction envelope")
        adjacent.append({"positions": [i + 1, i + 2],
                         "relabeling_orbit_sizes": sorted(map(len, blocks), reverse=True),
                         "u_cycle_type": sorted(cycle_lengths(u), reverse=True),
                         "maps": maps(h=h, k=k, u=u, z=z)})
    h, k = fixed[0], fixed[2]
    if h is not None and k is not None:
        a, b, c, d = t
        p = compose(a, b)
        u, v = compose(a, inverse(h)), compose(c, inverse(k))
        require(compose(h, k) == compose(k, h), "end relabelings do not commute")
        require(all(compose(p, q) == compose(q, p) for q in (h, k)),
                "boundary does not commute with relabelings")
        require(compose(u, k) == compose(k, u) and compose(v, h) == compose(h, v),
                "root centralizer constraint")
        require(power(u, 2) == product(p, inverse(h), inverse(h)), "left root")
        require(power(v, 2) == product(inverse(p), inverse(k), inverse(k)), "right root")
        reconstructed = (compose(u, h), compose(inverse(h), product(u, h, h)),
                         compose(v, k), compose(inverse(k), product(v, k, k)))
        require(reconstructed == t, "end-root tuple reconstruction")
        ends = {"maps": maps(h=h, k=k, p=p, u=u, v=v),
                "relabeling_orbit_sizes": sorted(map(len, point_orbits([h, k])), reverse=True),
                "elementary_abelian_three_frame": power(h, 3) == power(k, 3) == identity
                    and p == product(power(h, 2), k),
                "cycle_types": {name: sorted(cycle_lengths(q), reverse=True)
                                for name, q in (("h", h), ("k", k), ("p", p))}}
    return {"fixed_mask": "".join(str(i + 1) for i, h in enumerate(fixed) if h) or "-",
            "adjacent": adjacent, "ends": ends}


def selftest():
    rng = random.Random(20261011)
    checked = 0
    for n in range(4, 8):
        for _ in range(16):
            first = []
            for _ in range(3):
                p = list(range(n))
                rng.shuffle(p)
                first.append(tuple(p))
            t = tuple(first + [inverse(product(*first))])
            for offset in (0, 1):
                moved = t
                for mv in [offset, offset + 1] * 3:
                    moved = hurwitz(moved, mv)
                z = product(*t[offset:offset + 3])
                require(moved == tuple(conjugate(p, z) for p in t), "literal full twist")
                checked += 1
    # Compare root propagation to brute-force simultaneous relabeling in S4.
    t = parse_seed('[(1 2 3 4), (1 2), (2 3 4), (1)]', 4)
    for relabeling in itertools.permutations(range(4)):
        target = tuple(conjugate(p, relabeling) for p in t)
        require(tuple_conjugator(t, target) == relabeling, "root propagation control")
    wrong = (t[1], t[0], t[2], t[3])
    brute = [q for q in itertools.permutations(range(4))
             if tuple(conjugate(p, q) for p in t) == wrong]
    require(not brute and tuple_conjugator(t, wrong) is None, "nonconjugate control")
    # Bounded algebra controls, not a new degree-four weak search.
    perms = list(itertools.permutations(range(4)))
    order_three = [u for u in perms if power(u, 3) == tuple(range(4))]
    rigid_checks = 0
    for h in perms:
        for k in perms:
            if product(h, k, h) != product(k, h, k):
                continue
            for u in order_three:
                if compose(u, h) != compose(k, u):
                    continue
                a = product(u, k, h)
                b = conjugate(a, inverse(h))
                c = conjugate(b, inverse(k))
                t = (a, b, c, inverse(product(a, b, c)))
                require(tuple(conjugate(x, h) for x in t) == hurwitz(t, 0)
                        and tuple(conjugate(x, k) for x in t) == hurwitz(t, 1),
                        "rigid algebra control")
                rigid_checks += 1
    end_checks = 0
    perms = list(itertools.permutations(range(3)))
    for h in perms:
        for k in perms:
            if compose(h, k) != compose(k, h):
                continue
            for p in perms:
                if any(compose(p, q) != compose(q, p) for q in (h, k)):
                    continue
                left = [u for u in perms if power(u, 2) == product(p, inverse(h), inverse(h))
                        and compose(u, k) == compose(k, u)]
                right = [v for v in perms if power(v, 2) == product(inverse(p), inverse(k), inverse(k))
                         and compose(v, h) == compose(h, v)]
                for u, v in itertools.product(left, right):
                    t = (compose(u, h), product(inverse(h), u, h, h),
                         compose(v, k), product(inverse(k), v, k, k))
                    require(product(*t) == tuple(range(3)), "end-root product control")
                    require(tuple(conjugate(x, h) for x in t) == hurwitz(t, 0)
                            and tuple(conjugate(x, k) for x in t) == hurwitz(t, 2),
                            "end-root fixed-move control")
                    end_checks += 1
    counts = dict(full_twist=checked, conjugator=25, rigid=rigid_checks, end_roots=end_checks)
    print('algebra selftests passed:', counts)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificates', type=Path, default=Path('evidence/certificates'))
    parser.add_argument('--out', type=Path)
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        counts = selftest()
        if args.out is not None:
            with args.out.open('x') as stream:
                json.dump(dict(counts=counts, implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                               scope='bounded convention/algebra controls, not weak discovery'), stream, indent=2)
                stream.write('\n')
        return
    if args.out is None:
        parser.error('--out is required (new file only)')
    sources, cases = [], []
    for file in sorted(args.certificates.glob('CERTIFIED_*_weak.txt')):
        raw = file.read_text()
        fields = dict(line.split('=', 1) for line in raw.splitlines() if '=' in line)
        n = int(fields['d'])
        require(4 <= n <= 25, "only bit-packed degree 4--25 supported")
        bits = 4 if n <= 16 else 5
        rows = []
        for line in raw.splitlines():
            if not line.startswith('STATE\t'):
                continue
            words = line.split('\t')[1:]
            require(len(words) == 3, "invalid state record")
            first = [tuple((int(w, 16) >> (bits * i)) & ((1 << bits) - 1)
                           for i in range(n)) for w in words]
            rows.append(audit_state(tuple(first + [inverse(product(*first))])))
        require(len(rows) == int(fields['states']), "certificate state-count mismatch")
        seed = parse_seed(fields['seed'], n)
        seed_profile = audit_state(seed)
        graph = incidence(seed)
        if file.name == 'CERTIFIED_S18_XXXX_genus1_weak.txt':
            require(graph['left_sizes'] == graph['right_sizes'] == [6, 6, 6], "seed partitions")
            require(graph['vertex_degrees'] == [2] * 6 and len(graph['edges']) == 6,
                    "ten-state S18 seed is not a six-cycle")
            require(all(len(e['points']) == 3 for e in graph['edges']), "edge weights")
        sources.append({"path": str(file), "sha256": hashlib.sha256(file.read_bytes()).hexdigest()})
        cases.append({"certificate": file.name, "group": fields['group'], "degree": n,
                      "states": len(rows), "fixed_masks": dict(sorted(Counter(r['fixed_mask'] for r in rows).items())),
                      "adjacent_fixed_states": sum(bool(r['adjacent']) for r in rows),
                      "end_fixed_states": sum(r['ends'] is not None for r in rows),
                      "seed": seed_profile, "seed_incidence": graph, "state_checks": rows})
    require(bool(cases), "empty catalog")
    document = {"scope": "literal fixed-braid algebra on existing certificates, not new discovery or independent generation/closure verification",
                "quotient": "simultaneous S_d conjugation", "human_mathematical_review": "pending",
                "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "permutation_helper_sha256": hashlib.sha256(Path(__file__).with_name('odd_order_spin.py').read_bytes()).hexdigest(),
                "sources": sources, "cases": cases,
                "total_states": sum(c['states'] for c in cases),
                "six_three_cycle_square_roots": {
                    "by_paired_cycle_count": [math.factorial(6) //
                        (math.factorial(6 - 2*j) * 2**j * math.factorial(j)) * 3**j
                        for j in range(4)],
                    "meaning": "labeled root counts, not distinct generating states or weak components"}}
    with args.out.open('x') as stream:
        json.dump(document, stream, indent=2)
        stream.write('\n')
    for c in cases:
        print(c['certificate'], c['fixed_masks'], 'adjacent=', c['adjacent_fixed_states'], 'ends=', c['end_fixed_states'])


if __name__ == '__main__':
    main()
