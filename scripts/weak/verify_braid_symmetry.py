"""Independent SymPy group-order check of Rust's one-state braid certificate.

Run with uv run --isolated --no-project --with sympy==1.14.0 python ...
Array composition below is explicitly right-to-left; SymPy is used for
independently implemented exact subgroup orders, orbits and primitivity.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sympy
from sympy.combinatorics import Permutation, PermutationGroup


def compose(a, b):
    return [a[x] for x in b]


def inverse(a):
    result = [0] * len(a)
    for i, x in enumerate(a):
        result[x] = i
    return result


def conjugate(a, h):
    return compose(h, compose(a, inverse(h)))


def check(data):
    d = data['degree']
    t, h = data['tuple'], data['conjugators']
    identity = list(range(d))
    if len(t) != 4 or len(h) != 3 or any(sorted(p) != identity for p in t + h):
        raise ValueError('invalid permutation arrays')
    if compose(compose(compose(t[0], t[1]), t[2]), t[3]) != identity:
        raise ValueError('tuple product is not one')
    for i in range(3):
        moved = list(t)
        moved[i:i+2] = [conjugate(t[i+1], t[i]), t[i]]
        if [conjugate(p, h[i]) for p in t] != moved:
            raise ValueError('braid conjugator identity failed')
    for i in range(2):
        if compose(compose(h[i], h[i+1]), h[i]) != compose(compose(h[i+1], h[i]), h[i+1]):
            raise ValueError('braid relation failed')
    if compose(h[0], h[2]) != compose(h[2], h[0]):
        raise ValueError('far commutativity failed')
    delta = compose(compose(h[0], h[1]), h[2])
    if compose(compose(delta, delta), compose(delta, delta)) != identity:
        raise ValueError('full twist failed')
    group = PermutationGroup([Permutation(p) for p in t])
    order = int(group.order())
    if not all(Permutation(p).is_even for p in t) or order != math.factorial(d)//2:
        raise ValueError('tuple does not generate A_d')
    symmetry = PermutationGroup([Permutation(p) for p in h])
    symmetry_order = int(symmetry.order())
    all_odd = all(not Permutation(p).is_even for p in h)
    if order != int(data['group_order']) or symmetry_order != int(data['symmetry_order']):
        raise ValueError('Rust and SymPy group orders disagree')
    pairs = []
    for i, j in itertools.combinations(range(4), 2):
        subgroup = PermutationGroup([Permutation(t[i]), Permutation(t[j])])
        orbits = sorted([sorted(o) for o in subgroup.orbits()])
        if len(orbits) == 1:
            raise ValueError('transitive pair contradicts the claimed obstruction')
        pairs.append(dict(pair=[i,j], order=int(subgroup.order()), orbits=orbits))
    return dict(sympy_version=sympy.__version__, degree=d, tuple_order=str(order),
                symmetry_order=str(symmetry_order),
                all_braid_conjugators_odd=all_odd,
                labeled_bipartite_by_sign=all_odd,
                labeled_minimum_eigenvalue=-6 if all_odd else None,
                labeled_normalized_absolute_gap=0 if all_odd else None,
                symmetry_orbits=[sorted(o) for o in symmetry.orbits()],
                symmetry_primitive=bool(symmetry.is_primitive(randomized=False)),
                pairs=pairs, verified=True)


def selftest():
    a, b = [1,2,0,3], [1,0,3,2]
    assert compose(a,inverse(a)) == list(range(4))
    assert PermutationGroup([Permutation(a),Permutation(b)]).order() == 12
    bad = dict(degree=4, tuple=[a]*4, conjugators=[b]*3)
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('corrupt data accepted')
    print('independent verifier selftest passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--selftest',action='store_true')
    args = parser.parse_args()
    selftest()
    if not args.selftest:
        if args.input is None or args.out is None:
            parser.error('--input and --out required')
        data = args.input.read_bytes()
        result = check(json.loads(data))
        result['input_sha256'] = hashlib.sha256(data).hexdigest()
        args.out.mkdir(parents=True,exist_ok=True)
        with (args.out/'INDEPENDENT_CHECK.json').open('x') as stream:
            json.dump(result,stream,indent=2)
            stream.write('\n')
        with (args.out/'SUMMARY.txt').open('x') as stream:
            stream.write('experiment = independent SymPy braid verification\nsearch complete = true\n')
            stream.write(f"tuple order = {result['tuple_order']}\nsymmetry order = {result['symmetry_order']}\n")
            stream.write(f"labeled bipartite by sign = {result['labeled_bipartite_by_sign']}\n")
        print(json.dumps(result,indent=2))
