#!/usr/bin/env python3
"""Exact braid-image/passport-fiber and characteristic-polynomial checks.

Consumes verified mechanism-audit tables. Reuses the project's documented
SymPy 1.14 environment; no new numerical or group implementation is introduced.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import defaultdict
import csv
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup


def compose(a, b):
    return [a[i] for i in b]


def inverse(p):
    result = [0] * len(p)
    for i, j in enumerate(p):
        result[j] = i
    return result


def equivariant_map(source, target):
    """Exact move-labelled factor search; a root image determines the map."""
    n, m = len(source[0]), len(target[0])
    if n % m:
        return None
    for root_image in range(m):
        images = {0: root_image}
        queue = [0]
        valid = True
        for i in queue:
            for g, h in zip(source, target):
                dest, image = g[i], h[images[i]]
                if dest in images and images[dest] != image:
                    valid = False
                    break
                if dest not in images:
                    images[dest] = image
                    queue.append(dest)
            if not valid:
                break
        if valid and len(images) == n and len(set(images.values())) == m:
            return [images[i] for i in range(n)]
    return None


def analyze(generators, passports):
    n = len(passports)
    identity = list(range(n))
    if len(generators) != 3 or any(sorted(g) != identity for g in generators):
        raise ValueError('invalid generator permutations')
    def word(w):
        p = identity
        for i in w:
            p = compose(generators[i], p)
        return p
    if word([0, 1, 0]) != word([1, 0, 1]) or word([1, 2, 1]) != word([2, 1, 2]) or word([0, 2]) != word([2, 0]):
        raise ValueError('braid relations fail')
    if word([0, 1] * 3) != identity or word([0, 1, 2] * 4) != identity:
        raise ValueError('product-one conjugacy-quotient central relations fail')
    fibers = defaultdict(list)
    for i, passport in enumerate(passports):
        fibers[passport].append(i)
    keys = sorted(fibers)
    owner = {i: keys.index(p) for p, points in fibers.items() for i in points}
    quotient = []
    for g in generators:
        row = []
        for p in keys:
            images = {owner[g[i]] for i in fibers[p]}
            if len(images) != 1:
                raise ValueError('passport partition is not move-equivariant')
            row.append(images.pop())
        quotient.append(row)
    if len(set(map(len, fibers.values()))) != 1:
        raise ValueError('passport fibers have unequal cardinalities')
    transports = {0: identity}
    queue = [0]
    for src in queue:
        for g, q in zip(generators, quotient):
            dest = q[src]
            if dest not in transports:
                transports[dest] = compose(g, transports[src])
                queue.append(dest)
    if len(transports) != len(keys):
        raise ValueError('disconnected passport action')
    base = fibers[keys[0]]
    base_index = {value: i for i, value in enumerate(base)}
    loops = []
    for src in queue:
        for g, q in zip(generators, quotient):
            loop = compose(inverse(transports[q[src]]), compose(g, transports[src]))
            loops.append([base_index[loop[i]] for i in base])
    group = PermutationGroup([Permutation(g) for g in generators])
    if len(group.orbits()) != 1:
        raise ValueError('disconnected braid action')
    quotient_group = PermutationGroup([Permutation(g) for g in quotient])
    fiber_group = PermutationGroup([Permutation(g) for g in loops])
    def adjacency(perms):
        size = len(perms[0])
        matrix = sp.zeros(size)
        for p in perms:
            for move in (p, inverse(p)):
                for i, j in enumerate(move):
                    matrix[i, j] += 1
        if matrix != matrix.T or any(sum(matrix.row(i)) != 6 for i in range(size)):
            raise ValueError('six-move adjacency is not symmetric regular')
        return matrix
    variable = sp.Symbol('lambda')
    polynomial = adjacency(generators).charpoly(variable).as_expr()
    position = adjacency(quotient).charpoly(variable).as_expr()
    remainder, residue = sp.div(polynomial, position, variable)
    if residue != 0:
        raise ValueError('passport factor does not divide adjacency characteristic polynomial')
    return dict(states=n, ordered_passports=len(keys), states_per_passport=len(base),
                braid_image_order=str(group.order()), passport_action_order=str(quotient_group.order()),
                passport_stabilizer_order=str(group.order() // len(keys)),
                passport_fiber_image_order=str(fiber_group.order()),
                generator_cycle_types=[sorted(Permutation(g).cycle_structure.items()) for g in generators],
                adjacency_characteristic_polynomial=str(sp.factor(polynomial)),
                passport_characteristic_polynomial=str(sp.factor(position)),
                complementary_characteristic_polynomial=str(sp.factor(remainder)),
                braid_relations_verified=True, central_relations_verified=True,
                passport_spectral_divisibility_exact=True,
                quotient='simultaneous S_d conjugacy', spectral_arithmetic='exact integer/rational polynomial; no numerical gap claim')


def selftest():
    # The regular adjacent-transposition action of S4 with distinct entry types.
    states = list(itertools.permutations(range(4)))
    index = {p: i for i, p in enumerate(states)}
    generators = []
    for slot in range(3):
        row = []
        for p in states:
            q = list(p)
            q[slot], q[slot + 1] = q[slot + 1], q[slot]
            row.append(index[tuple(q)])
        generators.append(row)
    result = analyze(generators, list(map(str, states)))
    assert (result['braid_image_order'], result['states_per_passport']) == ('24', 1)
    assert result['passport_fiber_image_order'] == '1'
    assert result['complementary_characteristic_polynomial'] == '1'
    target = [[1, 0, 2, 3], [0, 2, 1, 3], [0, 1, 3, 2]]
    factor = equivariant_map(generators, target)
    assert factor is not None and sorted({factor[i] for i in range(24)}) == list(range(4))
    assert equivariant_map(target, generators) is None
    broken = [list(g) for g in generators]
    broken[0][0] = broken[0][1]
    try:
        analyze(broken, list(map(str, states)))
    except ValueError:
        pass
    else:
        raise AssertionError('corrupt permutation accepted')
    print('weak-mechanism exact checks passed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--selftest', action='store_true')
    parser.add_argument('--extended', action='store_true', help='reduced curve genus, defect transport, exact moments and coupled fibers')
    parser.add_argument('--null-replicates', type=int, default=0)
    parser.add_argument('--seed', type=int, default=20261004)
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    if args.audit is None or args.out is None:
        parser.error('--audit and --out required')
    if not 0<=args.null_replicates<=1000:
        parser.error('null-replicates must be 0..1000')
    if args.null_replicates and not args.extended:
        parser.error('null-replicates requires extended')
    def table(name):
        with (args.audit / name).open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))
    edges, decorations = table('EDGES.tsv'), table('DECORATIONS.tsv')
    symmetries=table('SYMMETRIES.tsv') if args.extended and (args.audit/'SYMMETRIES.tsv').exists() else []
    summaries=table('SUMMARY.tsv') if args.extended and (args.audit/'SUMMARY.tsv').exists() else []
    cases = sorted({r['case'] for r in edges})
    output, actions, passports = [], {}, {}
    for case in cases:
        rows = sorted((r for r in edges if r['case'] == case), key=lambda r: int(r['state']))
        data = sorted((r for r in decorations if r['case'] == case), key=lambda r: int(r['state']))
        if [int(r['state']) for r in rows] != list(range(len(rows))) or [r['state'] for r in data] != [r['state'] for r in rows]:
            raise ValueError('missing, duplicated or misaligned state')
        generators = [[int(r[f'sigma{i + 1}']) for r in rows] for i in range(3)]
        passports[case] = [r['ordered_passport'] for r in data]
        actions[case] = generators
        result = analyze(generators, passports[case])
        result.update(case=case, distinct_defects=sorted({int(r['defect']) for r in data}),
                      pair_cusps=len({r['cusp'] for r in data}))
        if args.extended:
            import weak_graph_extensions
            result['extended']=weak_graph_extensions.analyze(generators,passports[case],
                [int(r['defect']) for r in data],result['braid_image_order'],args.null_replicates,
                args.seed^int(hashlib.sha256(case.encode()).hexdigest()[:16],16))
            relabelings={int(r['move_zero_based']):list(map(int,r['images_zero_based'].split(',')))
                         for r in symmetries if r['case']==case}
            if 0 in relabelings and 1 in relabelings:
                result['extended']['relabeling_common_quotient']=weak_graph_extensions.symmetry_common_quotient([relabelings[0],relabelings[1]])
            metadata=[r for r in summaries if r['case']==case]
            if len(metadata)==1:
                meta=metadata[0];g=int(meta['genus']);d=int(meta['degree'])
                weak=int(meta['all_six_pairs_intransitive_states'])==len(rows)
                group=Path(case).name.split('_')[1][0] if Path(case).name.startswith('CERTIFIED_') else 'unspecified'
                defects={int(r['defect']) for r in data}
                result['extended']['rare_feature_targets']=dict(weak_component=weak,degree_below_18=weak and d<18,
                    alternating_odd_defect=weak and group=='A' and any(delta%2 for delta in defects),
                    symmetric_genus_zero=weak and group=='S' and g==0,
                    positive_genus_identically_zero_defect=weak and g>0 and defects=={0},
                    variable_defect=weak and len(defects)>1,
                    no_generator_fixed_anywhere=all(int(meta[f'fixed_sigma{i}'])==0 for i in (1,2,3)))
        output.append(result)
        print(Path(case).name, result['braid_image_order'], result['passport_fiber_image_order'], flush=True)
    factors = []
    for source in cases:
        for target in cases:
            if source == target or len(actions[target][0]) == 1:
                continue
            images = equivariant_map(actions[source], actions[target])
            if images is not None:
                fibers = len(set(passports[source]))
                factors.append(dict(source=source, target=target, images=images,
                                    synchronous_product_with_passports=(len(images) == fibers * len(actions[target][0])
                                        and len(set(zip(passports[source], images))) == len(images))))
    sources = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
               [args.audit / 'EDGES.tsv', args.audit / 'DECORATIONS.tsv', Path(__file__)]}
    if args.extended:
        extension=Path(__file__).with_name('weak_graph_extensions.py')
        sources[str(extension)]=hashlib.sha256(extension.read_bytes()).hexdigest()
        if symmetries:
            source=args.audit/'SYMMETRIES.tsv';sources[str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
        if summaries:
            source=args.audit/'SUMMARY.tsv';sources[str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
    document = dict(sympy_version=sp.__version__, components=output, equivariant_factors=factors, input_sha256=sources,
                    scope=f'{len(cases)} supplied closed components; exact finite braid images and spectra, not a degree census')
    with args.out.open('x') as stream:
        json.dump(document, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
