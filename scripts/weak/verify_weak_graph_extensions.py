#!/usr/bin/env python3
"""Replay saved structural claims from raw permutations, not analysis helpers.

Exact group orders still use pinned SymPy, as in the analysis. Explicit maps,
quotient tables, stable defect partitions and Newton trace identities are
checked by separate code. This verifies supplied finite graphs, not discovery
coverage, a literature claim, or a theorem about all weak components.
"""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup


def require(condition, message):
    if not condition:
        raise ValueError(message)


def multiply(p, q):
    return tuple(p[q[i]] for i in range(len(q)))


def invert(p):
    result = [0]*len(p)
    for i, j in enumerate(p):
        result[j] = i
    return tuple(result)


def cycle_lengths(p):
    seen, lengths = set(), []
    for i in range(len(p)):
        if i in seen:
            continue
        length, x = 0, i
        while x not in seen:
            seen.add(x); length += 1; x = p[x]
        lengths.append(length)
    return sorted(lengths)


def parity(p):
    return (-1)**(len(p)-len(cycle_lengths(p)))


def all_moves(generators):
    return list(generators)+[invert(g) for g in generators]


def verify_reduction(generators, saved):
    q1, q2, q3 = generators
    identity = tuple(range(len(q1)))
    r = multiply(q1, invert(q3)); s = multiply(multiply(q1, q2), q3)
    t = multiply(s, s)
    klein = {identity, r, t, multiply(r, t)}
    require(multiply(r, t) == multiply(t, r) and all(multiply(k, k) == identity for k in klein), 'Klein presentation')
    require(all(multiply(multiply(g, k), invert(g)) in klein for g in generators for k in klein), 'Klein normality')
    blocks = sorted({tuple(sorted(k[x] for k in klein)) for x in range(len(q1))})
    # Deduplicate each point orbit; a nonfaithful Klein image has stabilizers.
    blocks = sorted({tuple(sorted(set(block))) for block in blocks})
    require([list(b) for b in blocks] == saved['blocks'], 'reduction blocks')
    require(len(klein) == saved['reduction_image_order'], 'reduction image order')
    owner = {x: i for i, block in enumerate(blocks) for x in block}
    induced = [tuple(owner[g[block[0]]] for block in blocks) for g in generators]
    require(all({owner[g[x]] for x in block} == {owner[g[block[0]]]} for g in generators for block in blocks), 'reduction equivariance')
    a, b, _ = induced
    branch = [multiply(a, b), multiply(multiply(a, b), a), b]
    rid = tuple(range(len(blocks)))
    require(multiply(multiply(branch[0], branch[0]), branch[0]) == rid, 'order-three branch')
    require(multiply(branch[1], branch[1]) == rid, 'order-two branch')
    require(multiply(multiply(branch[0], branch[1]), branch[2]) == rid, 'branch product')
    reached, queue = {0}, [0]
    for x in queue:
        for g in branch:
            if g[x] not in reached:
                reached.add(g[x]); queue.append(g[x])
    require(len(reached) == len(blocks), 'branch transitivity')
    types = [cycle_lengths(p) for p in branch]
    twice_genus = 2+len(blocks)-sum(map(len, types))
    require(twice_genus >= 0 and twice_genus % 2 == 0, 'curve genus integrality')
    require([list(g) for g in branch] == saved['branch_cycles'] and types == saved['branch_cycle_types'], 'branch cycles')
    require(len(blocks) == saved['reduced_states'] and twice_genus//2 == saved['reduced_curve_genus'], 'reduced genus')


def verify_moments(generators, saved):
    """Newton sums from charpoly, independent of rooted walk propagation."""
    n = len(generators[0]); adjacency = sp.zeros(n)
    for g in all_moves(generators):
        for i, j in enumerate(g):
            adjacency[i, j] += 1
    require(adjacency == adjacency.T and all(sum(adjacency.row(i)) == 6 for i in range(n)), 'regular symmetric adjacency')
    coefficients = [int(c) for c in adjacency.charpoly().all_coeffs()]
    traces = [n]
    for k in range(1, len(saved['adjacency_traces'])):
        total = sum(coefficients[j]*traces[k-j] for j in range(1, min(k-1, n)+1))
        if k <= n:
            total += k*coefficients[k]
        traces.append(-total)
    require(traces == saved['adjacency_traces'], 'Newton trace identities')
    require([str(Fraction(value, n*6**k)) for k, value in enumerate(traces)] == saved['mean_random_walk_return_probabilities'], 'rational return probabilities')


def verify_defect_partition(generators, defects, saved):
    n = len(defects); moves = all_moves(generators)
    labels = list(defects)
    while True:
        signatures = [(labels[i], tuple(sorted(Counter(labels[g[i]] for g in moves).items()))) for i in range(n)]
        keys = {value: j for j, value in enumerate(sorted(set(signatures)))}
        refined = [keys[value] for value in signatures]
        if len(set(refined)) == len(set(labels)):
            break
        labels = refined
    blocks = [[i for i in range(n) if labels[i] == label] for label in sorted(set(labels))]
    require(sorted(blocks) == sorted(saved['refined_blocks']), 'coarsest defect refinement')
    owner = {x: i for i, block in enumerate(saved['refined_blocks']) for x in block}
    require(sorted(owner) == list(range(n)), 'partition coverage')
    for j, block in enumerate(saved['refined_blocks']):
        require(len({defects[i] for i in block}) == 1, 'refinement respects defect')
        for i in block:
            counts = Counter(owner[g[i]] for g in moves)
            require([counts[k] for k in range(len(saved['refined_blocks']))] == saved['six_move_lumped_adjacency'][j], 'strong Markov quotient row')
    grouped = defaultdict(list)
    for i, defect in enumerate(defects):
        grouped[defect].append(i)
    lumpable = all(len({tuple(sorted(Counter(defects[g[i]] for g in moves).items())) for i in block}) == 1 for block in grouped.values())
    labelled = all(len({defects[g[i]] for i in block}) == 1 for g in moves for block in grouped.values())
    require(lumpable == saved['defect_strongly_lumpable'] and labelled == saved['defect_move_labelled_equivariant'], 'defect lumpability status')
    histogram = {int(k): v for k, v in saved['exact_defect_histogram'].items()}
    require(dict(Counter(defects)) == histogram, 'defect histogram')
    transport = Counter((defects[i], defects[g[i]]) for i in range(n) for g in moves)
    require(transport == Counter({(r['source'], r['target']): r['count'] for r in saved['directed_defect_transport']}), 'directed defect transport')


def verify_coupled_embedding(generators, passports, saved, group_order):
    if not saved.get('identified'):
        return False
    keys = sorted(set(passports)); fibers = saved['fiber_points_in_coordinate_order']
    require(keys == saved['ordered_passport_keys'] and len(fibers) == 6, 'coupled passport keys')
    require(all(len(b) == 10 and {x for x, p in enumerate(passports) if p == key} == set(b) for key, b in zip(keys, fibers)), 'fiber coordinate bijections')
    owner = {x: (i, j) for i, block in enumerate(fibers) for j, x in enumerate(block)}
    pairs = saved['kernel_equivalent_fiber_pairs']
    require(len(pairs) == 3 and sorted(i for pair in pairs for i in pair) == list(range(6)), 'complementary fiber pairing')
    pair_owner = {fiber: i for i, pair in enumerate(pairs) for fiber in pair}
    quotient_generators = []
    for g, coordinate in zip(generators, saved['generator_coordinates']):
        q = coordinate['passport_permutation']; pair_perm = coordinate['pair_permutation']; internal = coordinate['internal_permutations']
        require(sorted(q) == list(range(6)) and sorted(pair_perm) == list(range(3)), 'coordinate quotient permutations')
        require(all(sorted(p) == list(range(10)) for p in internal) and len(internal) == 3, 'coordinate internal permutations')
        for fiber, block in enumerate(fibers):
            require(pair_owner[q[fiber]] == pair_perm[pair_owner[fiber]], 'paired quotient action')
            for j, x in enumerate(block):
                require(owner[g[x]] == (q[fiber], internal[pair_owner[fiber]][j]), 'literal sixty-point reconstruction')
        require(math.prod(parity(p) for p in internal) == parity(pair_perm), 'coordinate parity character')
        quotient_generators.append(Permutation(q))
    require(len(saved['generator_coordinates']) == 3, 'coordinate generator count')
    require(PermutationGroup(quotient_generators).order() == 24, 'faithful S4 passport action')
    expected = 12*math.factorial(10)**3
    require(group_order == expected == int(saved['candidate_order']), 'index-two subgroup order')
    return True


def verify_common_quotient(generators, saved):
    if saved.get('status') != 'computed':
        return False
    group = PermutationGroup([Permutation(g) for g in generators])
    require(group.order() <= 65536, 'explicit quotient verification cap')
    elements = [tuple(p.array_form) for p in group.generate_schreier_sims()]
    blocks = sorted(tuple(sorted(block)) for block in group.orbits())
    require([list(b) for b in blocks] == saved['point_blocks'] and len(blocks) == 2, 'symmetry point blocks')
    def restrict(p, block):
        lookup = {x: i for i, x in enumerate(block)}
        return tuple(lookup[p[x]] for x in block)
    pairs = [(restrict(p, blocks[0]), restrict(p, blocks[1])) for p in elements]
    left, right = set(a for a, _ in pairs), set(b for _, b in pairs)
    lid, rid = tuple(range(len(blocks[0]))), tuple(range(len(blocks[1])))
    lk, rk = {a for a, b in pairs if b == rid}, {b for a, b in pairs if a == lid}
    def coset(a):
        return min(multiply(a, k) for k in lk)
    representatives = sorted({coset(a) for a in left}); lookup = {a: i for i, a in enumerate(representatives)}
    table = [[lookup[coset(multiply(a, b))] for b in representatives] for a in representatives]
    require(table == saved['quotient_multiplication_table'], 'common quotient multiplication')
    require([len(left), len(right)] == saved['block_image_orders'] and [len(lk), len(rk)] == saved['block_common_quotient_kernel_orders'], 'Goursat images and kernels')
    require(len(elements) == int(saved['image_order']) and len(representatives) == saved['common_quotient_order'], 'Goursat orders')
    require(len(elements)*len(representatives) == len(left)*len(right), 'Goursat order identity')
    identity = lookup[coset(lid)]
    images = [lookup[coset(restrict(g, blocks[0]))] for g in generators]
    require(identity == saved['quotient_identity'] and images == saved['generator_quotient_images'], 'quotient generator images')
    orders = []
    for image in images:
        value = identity
        for order in range(1, len(table)+1):
            value = table[value][image]
            if value == identity:
                orders.append(order); break
        else:
            raise ValueError('invalid quotient element order')
    abelian = all(table[i][j] == table[j][i] for i in range(len(table)) for j in range(len(table)))
    require(orders == saved['generator_quotient_orders'] and abelian == saved['common_quotient_abelian'], 'quotient orders and commutation')
    if len(table) == 6:
        require(saved['common_quotient'] == ('C6' if abelian else 'S3'), 'named order-six quotient')
    return True


def selftest():
    # C3 with three equal generators: eigenvalues 6,-3,-3.
    generators = [(1, 2, 0)]*3
    saved = dict(adjacency_traces=[3]+[6**k+2*(-3)**k for k in range(1, 13)],
                 mean_random_walk_return_probabilities=[str(Fraction((3 if k == 0 else 6**k+2*(-3)**k), 3*6**k)) for k in range(13)])
    verify_moments(generators, saved)
    saved['adjacency_traces'][3] += 1
    try:
        verify_moments(generators, saved)
    except ValueError:
        pass
    else:
        raise AssertionError('corrupt moment accepted')
    print('separate structural verifier checks passed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path); parser.add_argument('--analysis', type=Path)
    parser.add_argument('--out', type=Path); parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        selftest(); return
    if any(p is None for p in (args.audit, args.analysis, args.out)):
        parser.error('audit, analysis and out required')
    def rows(name):
        with (args.audit/name).open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))
    edges, decorations, symmetries = rows('EDGES.tsv'), rows('DECORATIONS.tsv'), rows('SYMMETRIES.tsv')
    document = json.loads(args.analysis.read_text()); results = []
    require(len(document['components']) == len({r['case'] for r in edges}), 'supplied component coverage')
    for saved in document['components']:
        case = saved['case']
        edge_rows = sorted((r for r in edges if r['case'] == case), key=lambda r: int(r['state']))
        data = sorted((r for r in decorations if r['case'] == case), key=lambda r: int(r['state']))
        require([int(r['state']) for r in edge_rows] == list(range(saved['states'])) and [r['state'] for r in data] == [r['state'] for r in edge_rows], 'raw state alignment')
        generators = [tuple(int(r[f'sigma{i}']) for r in edge_rows) for i in (1, 2, 3)]
        require(all(sorted(g) == list(range(saved['states'])) for g in generators), 'raw graph permutations')
        passports = [r['ordered_passport'] for r in data]; defects = [int(r['defect']) for r in data]
        order = int(PermutationGroup([Permutation(g) for g in generators]).order())
        require(order == int(saved['braid_image_order']), 'exact braid image order')
        extension = saved['extended']
        verify_reduction(generators, extension['reduced_hurwitz'])
        verify_moments(generators, extension['trace_moments'])
        verify_defect_partition(generators, defects, extension['defect_transport'])
        coupled = verify_coupled_embedding(generators, passports, extension['coupled_fibers'], order)
        common = False
        if 'relabeling_common_quotient' in extension:
            relabelings = {int(r['move_zero_based']): tuple(map(int, r['images_zero_based'].split(','))) for r in symmetries if r['case'] == case}
            common = verify_common_quotient([relabelings[0], relabelings[1]], extension['relabeling_common_quotient'])
        results.append(dict(case=case, states=len(data), reduced_curve_verified=True,
                            exact_moments_verified=True, defect_refinement_verified=True,
                            coupled_embedding_verified=coupled, common_quotient_verified=common))
        print(Path(case).name, 'structural maps verified', flush=True)
    sources = [args.analysis, args.audit/'EDGES.tsv', args.audit/'DECORATIONS.tsv', args.audit/'SYMMETRIES.tsv', Path(__file__)]
    output = dict(sympy_version=sp.__version__, components=results,
                  source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  scope='separate explicit finite-map/partition/Newton checks; group orders reuse exact SymPy, not a separate group implementation; human review pending')
    with args.out.open('x') as stream:
        json.dump(output, stream, indent=2); stream.write('\n')


if __name__ == '__main__':
    main()
