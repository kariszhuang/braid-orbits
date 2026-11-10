#!/usr/bin/env python3
"""Exact saved-graph investigations; no discovery search or floating gaps.

The reduced-curve calculation verifies its finite permutation presentation.
It is not the Riemann-Hurwitz genus of the original degree-d cover.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from collections import Counter, defaultdict
import argparse
from fractions import Fraction
import itertools
import math
import random
import hashlib
import json
from pathlib import Path

from analyze_weak_mechanisms import compose, inverse, equivariant_map


def power(p, n):
    result = list(range(len(p)))
    while n:
        if n & 1:
            result = compose(p, result)
        p = compose(p, p)
        n >>= 1
    return result


def cycles(p):
    unseen = set(range(len(p)))
    result = []
    while unseen:
        root = min(unseen)
        cycle, current = [], root
        while current in unseen:
            unseen.remove(current)
            cycle.append(current)
            current = p[current]
        result.append(cycle)
    return result


def generated_elements(generators, cap):
    """Explicit small image only. A capped image is never called complete."""
    identity = tuple(range(len(generators[0])))
    seen, queue = {identity}, [identity]
    for p in queue:
        for g in generators:
            q = tuple(compose(g, p))
            if q not in seen:
                if len(queue) == cap:
                    return None
                seen.add(q)
                queue.append(q)
    return queue


def induced(generators, blocks):
    owner = {point: i for i, block in enumerate(blocks) for point in block}
    output = []
    for g in generators:
        row = []
        for block in blocks:
            images = {owner[g[x]] for x in block}
            if len(images) != 1:
                raise ValueError('partition is not invariant')
            row.append(images.pop())
        if sorted(row) != list(range(len(blocks))):
            raise ValueError('induced action is not a permutation')
        output.append(row)
    return output


def reduced_hurwitz(generators):
    q1, q2, q3 = generators
    n, identity = len(q1), list(range(len(q1)))
    r = compose(q1, inverse(q3))
    s = compose(compose(q1, q2), q3)
    klein = generated_elements([r, power(s, 2)], 4)
    if klein is None or len(klein) not in (1, 2, 4):
        raise ValueError('reduction subgroup is not an image of Klein four')
    members = set(klein)
    for k in klein:
        if power(list(k), 2) != identity:
            raise ValueError('reduction subgroup is not elementary abelian')
        for g in generators:
            if tuple(compose(compose(g, k), inverse(g))) not in members:
                raise ValueError('reduction subgroup is not normal')
    unseen, blocks = set(range(n)), []
    while unseen:
        root = min(unseen)
        block = sorted({k[root] for k in klein})
        blocks.append(block)
        unseen.difference_update(block)
    a, b, _ = induced(generators, blocks)
    gamma0, gamma1, gamma_inf = compose(a, b), compose(compose(a, b), a), b
    m, reduced_identity = len(blocks), list(range(len(blocks)))
    if power(gamma0, 3) != reduced_identity or power(gamma1, 2) != reduced_identity:
        raise ValueError('reduced modular presentation fails')
    if compose(compose(gamma0, gamma1), gamma_inf) != reduced_identity:
        raise ValueError('reduced branch-cycle product is not identity')
    reached, queue = {0}, [0]
    for x in queue:
        for g in (gamma0, gamma1, gamma_inf):
            if g[x] not in reached:
                reached.add(g[x]); queue.append(g[x])
    if len(reached) != m:
        raise ValueError('reduced branch cycles are not transitive')
    cycle_types = [sorted(map(len, cycles(g))) for g in (gamma0, gamma1, gamma_inf)]
    twice_genus = 2 - 2*m + sum(m - len(c) for c in cycle_types)
    if twice_genus < 0 or twice_genus % 2:
        raise ValueError('invalid reduced curve genus')
    return dict(reduction_image_order=len(klein), reduced_states=m, blocks=blocks,
                branch_cycles=[gamma0, gamma1, gamma_inf], branch_cycle_types=cycle_types,
                reduced_curve_genus=twice_genus//2, presentation_verified=True,
                convention='Klein-four reduction of the supplied conjugacy-quotient action; curve genus, not cover genus')


def six_moves(generators):
    return list(generators) + [inverse(g) for g in generators]


def defect_lumping(generators, defects):
    """Coarsest strong Markov refinement respecting initial defect classes.

    This averages six equally probable moves. Move-labelled equivariance is a
    stronger, different property and is tested separately.
    """
    n, moves = len(defects), six_moves(generators)
    if any(len(g) != n for g in moves):
        raise ValueError('defect/action cardinality mismatch')
    def colors(values):
        labels = {value: i for i, value in enumerate(sorted(set(values)))}
        return [labels[value] for value in values]
    initial, current = colors(defects), colors(defects)
    def signatures(labels):
        width = max(labels)+1
        return [(labels[i], tuple(Counter(labels[g[i]] for g in moves).get(j,0) for j in range(width))) for i in range(n)]
    first = signatures(current)
    lumpable = all(len({first[i][1] for i in range(n) if current[i] == block}) == 1 for block in set(current))
    rounds = 0
    while True:
        new = colors(signatures(current))
        if len(set(new)) == len(set(current)):
            break
        current, rounds = new, rounds+1
        if rounds > n:
            raise AssertionError('partition refinement did not stabilize')
    blocks = [[i for i in range(n) if current[i] == j] for j in range(max(current)+1)]
    rows = [list(signatures(current)[block[0]][1]) for block in blocks]
    labelled = all(len({initial[g[i]] for i in range(n) if initial[i] == block}) == 1
                   for g in moves for block in set(initial))
    hist = Counter((defects[i], defects[g[i]]) for i in range(n) for g in moves)
    return dict(defect_strongly_lumpable=lumpable, defect_move_labelled_equivariant=labelled,
                refinement_rounds=rounds, refined_blocks=blocks, six_move_lumped_adjacency=rows,
                exact_defect_histogram=dict(sorted(Counter(defects).items())),
                directed_defect_transport=[dict(source=a,target=b,count=count) for (a,b),count in sorted(hist.items())],
                denominator='states for histograms; six directed move occurrences per state for transport')


def trace_moments(generators, maximum=12):
    n, moves = len(generators[0]), six_moves(generators)
    traces = [0]*(maximum+1)
    for root in range(n):
        vector = [0]*n; vector[root] = 1
        for k in range(maximum+1):
            traces[k] += vector[root]
            next_vector = [0]*n
            for x, count in enumerate(vector):
                if count:
                    for g in moves:
                        next_vector[g[x]] += count
            vector = next_vector
    radial, tree = {0:1}, []
    for _ in range(maximum+1):
        tree.append(radial.get(0,0))
        next_radial = defaultdict(int)
        for distance, count in radial.items():
            next_radial[distance+1] += count*(6 if distance == 0 else 5)
            if distance:
                next_radial[distance-1] += count
        radial = next_radial
    normalized = [str(Fraction(value,n*6**k)) for k,value in enumerate(traces)]
    return dict(adjacency_traces=traces, mean_random_walk_return_probabilities=normalized,
                infinite_six_regular_tree_moments=tree,
                arithmetic='exact integers and rational probabilities; loops and parallel moves included')


def relation_profile(generators):
    q1,q2,q3=generators; n=len(q1)
    relations = dict(braid12=[0,1,0,4,3,4], braid23=[1,2,1,5,4,5],
                     far_commute=[0,2,3,5], first_three_center=[0,1]*3,
                     four_strand_center=[0,1,2]*4)
    moves=six_moves(generators)
    result={}
    for name,word in relations.items():
        p=list(range(n))
        for mv in word:p=compose(moves[mv],p)
        result[name]=dict(word=word,fixed_points=sum(i==j for i,j in enumerate(p)),identity=p==list(range(n)))
    return result


def cycle_matched_controls(generators,replicates,seed,maximum=8):
    """Condition on generator cycle types, but do NOT condition on braid relations."""
    rng=random.Random(seed); n=len(generators[0]); results=[]
    for replicate in range(replicates):
        sample=[]
        for g in generators:
            q=list(range(n));rng.shuffle(q);sample.append(compose(compose(q,g),inverse(q)))
        results.append(dict(replicate=replicate,trace_moments=trace_moments(sample,maximum),
                            relations=relation_profile(sample)))
    return dict(seed=seed,replicates=results,
                population='independent uniform relabelings of each observed generator; not a braid-relation-matched null')


def sign(p):
    return (-1)**(len(p)-len(cycles(p)))


def coupled_fiber_hypothesis(generators, passports, braid_order):
    """Attempt an explicit parity-constrained S10^3 semidirect S4 identification.

    Group order alone is insufficient. Return an unsuccessful hypothesis test
    when a faithful coordinate embedding or parity character cannot be verified.
    """
    fibers=defaultdict(list)
    for i,p in enumerate(passports):fibers[p].append(i)
    keys=sorted(fibers);blocks=[fibers[k] for k in keys]
    if len(blocks)!=6 or {len(b) for b in blocks}!={10}:
        return dict(applicable=False,reason='requires six ordered-passport fibers of size ten')
    quotient=induced(generators,blocks)
    identity=list(range(len(passports)))
    qid=tuple(range(6));transports={qid:identity};queue=[qid]
    for p in queue:
        for g,q in zip(generators,quotient):
            dest=tuple(compose(q,p))
            if dest not in transports:
                if len(queue)>=24:raise ValueError('passport image larger than S4')
                transports[dest]=compose(g,transports[p]);queue.append(dest)
    if len(queue)!=24:return dict(applicable=True,identified=False,reason='passport image is not S4')
    loops=[]
    for p in queue:
        for g,q in zip(generators,quotient):
            dest=tuple(compose(q,p));loops.append(compose(inverse(transports[dest]),compose(g,transports[p])))
    actions=[]
    for block in blocks:
        lookup={x:i for i,x in enumerate(block)}
        actions.append([[lookup[g[x]] for x in block] for g in loops])
    edges=[];alignments={}
    for i in range(6):
        for j in range(i+1,6):
            mapping=equivariant_map(actions[i],actions[j])
            if mapping is not None:
                edges.append((i,j));alignments[i,j]=mapping
    if len(edges)!=3 or Counter(x for edge in edges for x in edge)!=Counter(range(6)):
        return dict(applicable=True,identified=False,kernel_equivalent_fiber_pairs=edges,
                    reason='kernel actions do not identify three unique complementary fiber pairs')
    owner={x:i for i,pair in enumerate(edges) for x in pair}
    aligned={}
    for i,j in edges:
        aligned[i]={x:k for k,x in enumerate(blocks[i])}
        mapping=alignments[i,j]
        aligned[j]={blocks[j][mapping[k]]:k for k in range(10)}
    generator_coordinates=[]
    for g,q in zip(generators,quotient):
        pair_perm=[];internal=[]
        for i,j in edges:
            if owner[q[i]]!=owner[q[j]]:
                raise ValueError('kernel-equivalence pairing not invariant')
            pair_perm.append(owner[q[i]])
            first=[aligned[q[i]][g[x]] for x in blocks[i]]
            # blocks[j] need sorting by its newly aligned coordinates.
            second=[None]*10
            for x in blocks[j]:second[aligned[j][x]]=aligned[q[j]][g[x]]
            if first!=second:
                return dict(applicable=True,identified=False,kernel_equivalent_fiber_pairs=edges,
                            reason='paired orientations have different internal generator actions in these coordinates')
            internal.append(first)
        parity=math.prod(sign(p) for p in internal)==sign(pair_perm)
        generator_coordinates.append(dict(passport_permutation=q,pair_permutation=pair_perm,
                                          internal_permutations=internal,parity_constraint=parity))
    expected=12*math.factorial(10)**3
    identified=all(r['parity_constraint'] for r in generator_coordinates) and int(braid_order)==expected
    return dict(applicable=True,identified=identified,kernel_equivalent_fiber_pairs=edges,
                ordered_passport_keys=keys,
                fiber_points_in_coordinate_order=[sorted(aligned[i],key=aligned[i].get) for i in range(6)],
                generator_coordinates=generator_coordinates,candidate_order=str(expected),
                faithful_coordinate_embedding_verified=True,
                identification=('parity-constrained subgroup of S10^3 semidirect S4 on three unoriented partitions and their two orientations'
                                if identified else 'coordinate embedding verified; proposed parity/order identification failed'))


def symmetry_common_quotient(generators, cap=65536):
    """Exact Goursat quotient for a relabeling action with two point blocks."""
    n=len(generators[0]);unseen=set(range(n));blocks=[]
    while unseen:
        root=min(unseen);block={root};queue=[root]
        for x in queue:
            for g in generators:
                if g[x] not in block:block.add(g[x]);queue.append(g[x])
        points=sorted(block);blocks.append(points);unseen.difference_update(points)
    if len(blocks)!=2:
        return dict(applicable=False,point_orbits=list(map(len,blocks)))
    elements=generated_elements(generators,cap)
    if elements is None:return dict(applicable=True,status='capped',element_cap=cap)
    def restrict(p,block):
        labels={x:i for i,x in enumerate(block)}
        return tuple(labels[p[x]] for x in block)
    pairs=[tuple(restrict(p,b) for b in blocks) for p in elements]
    left,right=({p[i] for p in pairs} for i in range(2))
    lid,rid=tuple(range(len(blocks[0]))),tuple(range(len(blocks[1])))
    lk={a for a,b in pairs if b==rid};rk={b for a,b in pairs if a==lid}
    def coset(a):return tuple(min(compose(a,k) for k in lk))
    representatives=sorted({coset(a) for a in left});lookup={p:i for i,p in enumerate(representatives)}
    table=[[lookup[coset(compose(a,b))] for b in representatives] for a in representatives]
    abelian=all(table[i][j]==table[j][i] for i in range(len(table)) for j in range(len(table)))
    quotient_order=len(representatives)
    if quotient_order*len(elements)!=len(left)*len(right):raise ValueError('Goursat order identity failed')
    quotient_identity=lookup[coset(lid)]
    generator_images=[lookup[coset(restrict(g,blocks[0]))] for g in generators]
    generator_orders=[]
    for image in generator_images:
        value=quotient_identity
        for order in range(1,quotient_order+1):
            value=table[value][image]
            if value==quotient_identity:
                generator_orders.append(order);break
        else:raise ValueError('invalid quotient generator order')
    return dict(applicable=True,status='computed',point_orbits=list(map(len,blocks)),image_order=str(len(elements)),
                point_blocks=blocks,
                block_image_orders=[len(left),len(right)],block_common_quotient_kernel_orders=[len(lk),len(rk)],
                common_quotient_order=quotient_order,common_quotient_abelian=abelian,
                common_quotient=('C6' if abelian else 'S3') if quotient_order==6 else 'not named',
                quotient_multiplication_table=table,quotient_identity=quotient_identity,
                generator_quotient_images=generator_images,generator_quotient_orders=generator_orders,
                arithmetic='explicit finite group and exact cosets; no order-only identification')


def analyze(generators,passports,defects,braid_order,replicates=0,seed=20261004):
    return dict(reduced_hurwitz=reduced_hurwitz(generators),
                defect_transport=defect_lumping(generators,defects),
                trace_moments=trace_moments(generators),relations=relation_profile(generators),
                coupled_fibers=coupled_fiber_hypothesis(generators,passports,braid_order),
                cycle_matched_controls=cycle_matched_controls(generators,replicates,seed))


def selftest():
    states=list(itertools.permutations(range(4)));lookup={p:i for i,p in enumerate(states)}
    generators=[]
    for j in range(3):
        row=[]
        for p in states:
            q=list(p);q[j],q[j+1]=q[j+1],q[j];row.append(lookup[tuple(q)])
        generators.append(row)
    reduced=reduced_hurwitz(generators)
    assert reduced['reduced_states']==6 and reduced['reduced_curve_genus']==0
    assert all(r['identity'] for r in relation_profile(generators).values())
    trivial=defect_lumping(generators,[0]*24)
    assert trivial['defect_strongly_lumpable'] and trivial['six_move_lumped_adjacency']==[[6]]
    refined=defect_lumping(generators,[int(i==0) for i in range(24)])
    assert not refined['defect_strongly_lumpable'] and len(refined['refined_blocks'])>2
    traces=trace_moments(generators,6)['adjacency_traces']
    assert traces[0]==24 and traces[1]==0 and traces[2]==288
    assert cycle_matched_controls(generators,2,1)==cycle_matched_controls(generators,2,1)
    diagonal=[[1,0,2,4,3,5],[1,2,0,4,5,3]]
    common=symmetry_common_quotient(diagonal)
    assert common['image_order']=='6' and common['common_quotient']=='S3' and common['block_common_quotient_kernel_orders']==[1,1]
    print('saved weak-graph extensions passed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state-action',type=Path);parser.add_argument('--out',type=Path)
    parser.add_argument('--selftest',action='store_true');args=parser.parse_args()
    if args.selftest or args.state_action is None:
        if args.out is not None:parser.error('state-action required with out')
        selftest();return
    if args.out is None:parser.error('out required with state-action')
    text=args.state_action.read_text();lines=text.splitlines()
    if not lines or lines[0]!='HURWITZ_BRAID_STATE_ACTION_V1':raise ValueError('expected a STATE action, not point relabelings')
    fields=dict(line.split('=',1) for line in lines[1:]);n=int(fields['degree'])
    generators=[list(map(int,fields[f'sigma{i}'].split(','))) for i in range(1,4)]
    if any(sorted(g)!=list(range(n)) for g in generators):raise ValueError('invalid state-action permutation')
    if not all(r['identity'] for r in relation_profile(generators).values()):raise ValueError('invalid abstract braid relations')
    document=dict(state_count=n,reduced_hurwitz=reduced_hurwitz(generators),trace_moments=trace_moments(generators),
                  no_generator_fixed_anywhere=all(all(i!=j for i,j in enumerate(g)) for g in generators),
                  source_sha256=hashlib.sha256(args.state_action.read_bytes()).hexdigest(),
                  scope='abstract finite braid-state action; no tuple realization, covering genus or weak failure asserted')
    with args.out.open('x') as stream:json.dump(document,stream,indent=2);stream.write('\n')


if __name__=='__main__':main()
