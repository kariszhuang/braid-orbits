#!/usr/bin/env python3
"""Independent S5/S6 marked parameter coverage; no Rust/engine imports.

This is the specified marked full/partial-fibre family, NOT an entire degree. Ordinary
negative words are audited separately by smooth_negative_audit. Here we verify
that every qualifying parameter was actually submitted with its literal tuple.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

from verify_smooth_geometry import parse_tuple


def compose(a, b):
    return tuple(a[x] for x in b)


def inverse(a):
    b = [0] * len(a)
    for x, y in enumerate(a):
        b[y] = x
    return tuple(b)


def index(a):
    unseen = set(range(len(a))); cycles = 0
    while unseen:
        cycles += 1; x = min(unseen)
        while x in unseen:
            unseen.remove(x); x = a[x]
    return len(a) - cycles


def audit(folder):
    folder = Path(folder)
    coverage = json.loads((folder / 'FIBER_COVERAGE.json').read_text())
    kind = coverage['kind']; r = 3 if kind in ['C3','C3P31','C3P22'] else 4
    k=coverage['inner_degree']
    if kind not in ['C3', 'C3P31', 'C3P22', 'C4', 'S4'] or k not in [5,6]:
        raise ValueError('only the independently specified S5/S6 families are supported')
    n = r * k
    permutations = list(itertools.permutations(range(k))); size=len(permutations)
    ranks = {p: i for i, p in enumerate(permutations)}
    mul = [[ranks[compose(a, b)] for b in permutations] for a in permutations]
    inv = [ranks[inverse(a)] for a in permutations]
    conj = [[mul[mul[h][p]][inv[h]] for h in range(size)] for p in range(size)]
    trans = [i for i, p in enumerate(permutations) if index(p) == 1]
    cycle3 = ranks[(1, 2, 0)+tuple(range(3,k))]; double = ranks[(1, 0, 3, 2)+tuple(range(4,k))]
    swap = ranks[(1, 0)+tuple(range(2,k))]
    if kind in ['C3P31','C3P22']:
        first_index,second_index=(3,1) if kind=='C3P31' else (2,2)
        firsts={}
        for i,p in enumerate(permutations):
            if index(p)==first_index:
                unseen=set(range(k));shape=[]
                while unseen:
                    start=min(unseen);x=start;length=0
                    while x in unseen:unseen.remove(x);length+=1;x=p[x]
                    shape.append(length)
                firsts.setdefault(tuple(sorted(shape)),i)
        seconds=[i for i,p in enumerate(permutations) if index(p)==second_index]
        locals_raw=((a,b,c) for a in firsts.values() for other in seconds
                    for b,c in [(other,0),(0,other)])
    else:
        locals_raw = (itertools.product([cycle3, double], trans, trans) if r == 3
                      else itertools.product([swap], trans, trans, trans))
    representatives = {}
    for local in locals_raw:
        # Return of A*beta on block zero, evaluated rightmost first.
        product = local[1]
        for j in range(2, r):
            product = mul[local[j]][product]
        product = mul[local[0]][product]
        if index(permutations[product]) != 4:
            continue
        key = min(tuple(y for p in local for y in permutations[conj[p][h]])
                  for h in range(size))
        representatives.setdefault(key, local)
    frames = [representatives[k] for k in sorted(representatives)]
    if coverage['local_frames'] != len(frames):
        raise ValueError('local-frame normalization mismatch')
    total = len(frames) * size**2
    if any(coverage[k] != v for k, v in [('degree', n), ('start', 0),
            ('next', total), ('total', total), ('internal_permutations', size)]):
        raise ValueError('incomplete or inconsistent parameter coverage')
    if coverage['producer_complete'] is not True:
        raise ValueError('producer did not finish')
    observed = {}
    with (folder / 'LEADS.tsv').open() as stream:
        for row in csv.reader(stream, delimiter='\t'):
            rank = int(row[0])
            if rank in observed:
                raise ValueError('duplicate literal candidate rank')
            observed[rank] = parse_tuple(row[1], n)
    expected = set(); outside_transitive=0
    for fi, a in enumerate(frames):
        for u in range(size):
            if kind in ['C3','C3P31','C3P22']:
                x = mul[a[2]][inv[u]]; y = inv[mul[mul[a[0]][u]][a[1]]]
            elif kind == 'C4':
                x = mul[a[3]][inv[u]]; y = inv[mul[mul[mul[a[0]][u]][a[2]]][a[1]]]
            else:
                x = mul[a[2]][inv[u]]; y = inv[mul[mul[mul[a[0]][a[3]]][u]][a[1]]]
            for v in range(size):
                if conj[x][v] != y:
                    continue
                rank = fi * size**2 + u * size + v + 1
                if kind == 'C4':
                    v1 = mul[mul[inv[a[3]]][inv[v]]][inv[a[1]]]
                    voltage = [mul[mul[inv[v1]][inv[u]]][inv[v]], v1, u, v]
                else:
                    voltage = [mul[inv[u]][inv[v]], u, v] + ([0] if r == 4 else [])
                gamma = list(range(1, r)) + [0]
                if kind == 'S4':
                    gamma = [1, 2, 0, 3]
                aa = tuple(j*k + permutations[a[j]][i] for j in range(r) for i in range(k))
                bb = tuple(((j+1) % r)*k+i for j in range(r) for i in range(k))
                cc = tuple(gamma[j]*k+permutations[voltage[j]][i] for j in range(r) for i in range(k))
                dd = inverse(compose(compose(aa, bb), cc))
                seen={0};todo=[0]
                while todo:
                    point=todo.pop()
                    for perm in [aa,bb,cc]:
                        if perm[point] not in seen:seen.add(perm[point]);todo.append(perm[point])
                if len(seen)!=n:
                    outside_transitive+=1;continue
                expected.add(rank)
                if observed.get(rank) != [aa, bb, cc, dd]:
                    # parse_tuple returns a tuple or list depending on caller;
                    # compare immutable permutations, never displayed notation.
                    if tuple(observed.get(rank, ())) != (aa, bb, cc, dd):
                        raise ValueError(f'literal parameter mismatch at {rank}')
    if set(observed) != expected:
        raise ValueError('missing or extra qualifying parameters')
    with (folder / 'COUNTS.tsv').open() as stream:
        counts = next(csv.DictReader(stream, delimiter='\t'))
    if int(counts['constructed']) != total or int(counts['scope']) != len(expected):
        raise ValueError('construction counters do not reconcile')
    if coverage.get('outside_transitive',0)!=outside_transitive:
        raise ValueError('disconnected parameter counter mismatch')
    native = json.loads((folder / 'CHECKS.json').read_text())['independent_native']
    if native['records'] != len(expected) or native['unresolved']:
        raise ValueError('qualifying parameters lack completed independent resolutions')
    return dict(passed=True, kind=kind, degree=n, inner_degree=k,local_frames=len(frames),parameter_count=total,
        qualifying_parameters=len(expected), nonsmooth=native['nonsmooth'],
        closed=native['closed'], unresolved=0,outside_transitive=outside_transitive,
        literal_journal_sha256=hashlib.sha256((folder/'LEADS.tsv').read_bytes()).hexdigest(),
        scope=f'complete normalized marked {kind} S{k} parameter family; NOT a degree census')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.folder), indent=2))
