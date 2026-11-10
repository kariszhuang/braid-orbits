#!/usr/bin/env python3
"""Literal positive-word audit of closed exported proof tables; stdlib only.

The inputs select difficult leads or bounded successful probe walks. Recorded
seed strings are deduplicated, not counted as independent braid components.
No prevalence, shortest-word, full-group or nonexistence theorem is inferred.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path

from odd_order_spin import compose, cycle_lengths, hurwitz, parse_seed


def transitive(perms):
    found={0};queue=[0]
    for point in queue:
        for p in perms:
            if p[point] not in found:found.add(p[point]);queue.append(p[point])
    return len(found)==len(perms[0])


def free_reduce(word):
    result=[]
    for move in word:
        if result and (result[-1]+3)%6==move:result.pop()
        else:result.append(move)
    return tuple(result)


def pair_partition(perms):
    unseen=set(range(len(perms[0])));sizes=[]
    while unseen:
        root=min(unseen);found={root};queue=[root]
        for point in queue:
            for p in perms:
                if p[point] not in found:found.add(p[point]);queue.append(p[point])
        unseen-=found;sizes.append(len(found))
    return tuple(sorted(sizes,reverse=True))


def decoration(t):
    d=len(t[0]);left=pair_partition(t[:2]);right=pair_partition(t[2:])
    indices=[d-len(cycle_lengths(p)) for p in t]
    product_cycles=len(cycle_lengths(compose(t[0],t[1])))
    delta=product_cycles-len(left)-len(right)+1
    gl=len(left)-d+(indices[0]+indices[1]+d-product_cycles)//2
    gr=len(right)-d+(indices[2]+indices[3]+d-product_cycles)//2
    genus=1-d+sum(indices)//2
    if min(delta,gl,gr)<0 or delta+gl+gr!=genus:raise ValueError('pair-genus/defect identity failed')
    return left,right,delta,gl,gr


def checked_word(d,raw,moves):
    if not 4<=d<=34 or len(moves)>4096 or any(m not in range(6) for m in moves):
        raise ValueError('invalid degree or braid word')
    t=parse_seed(raw,d);identity=tuple(range(d));product=identity
    for p in t:product=compose(product,p)
    if product!=identity or not transitive(t):raise ValueError('seed must be product one and point transitive')
    word=free_reduce(moves);u=t
    for move in word:u=hurwitz(u,move)
    if not transitive(u[:2]):raise ValueError('false transitive-first-pair witness')
    return t,word


def analyze(paths):
    sources={};seeds={};occurrences=0
    for path in paths:
        path=Path(path)
        if path.name=='POSITIVE_WITNESSES.tsv' and not (path.parent/'VERIFIED').exists():
            raise ValueError('attempt output needs VERIFIED before analysis')
        payload=path.read_bytes();sources[str(path)]=hashlib.sha256(payload).hexdigest()
        reader=csv.DictReader(io.StringIO(payload.decode()),delimiter='\t')
        if reader.fieldnames!=['degree','moves','seed']:raise ValueError('invalid proof table header')
        for row in reader:
            d=int(row['degree']);moves=tuple(map(int,row['moves'].split(','))) if row['moves'] else ()
            t,word=checked_word(d,row['seed'],moves);occurrences+=1
            key=(d,row['seed'])
            if key not in seeds or len(word)<len(seeds[key][1]):seeds[key]=(t,word)
    words=Counter();lengths=Counter();passports=Counter();groups=Counter();structures=Counter()
    for (d,_),(t,word) in seeds.items():
        indices=[d-len(cycle_lengths(p)) for p in t]
        group='S' if any(i%2 for i in indices) else 'A'
        genus=1-d+sum(indices)//2
        passport='/'.join(sorted('.'.join(map(str,sorted(cycle_lengths(p)))) for p in t))
        words[(d,group,word)]+=1;lengths[(d,group,len(word))]+=1
        passports[(d,group,genus,passport,word)]+=1;groups[(d,group)]+=1
        end=t
        for move in word:end=hurwitz(end,move)
        structures[(d,group,genus,decoration(t),decoration(end),word)]+=1
    return dict(source_sha256=sources,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        permutation_helper_sha256=hashlib.sha256(Path(__file__).with_name('odd_order_spin.py').read_bytes()).hexdigest(),
        literal_rows_checked=occurrences,distinct_recorded_seeds=len(seeds),
        groups=[dict(degree=d,group=g,seeds=n) for (d,g),n in sorted(groups.items())],
        word_lengths=[dict(degree=d,group=g,length=l,seeds=n) for (d,g,l),n in sorted(lengths.items())],
        recurring_words=[dict(degree=d,group=g,moves=list(w),seeds=n) for (d,g,w),n in words.most_common()],
        passport_words=[dict(degree=d,group=g,genus=genus,passport=p,moves=list(w),seeds=n)
                        for (d,g,genus,p,w),n in passports.most_common()],
        escape_structures=[dict(degree=d,parity_ambient=g,genus=genus,
            start=dict(zip(('left_orbits','right_orbits','defect','left_genus','right_genus'),before)),
            end=dict(zip(('left_orbits','right_orbits','defect','left_genus','right_genus'),after)),moves=list(w),seeds=n)
            for (d,g,genus,before,after,w),n in structures.most_common()],
        scope='Literal replay on selected exported proofs; recorded seed-string deduplication, not full conjugacy/component deduplication; group labels are parity-compatible ambient targets, not verified full generation; early construction proofs may generate proper groups; words need not be shortest; no uniform frequencies or universal lemma')


def selftest():
    seed='[(2 3)(5 7)(6 8)(9 10), (3 4)(5 6)(7 9)(8 10), (1 6 10)(2 9)(3 8)(4 5), (1 10 7)(2 5)(3 9)(4 8)]'
    _,word=checked_word(10,seed,(0,3,1,3,1));assert word==(1,3,1)
    for moves in [(),(6,)]:
        try:checked_word(10,seed,moves)
        except ValueError:pass
        else:raise AssertionError('false witness accepted')
    assert free_reduce((0,1,4,3))==()
    t,word=checked_word(10,seed,(1,3,1));before=decoration(t)
    end=t
    for move in word:end=hurwitz(end,move)
    assert pair_partition(end[:2])==(10,)
    assert sum(before[2:])==sum(decoration(end)[2:])
    print('literal escape-word checks passed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table',type=Path,action='append',default=[])
    parser.add_argument('--out',type=Path);parser.add_argument('--selftest',action='store_true')
    args=parser.parse_args()
    if args.selftest:selftest();return
    if not args.table or args.out is None:parser.error('table(s) and out required')
    result=analyze(args.table)
    with args.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:result[k] for k in ('literal_rows_checked','distinct_recorded_seeds','groups')}))


if __name__=='__main__':main()
