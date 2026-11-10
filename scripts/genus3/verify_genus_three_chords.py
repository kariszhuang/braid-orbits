#!/usr/bin/env python3
"""Independent finite coverage reconstruction and literal negative-word replay.

The producer grows vertex colours; this verifier partitions subsets of the
point set, then independently orders each cycle. No producer modules imported.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
from itertools import combinations, permutations, product
import json
from math import comb, gcd
from pathlib import Path
import time

from audit_guided_native import audit_guided_native
from run_genus_three import sha, write
from verify_genus_three_finite import mul, inv, profile, braid
from verify_smooth_geometry import parse_tuple

if not __debug__:
    raise RuntimeError('Assertions are audit gates; do not use Python -O')


def is_smooth(t):
    n=len(t[0]); ix=[n-len(profile(p)) for p in t]
    return n-len(profile(mul(t[0],t[1]))) in (sum(ix[:2]),sum(ix[2:]))


def frames(n, mode):
    identity=tuple(range(n)); p=tuple((x+1)%n for x in range(n))
    if mode=='transposition':
        out=[]
        for r in range(1,n//2+1):
            a=list(identity);a[0],a[r]=r,0;out.append(tuple(a))
        return out
    # Deliberately enumerate ALL matchings, including crossing chords, and
    # select geodesics by literal cycle count, unlike the producer's two forms.
    out=set(); raw=0
    for vertices in combinations(range(n),4):
        i=vertices[0]
        for j in vertices[1:]:
            k,l=[x for x in vertices if x not in (i,j)]
            a=list(identity);a[i],a[j]=j,i;a[k],a[l]=l,k;a=tuple(a)
            if len(profile(mul(a,p)))!=3:continue
            raw+=1
            rotated=[]
            for shift in range(n):
                u=[0]*n
                for x in range(n):u[(x+shift)%n]=(a[x]+shift)%n
                rotated.append(tuple(u))
            assert len(set(rotated))==n
            out.add(min(rotated))
    assert raw==2*comb(n,4) and len(out)*n==raw
    return sorted(out)


def admissible_partitions(a):
    """Unlabelled THREE-block partitions, checked against every full twist."""
    n=len(a); full=(1<<n)-1
    edges=[(i,j) for i,j in enumerate(a) if i<j]
    shifts=[[( (i+k)%n,(j+k)%n) for i,j in edges] for k in range(n)]
    independent=[]
    for bits in range(1,full+1):
        if all(not(bits>>i&1 and bits>>j&1) for pairs in shifts for i,j in pairs):independent.append(bits)
    allowed=set(independent)
    for first in independent:
        if not first&1:continue
        rest=full^first
        if not rest:continue
        root=rest&-rest
        for second in independent:
            if second&rest!=second or not second&root:continue
            third=rest^second
            if third not in allowed:continue
            blocks=[tuple(x for x in range(n) if bits>>x&1) for bits in (first,second,third)]
            if len(edges)==2:
                labels={x:k for k,block in enumerate(blocks) for x in block}
                if any(len({labels[x] for pair in pairs for x in pair})!=3 for pairs in shifts):continue
            yield blocks


def cycle_permutations(blocks,n):
    for tails in product(*(list(permutations(block[1:])) for block in blocks)):
        p=list(range(n))
        for block,tail in zip(blocks,tails):
            cycle=(block[0],)+tail
            for x,y in zip(cycle,cycle[1:]+cycle[:1]):p[x]=y
        yield tuple(p)


def audit(folder,auditor):
    start=time.perf_counter();folder=Path(folder);coverage=json.loads((folder/'COVERAGE.json').read_text())
    n=coverage['degree'];mode=coverage['mode'];assert n in (7,9,11) and mode in ('transposition','double')
    rank=1 if mode=='transposition' else 2
    assert coverage['ordered_indices']==[rank,n-rank-1,n-3,n-3]
    assert coverage['genus']==(n-5)//2
    p=tuple((x+1)%n for x in range(n));f=frames(n,mode)
    assert [tuple(row['special']) for row in coverage['frames']]==f
    counts=Counter(); saved={}; words=Counter(); max_word=0
    for name in ('RESOLUTIONS.tsv','SURVIVORS.tsv'):
        with (folder/name).open()as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                t=tuple(tuple(q) for q in parse_tuple(row['seed'],n));a,e,b,c=t
                assert all(sorted(q)==list(range(n)) for q in t)
                assert mul(a,e)==p and mul(p,mul(b,c))==tuple(range(n))
                key=(a,b);assert key not in saved,'duplicate normalized state'
                assert is_smooth(t)
                status='pending' if name=='SURVIVORS.tsv' else 'negative'
                saved[key]=status;counts[status]+=1
                if status=='negative':
                    assert row['status']=='nonsmooth' and not row['folder']
                    word=tuple(map(int,row['word'].split(',')));assert word and all(0<=m<6 for m in word)
                    u=t
                    for m in word:u=braid(u,m)
                    assert not is_smooth(u),'false negative word'
                    max_word=max(max_word,len(word));words[len(word)]+=1
    native=audit_guided_native(folder,auditor,native_components=True)
    assert native['independent_native']['records']==counts['negative']
    reports=[];seen=set()
    for number,a in enumerate(f):
        partitions=constructed=qualified=negative=pending=0
        for blocks in admissible_partitions(a):
            partitions+=1
            for b in cycle_permutations(blocks,n):
                constructed+=1
                c=inv(mul(p,b))
                if len(profile(c))!=3:continue
                qualified+=1;key=(a,b);assert key not in seen;seen.add(key)
                status=saved.pop(key);negative+=status=='negative';pending+=status=='pending'
        reported=coverage['frames'][number]
        expected=dict(frame=number,special=list(a),partitions=partitions,constructed=constructed,qualified=qualified,negative=negative,pending=pending)
        assert reported==expected,('frame coverage',number,reported,expected)
        if mode=='transposition':
            r=a[0];g=gcd(n,r);length=n//g
            assert partitions==((2**length-2)**g)//6
        reports.append(expected)
    assert not saved,'extra/out-of-scope normalized states'
    for k in ('constructed','qualified','negative','pending'):assert coverage[k]==sum(r[k] for r in reports)
    assert len(seen)==coverage['qualified']==counts['negative']+counts['pending']
    return dict(passed=True,status='complete_negative_slice' if counts['pending']==0 else 'pending_local_survivors',coverage=coverage,
        maximum_word=max_word,word_length_histogram=dict(words),independent_seconds=time.perf_counter()-start,
        native=native,files={q.name:sha(q)for q in folder.iterdir()if q.is_file()and q.name not in ('EXACT_AUDIT.json',)},
        scope='complete full-product slice after necessary twist colouring; counts are S_n states in the selected layout,not braid components')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('auditor',type=Path);a=p.parse_args()
    result=audit(a.folder,a.auditor);write(a.folder/'EXACT_AUDIT.json',result)
    print(json.dumps({k:result[k]for k in ('status','coverage','independent_seconds')}))


if __name__=='__main__':main()
