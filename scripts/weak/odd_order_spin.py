#!/usr/bin/env python3
"""Bounded exact odd-order lift invariant in the Clifford model of 2.A_n.

This is an interpretation experiment, not a weak-counterexample test. Even-
order entries have no canonical odd-order lift and are explicitly inapplicable.
Arbitrary-precision integers avoid numerical sign ambiguity. Degree and work
caps return 'capped', never an invented invariant.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import random
import re


class Capped(Exception):
    pass


def compose(a,b):
    return tuple(a[i] for i in b)


def inverse(p):
    result=[0]*len(p)
    for i,j in enumerate(p):result[j]=i
    return tuple(result)


def cycle_lengths(p):
    seen=set();result=[]
    for root in range(len(p)):
        if root in seen:continue
        x,length=root,0
        while x not in seen:seen.add(x);length+=1;x=p[x]
        result.append(length)
    return result


def normalize(coefficients,exponent):
    coefficients={mask:c for mask,c in coefficients.items() if c}
    while exponent and all(c%2==0 for c in coefficients.values()):
        coefficients={mask:c//2 for mask,c in coefficients.items()};exponent-=1
    return coefficients,exponent


def blade_sign(a,b):
    # e_I e_J: one minus sign for every reversed distinct-index pair;
    # e_i^2=+1 and coincident factors disappear from the blade mask.
    parity=0
    while a:
        bit=a&-a;a-=bit;parity^=(b&(bit-1)).bit_count()&1
    return -1 if parity else 1


def multiply(left,right):
    a,ae=left;b,be=right
    if len(a)*len(b)>2_000_000:raise Capped('Clifford multiplication exceeds the operation cap')
    result={}
    for i,x in a.items():
        for j,y in b.items():
            mask=i^j;result[mask]=result.get(mask,0)+blade_sign(i,j)*x*y
    if len(result)>4096:raise Capped('Clifford support exceeds 4096 blades')
    return normalize(result,ae+be)


def power(lift,n):
    result=({0:1},0)
    while n:
        if n&1:result=multiply(result,lift)
        n>>=1
        if n:lift=multiply(lift,lift)
    return result


@lru_cache(maxsize=512)
def canonical_lift(p):
    """Bubble factorization into adjacent transpositions, normalized exactly.

    Each adjacent lift is (e_i-e_{i+1})/sqrt(2). Only even permutations
    enter this interface, so the accumulated denominator is an integer power
    of two. No choice of a lift is inferred for an even-order conjugacy class.
    """
    current=list(range(len(p)));word=[]
    for i,value in enumerate(p):
        j=current.index(value,i)
        while j>i:
            current[j-1],current[j]=current[j],current[j-1];word.append(j-1);j-=1
    if len(word)%2:raise ValueError('spin invariant requires alternating permutations')
    coefficients={0:1}
    for index in word:
        result={}
        for blade,c in coefficients.items():
            for point,weight in ((index,1),(index+1,-1)):
                bit=1<<point;mask=blade^bit
                sign=-1 if (blade>>(point+1)).bit_count()%2 else 1
                result[mask]=result.get(mask,0)+sign*weight*c
        coefficients={blade:c for blade,c in result.items() if c}
        if len(coefficients)>4096:raise Capped('canonical lift support cap')
    return normalize(coefficients,len(word)//2)


def central_sign(lift):
    coefficients,exponent=normalize(*lift)
    if exponent or coefficients not in ({0:1},{0:-1}):
        raise ValueError('expected a central plus/minus one lift')
    return coefficients[0]


def odd_lift(p):
    order=math.lcm(*cycle_lengths(p))
    if order%2==0:raise ValueError('even-order entry: canonical odd-order spin lift is inapplicable')
    lift=canonical_lift(p);value=central_sign(power(lift,order))
    return ({mask:value*c for mask,c in lift[0].items()},lift[1])


def invariant(t):
    d=len(t[0]);identity=tuple(range(d))
    if not 4<=d<=12:raise Capped('exact exploratory spin degree is limited to 4..12')
    if len(t)!=4 or any(len(p)!=d or sorted(p)!=list(identity) for p in t):raise ValueError('invalid tuple permutations')
    product=identity
    for p in t:product=compose(product,p)
    if product!=identity:raise ValueError('tuple is not product one')
    lifted=({0:1},0)
    for p in t:lifted=multiply(lifted,odd_lift(tuple(p)))
    return central_sign(lifted)


def hurwitz(t,mv):
    t=list(t);i=mv%3;a,b=t[i:i+2]
    if mv<3:t[i:i+2]=[compose(compose(a,b),inverse(a)),a]
    else:t[i:i+2]=[b,compose(compose(inverse(b),a),b)]
    return tuple(t)


def parse_seed(raw,d):
    parts=raw.strip()[1:-1].split(', ')
    if len(parts)!=4:raise ValueError('certificate seed does not have four entries')
    result=[]
    for part in parts:
        p=list(range(d));used=set()
        for cycle in re.findall(r'\(([^()]*)\)',part):
            points=[int(x)-1 for x in cycle.split()]
            if any(x<0 or x>=d or x in used for x in points) or len(set(points))!=len(points):raise ValueError('invalid certificate cycle')
            used.update(points)
            for a,b in zip(points,points[1:]+points[:1]):p[a]=b
        result.append(tuple(p))
    return tuple(result)


def selftest():
    # The three-cycle lift chosen by odd order has cube +1.
    p=(1,2,0,3)
    assert central_sign(power(odd_lift(p),3))==1
    # Disjoint transpositions lift to an element whose square is -1.
    assert central_sign(power(canonical_lift((1,0,3,2)),2))==-1
    rng=random.Random(20261004);checked=0
    for d in range(4,8):
        for _ in range(256):
            entries=[]
            for j in range(3):
                q=list(range(d));rng.shuffle(q);entries.append(tuple(q))
            entries.append(inverse(compose(compose(entries[0],entries[1]),entries[2])))
            if any(any(n%2==0 for n in cycle_lengths(p)) for p in entries):continue
            value=invariant(tuple(entries))
            for mv in range(6):assert invariant(hurwitz(tuple(entries),mv))==value
            checked+=1
    assert checked>0
    try:odd_lift((1,0,3,2))
    except ValueError:pass
    else:raise AssertionError('even-order lift accepted')
    print('bounded exact odd-order spin checks passed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path);parser.add_argument('--out',type=Path)
    parser.add_argument('--selftest',action='store_true');args=parser.parse_args()
    if args.selftest:selftest();return
    if args.certificate is None or args.out is None:parser.error('certificate and out required')
    fields=dict(line.split('=',1) for line in args.certificate.read_text().splitlines() if '=' in line)
    d=int(fields['d']);t=parse_seed(fields['seed'],d)
    document=dict(source=str(args.certificate),source_sha256=hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
                  implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),degree=d,
                  scope='exact bounded odd-order lift interpretation; not a weak-failure or group-generation verification')
    if any(any(n%2==0 for n in cycle_lengths(p)) for p in t):
        document.update(status='inapplicable_even_order',reason='at least one entry has even order')
    else:
        try:
            value=invariant(t);moves=[invariant(hurwitz(t,mv)) for mv in range(6)]
            if any(x!=value for x in moves):raise ValueError('spin lift failed literal braid invariance')
            document.update(status='computed',spin=value,six_braid_moves_verified=True)
        except Capped as error:document.update(status='capped',reason=str(error))
    with args.out.open('x') as stream:json.dump(document,stream,indent=2);stream.write('\n')


if __name__=='__main__':main()
