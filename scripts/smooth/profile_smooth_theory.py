"""Independent block/quotient invariants for certified smooth components.

No constructor imports. Classify one verified invariant block system, not all
factorizations, moduli loci or isotriviality. A proper invariant partition is a
proof of imprimitivity; equal orders are not identities of marked actions.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from collections import Counter
import json
from pathlib import Path

from sympy.combinatorics import Permutation,PermutationGroup
from verify_smooth_geometry import (parse_tuple,cycles,compose,forced_partition,orbits)
from run_smooth_degree9 import sha


def recover_general(t,support,r,k):
    base=forced_partition(t,[sorted(support)]);stack=[base];seen=set();anchor=min(support)
    while stack:
        p=stack.pop()
        if p in seen:continue
        seen.add(p)
        if len(seen)>10000:raise ValueError('block recovery unresolved at cap')
        block=next(x for x in p if anchor in x)
        if len(block)>k:continue
        if len(block)==k:
            if len(p)!=r or any(len(b)!=k for b in p):raise ValueError('invariant partition size')
            return p
        for other in p:
            if other!=block and len(block)+len(other)<=k:
                stack.append(forced_partition(t,[list(block+other)]))
    raise ValueError('no required unramified block system')


def profile(seed,n,r=None,k=None,compute_orders=True):
    t=parse_tuple(seed,n);supports=[{x for x in range(n) if p[x]!=x} for p in t]
    if len(orbits(t))!=1 or compose(compose(compose(t[0],t[1]),t[2]),t[3])!=tuple(range(n)):
        raise ValueError('profile product/transitivity')
    shapes=[(r,k)] if r is not None else [(r,n//r) for r in range(3,n//3+1) if n%r==0]
    errors=[]
    for r,k in shapes:
        local=[i for i,s in enumerate(supports) if 0<len(s)<=k]
        if len(local)!=1:continue
        s=local[0]
        try:blocks=recover_general(t,supports[s],r,k)
        except ValueError as error:errors.append(str(error));continue
        membership={x:j for j,b in enumerate(blocks) for x in b}
        q=[];valid=True
        for i,p in enumerate(t):
            image=[{membership[p[x]] for x in b} for b in blocks]
            if any(len(v)!=1 for v in image):valid=False;break
            qi=tuple(next(iter(v)) for v in image);q.append(qi)
            if i==s:
                valid&=qi==tuple(range(r)) and len({membership[x] for x in supports[i]})==1
            else:
                expected=sorted(len(c) for c in cycles(qi) for _ in range(k))
                valid&=expected==sorted(map(len,cycles(p)))
        if not valid:continue
        htotal=sum(r-len(cycles(p)) for p in q);total=sum(n-len(cycles(p)) for p in t)
        if htotal%2 or total%2:raise ValueError('genus parity')
        h=1-r+htotal//2;g=1-n+total//2;index=n-len(cycles(t[s]))
        if g!=1+k*(h-1)+index//2 or index%2 or h<1:raise ValueError('relative genus/parity')
        quotient_order=int(PermutationGroup([Permutation(list(p)) for p in q]).order()) if compute_orders else None
        order=int(PermutationGroup([Permutation(list(p)) for p in t]).order()) if compute_orders else None
        large_types=[sorted(map(len,cycles(q[i])),reverse=True) for i in range(4) if i!=s]
        local_block=membership[min(supports[s])]
        return dict(degree=n,blocks=r,block_size=k,block_partition=[[x+1 for x in b] for b in blocks],
            quotient_genus=h,genus=g,local_index=index,
            local_cycle_type=sorted([len(c) for c in cycles(t[s]) if membership[c[0]]==local_block],reverse=True),
            quotient_cycle_types=large_types,quotient_order=quotient_order,group_order=order,
            equal_three_quotient_types=len({tuple(x) for x in large_types})==1,
            all_three_full_quotient_cycles=all(x==[r] for x in large_types),
            imprimitive=True,relative_genus_checked=True,representative=seed,
            geometric_classification='not classified here; theorem review separate')
    raise ValueError('missing single-local/unramified structure: '+','.join(errors))


def structural_summary(entries):
    signatures=Counter()
    for e in entries:
        key=json.dumps([e['blocks'],e['block_size'],e['quotient_genus'],
            sorted(e['quotient_cycle_types']),e['local_index'],e['quotient_order'],e['states']],separators=(',',':'))
        signatures[key]+=1
    return dict(components=len(entries),states=sum(e['states'] for e in entries),
        by_degree=dict(Counter(e['degree'] for e in entries)),
        by_genus=dict(Counter(e['genus'] for e in entries)),
        by_quotient_genus=dict(Counter(e['quotient_genus'] for e in entries)),
        by_component_size=dict(Counter(e['states'] for e in entries)),
        unequal_three_quotient_types=sum(not e['equal_three_quotient_types'] for e in entries),
        mixed_cycle_quotients=sum(not e['all_three_full_quotient_cycles'] for e in entries),
        structural_signatures=[dict(signature=json.loads(k),components=v) for k,v in sorted(signatures.items())],
        warning='Sampled components; equal signatures do not prove marked-action isomorphism or moduli identity.')


def old_profiles():
    project=Path(__file__).resolve().parents[2]
    old=json.loads((project/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components']
    recent=json.loads(Path('data/discovery_design_catalog_20261013/CATALOG.json').read_text())['components']
    out=[]
    for e in old+recent:
        folder=Path(e['raw_folder'])
        if sha(folder/'STATES.tsv')!=e['table_sha256'] or sha(folder/'SUMMARY.tsv')!=e['summary_sha256']:
            raise ValueError('old source binding')
        p=profile(e['representative'],e['degree'],e['blocks'],e['block_size'])
        if p['group_order']!=int(e['group_order']) or p['genus']!=e['genus']:raise ValueError('old group/genus profile')
        out.append(dict(p,states=e['states'],raw_folder=e['raw_folder'],table_sha256=e['table_sha256']))
    return out
