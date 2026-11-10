"""Independent arithmetic audit of passport THEOREM exclusions.

No producer code, polynomial bitsets, braid closure or isotrivial classifier
is imported. A sampled theorem exclusion is not a literal negative word.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from functools import lru_cache
from itertools import combinations
from pathlib import Path

from verify_smooth_geometry import parse_tuple,compose,cycles
from run_genus_three import sha


def prime(n):return n>=2 and all(n%j for j in range(2,int(n**0.5)+1))


@lru_cache(maxsize=4096)
def arithmetic_rules(n,k):
    if not(prime(n)and 2<=k<n):return frozenset()
    rules=set()
    if k==3:rules.add(('prime-three-tiling',3,3))
    elif prime(k):rules.add(('prime-tile-valuation',k,k))
    if k==4:rules.add(('four-tile-prime-periodicity',2,4))
    for ell in range(2,k+1):
        if not(prime(ell)and k%ell==0):continue
        q=1;remaining=k
        while remaining%ell==0:q*=ell;remaining//=ell
        orbit={0,*[pow(ell,e,n)for e in range(n-1)]}
        run=any(all((start+j)%n in orbit for j in range(k))for start in range(n))
        if run and not any(modulus==q and method!='prime-tile-valuation'for method,p,modulus in rules):
            rules.add(('consecutive-root-run',ell,q))
    return frozenset(rules)


def violations(n,shapes):
    indices=[n-len(shape)for shape in shapes]
    if 0 in indices or sum(indices)!=2*n+4:return set()
    answer=set()
    for special,shape in enumerate(shapes):
        active=[length for length in shape if length!=1]
        if len(active)!=1:continue
        k=active[0];rules=arithmetic_rules(n,k)
        ordinary=[j for j in range(4)if j!=special]
        if any(indices[i]+indices[j]<=n-1 for i,j in combinations(ordinary,2)):continue
        for partner in ordinary:
            if len(shapes[partner])!=k:continue
            for method,p,q in rules:
                if method!='prime-tile-valuation'and n//k<=q:
                    answer.add((special,partner,k,p,q,method,'length-bound',None));continue
                for violator in ordinary:
                    if violator==partner:continue
                    bad=(any(length%p==0 for length in shapes[violator])if method=='prime-tile-valuation'else
                         len({length%q for length in shapes[violator]})>1)
                    if bad:answer.add((special,partner,k,p,q,method,'ordinary-cycle-divisibility'if method=='prime-tile-valuation'else'ordinary-cycle-congruence',violator))
    return answer


def check_row(n,row,product_one=True):
    t=parse_tuple(row['seed'],n);shapes=[[len(c)for c in cycles(p)]for p in t]
    if product_one:
        if compose(compose(t[0],t[1]),compose(t[2],t[3]))!=tuple(range(n)):raise ValueError('passport theorem sample product')
        if any(len(shape)==n for shape in shapes)or sum(n-len(shape)for shape in shapes)!=2*n+4:raise ValueError('passport theorem sample genus/identity scope')
    bad=violations(n,shapes)
    if row['status']=='not-rejected':
        if bad:raise ValueError('native passport classifier missed an independent rejection')
        return
    if row['status']!='theorem-excluded':raise ValueError('passport theorem status')
    item=tuple(int(row[key])for key in('special','partner','k','prime','modulus'))+(row['method'],row['reason'],int(row['violator'])if row['violator']else None)
    if item not in bad:raise ValueError('passport theorem reported unsupported hypothesis/conclusion')


def verify_saturated_exclusions(folder):
    folder=Path(folder)
    with(folder/'COUNTS.tsv').open()as stream:counts=next(csv.DictReader(stream,delimiter='\t'))
    with(folder/'SATURATED_SUMMARY.tsv').open()as stream:summary=next(csv.DictReader(stream,delimiter='\t'))
    with(folder/'SATURATED_THEOREMS.tsv').open()as stream:theorems=list(csv.DictReader(stream,delimiter='\t'))
    n=int(counts['degree']);enabled=summary['enabled']=='true';excluded=int(summary['theorem_excluded_occurrences']);sampled=int(summary['sampled_occurrences'])
    if summary['enabled']not in('true','false')or not 0<=excluded<=int(summary['checked_scope_occurrences'])<=int(counts['scope']):raise ValueError('passport theorem summary mass')
    if not enabled and excluded:raise ValueError('OFF passport flag excluded proposals')
    if sum(int(r['excluded_occurrences'])for r in theorems)!=excluded or sum(int(r['sampled_occurrences'])for r in theorems)!=sampled:raise ValueError('passport theorem aggregate mass')
    paths=['SATURATED_SUMMARY.tsv','SATURATED_THEOREMS.tsv'];rows=[]
    if enabled:
        paths.append('SATURATED_EXCLUSIONS.tsv')
        with(folder/paths[-1]).open()as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
        if not len(rows)==sampled<=min(64,excluded):raise ValueError('passport theorem bounded sample mass')
    seen=set()
    for row in rows:
        ident=int(row['candidate'])
        if ident in seen or not 0<ident<=int(counts['constructed']):raise ValueError('passport theorem candidate')
        seen.add(ident);check_row(n,row)
    traced=0
    if(folder/'SATURATED_TRACE.tsv').exists():
        paths.append('SATURATED_TRACE.tsv')
        with(folder/paths[-1]).open()as stream:
            for row in csv.DictReader(stream,delimiter='\t'):check_row(n,row);traced+=1
    if traced!=int(summary['trace_records']):raise ValueError('passport trace mass')
    return dict(passed=True,enabled=enabled,checked_scope_occurrences=int(summary['checked_scope_occurrences']),theorem_excluded_occurrences=excluded,
        sampled_occurrences=sampled,trace_records=traced,sources={name:sha(folder/name)for name in paths},
        qualification='goal-specific arithmetic THEOREM exclusions, not transitive literal negatives, isotrivial labels or population coverage')
