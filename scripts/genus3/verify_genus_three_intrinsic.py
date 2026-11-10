"""Independent exact coordinate/population audit; closure is checked separately."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import itertools
import json
from pathlib import Path
import time

from genus_three_intrinsic import IntrinsicFrame,TRANS,cycles,compose
from verify_smooth_geometry import parse_tuple
from run_genus_three import sha,write


def rows(path):
    with Path(path).open()as stream:return list(csv.DictReader(stream,delimiter='\t'))


def numbers(text):return tuple(map(int,text.split(',')))if text else()


def connected(*perms):
    seen={0};queue=[0]
    for x in queue:
        for p in perms:
            y=p[x]
            if y not in seen:seen.add(y);queue.append(y)
    return len(seen)==len(perms[0])


def equivalent(pair,other):
    n=len(pair[0])
    for image in range(n):
        mapping={0:image};queue=[0];good=True
        for x in queue:
            for p,h in zip(pair,other):
                y=p[x];z=h[mapping[x]]
                if y in mapping:
                    if mapping[y]!=z:good=False;break
                elif z in mapping.values():good=False;break
                else:mapping[y]=z;queue.append(y)
            if not good:break
        if good and len(mapping)==n:return True
    return False


def support_orbit(beta,gamma,mask):
    def image(p,m):
        result=0
        while m:
            bit=m&-m;m-=bit;result|=1<<p[bit.bit_length()-1]
        return result
    seen={mask};queue=[mask]
    for m in queue:
        for p in(beta,gamma):
            y=image(p,m)
            if y not in seen:seen.add(y);queue.append(y)
    cycle_sets=[sum(1<<x for x in c)for p in(beta,gamma,compose(beta,gamma))for c in cycles(p)]
    good=all((m&c).bit_count()<=2 for m in seen for c in cycle_sets)
    return sorted(seen),good


def elliptic_betas(q):
    def partitions(total,least=1):
        if total==0:yield();return
        for k in range(least,min(total,q//2)+1):
            for rest in partitions(total-k,k):yield(k,)+rest
    result=[]
    for parts in partitions(q):
        if len(parts)>q//3:continue
        beta=[];at=0
        for length in parts:beta.extend(at+(j+1)%length for j in range(length));at+=length
        result.append(tuple(beta))
    return result


def guard(deadline):
    if deadline is not None and time.time()>=deadline:raise TimeoutError('absolute research deadline in independent intrinsic audit')


def quotient_population(q,records,all_elliptic=False,deadline=None):
    expected_shape=(2,)*4 if q==8 else(3,)*(q//3)
    delta_shape=(4,4)if q==8 else(3,)*(q//3)
    pairs=[(numbers(r['beta']),numbers(r['gamma']))for r in records]
    betas=elliptic_betas(q)if all_elliptic else[numbers(records[0]['beta'])]
    labeled=0;class_mass=[0]*len(pairs)
    for beta in betas:
        cb=len(cycles(beta));indices=[i for i,pair in enumerate(pairs)if sorted(map(len,cycles(pair[0])))==sorted(map(len,cycles(beta)))]
        # Independent explicit C_Sq(beta) orbits, not the production rooted
        # canonicalizer or a quadratic comparison of quotient representatives.
        groups={}
        for row in cycles(beta):groups.setdefault(len(row),[]).append(row)
        partial=[tuple(range(q))]
        for length,blocks in sorted(groups.items()):
            next_partial=[]
            for order in itertools.permutations(range(len(blocks))):
                for offsets in itertools.product(range(length),repeat=len(blocks)):
                    for old in partial:
                        h=list(old)
                        for i,block in enumerate(blocks):
                            target=blocks[order[i]]
                            for j,x in enumerate(block):h[x]=target[(j+offsets[i])%length]
                        next_partial.append(tuple(h))
            partial=next_partial
        orbits={}
        for i in indices:
            if i%1024==0:guard(deadline)
            if pairs[i][0]!=beta:raise ValueError('noncanonical beta profile')
            for h in partial:
                gamma=[0]*q
                for x in range(q):gamma[h[x]]=h[pairs[i][1][x]]
                key=tuple(gamma)
                if key in orbits and orbits[key]!=i:raise ValueError('duplicate quotient class')
                orbits[key]=i
        seen=set()
        for ordinal,gamma in enumerate(itertools.permutations(range(q))):
            if ordinal%1024==0:guard(deadline)
            cg=cycles(gamma)
            if all_elliptic:
                if len(cg)<cb or len(cg)+2*cb>q or max(map(len,cg))>q//2:continue
                cd=cycles(compose(beta,gamma))
                if len(cd)<cb or cb+len(cg)+len(cd)!=q or max(map(len,cd))>q//2:continue
            else:
                if tuple(sorted(map(len,cg)))!=expected_shape:continue
                if tuple(sorted(map(len,cycles(compose(beta,gamma)))))!=delta_shape:continue
            if not connected(beta,gamma):continue
            if gamma not in orbits:raise ValueError('quotient population gap')
            seen.add(gamma);labeled+=1;class_mass[orbits[gamma]]+=1
        if len(seen)!=len(orbits):raise ValueError('extra inadmissible quotient class')
    if not all(class_mass):raise ValueError('empty quotient class')
    return dict(labeled_quotients=labeled,class_mass=class_mass)


def verify(folder,check_population=True,deadline=None):
    folder=Path(folder);population=rows(folder/'POPULATION.tsv')[0];q=int(population['quotient_degree']);n=3*q
    all_elliptic=population.get('all_elliptic')=='true'
    if(q not in(3,6,8,9)and not(all_elliptic and 6<=q<=11))or population['complete']!='true':raise ValueError('incomplete/out-of-scope intrinsic population')
    records=rows(folder/'QUOTIENTS.tsv');cases=iter(rows(folder/'CASES.tsv'));coords=iter(rows(folder/'COORDINATES.tsv'))
    if len(records)!=int(population['quotient_frames']):raise ValueError('quotient frame mass')
    expected_qualified={};attempted=solved=transitive=smooth=0;support_cases=0;eligible_cases=0
    for fi,row in enumerate(records):
        if fi%32==0:guard(deadline)
        if int(row['quotient'])!=fi:raise ValueError('quotient order')
        f=IntrinsicFrame(numbers(row['beta']),numbers(row['gamma']))
        for name in('primal','dual','handles'):
            if tuple(getattr(f,name))!=numbers(row[name]):raise ValueError('coordinate tree marking')
        visited=set();masks=[1]if q==3 else[m for m in range(1<<q)if m.bit_count()==4]
        for mask in masks:
            if mask in visited:continue
            orbit,eligible=([1],True)if q==3 else support_orbit(f.beta,f.gamma,mask);visited.update(orbit)
            case=next(cases,None)
            if case is None or(int(case['quotient']),int(case['support']),numbers(case['orbit']),case['eligible'])!=(fi,mask,tuple(orbit),str(eligible).lower()):raise ValueError('support orbit coverage')
            support_cases+=1;eligible_cases+=eligible;cparam=csolved=ctrans=csm=0
            support=[j for j in range(q)if mask&(1<<j)]
            if eligible:
                for code in range(1 if q==3 else 81):
                    local=[0]*q
                    for j in support:
                        if q==3:local[j]=3
                        else:local[j]=TRANS[code%3];code//=3
                    for h in itertools.product(range(6),repeat=2):
                        attempted+=1;cparam+=1;t=f.solve(local,h)
                        if t is None:continue
                        solved+=1;csolved+=1;tr=connected(*t[:3]);sm=(n-len(cycles(compose(t[0],t[1]))))in(n-len(cycles(t[0]))+n-len(cycles(t[1])),n-len(cycles(t[2]))+n-len(cycles(t[3])))
                        actual=next(coords,None)
                        if actual is None:raise ValueError('missing solved coordinate')
                        marking=(int(actual['candidate']),int(actual['quotient']),int(actual['support']),numbers(actual['local']),numbers(actual['handles']))
                        if marking!=(attempted,fi,mask,tuple(local),h):raise ValueError('coordinate identity/order')
                        if parse_tuple(actual['seed'],n)!=t or actual['transitive']!=str(tr).lower()or actual['smooth']!=str(sm).lower():raise ValueError('coordinate lift/status mismatch')
                        if tr:transitive+=1;ctrans+=1
                        if tr and sm:smooth+=1;csm+=1;expected_qualified[attempted]=actual['seed']
            if tuple(map(int,(case['parameters'],case['solved'],case['transitive'],case['smooth'])))!=(cparam,csolved,ctrans,csm):raise ValueError('support case counts')
    if next(cases,None)is not None or next(coords,None)is not None:raise ValueError('unexpected extra coordinate/case')
    for key,count in dict(parameters=attempted,solved=solved,transitive=transitive,smooth=smooth).items():
        if int(population[key])!=count:raise ValueError('population '+key)
    resolutions=rows(folder/'RESOLUTIONS.tsv')
    if {int(r['candidate']):r['seed']for r in resolutions}!=expected_qualified or len(resolutions)!=len(expected_qualified):raise ValueError('qualified resolution mass')
    complete_quotients=quotient_population(q,records,all_elliptic,deadline)if check_population else None
    return dict(passed=True,quotient_degree=q,all_elliptic=all_elliptic,source_genus=2 if q==3 else 3,quotient_frames=len(records),support_cases=support_cases,eligible_cases=eligible_cases,
        parameters=attempted,solved=solved,transitive=transitive,smooth=smooth,quotient_population=complete_quotients,
        audit_software={str(p):sha(p)for p in(Path(__file__),Path(__file__).with_name('genus_three_intrinsic.py'))},
        bindings={name:sha(folder/name)for name in('POPULATION.tsv','QUOTIENTS.tsv','CASES.tsv','COORDINATES.tsv','LEADS.tsv','RESOLUTIONS.tsv','COUNTS.tsv')},
        qualification='exact finite coordinate and quotient/support coverage; all-smoothness and geometric nonisotriviality require separate certificates')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('--skip-quotient-population',action='store_true');p.add_argument('--deadline-epoch',type=float,default=1792762200);a=p.parse_args()
    try:check=verify(a.folder,not a.skip_quotient_population,a.deadline_epoch)
    except Exception as error:
        write(a.folder/'COORDINATE_AUDIT_FAILED.json',dict(passed=False,error=repr(error),time=time.time()));raise
    write(a.folder/'COORDINATE_AUDIT.json',check)
    print(json.dumps({k:v for k,v in check.items()if k not in('bindings','quotient_population','qualification')}))


if __name__=='__main__':main()
