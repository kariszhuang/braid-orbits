"""Independent finite group/frame coverage and literal omitted-seed verification.

This is an exact subgroup population, not all primitive groups of a degree.
The producer's class/centralizer reduction is reconstructed with Python kernels.
Every identity-free positive-genus transitive occurrence not represented in the
durable journals must have a literal nonsmooth braid translate or known closure.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from analyze_fixed_braid_constraints import product,point_orbits,conjugate
from odd_order_spin import parse_seed,inverse,hurwitz,cycle_lengths
from run_s18_linked_roots import state_key
from run_smooth_degree9 import sha,rows


def index(p):return len(p)-len(cycle_lengths(p))
def parse_p(s,n):return parse_seed('['+s+', (), (), ()]',n)[0]


def independent_frames(group,n):
    group=sorted(group);indices={p:index(p) for p in group};seen=set();frames=[]
    for a in group:
        if a in seen:continue
        seen.update(conjugate(a,h) for h in group)
        if indices[a]==0 or indices[a]>=n-1:continue
        center=[h for h in group if product(a,h)==product(h,a)];b_seen=set()
        for b in group:
            if b in b_seen or indices[b]==0 or indices[a]+indices[b]>=n:continue
            if index(product(a,b))!=indices[a]+indices[b]:continue
            b_seen.update(conjugate(b,h) for h in center);frames.append((a,b))
    return frames


WORDS=[[1],[4],[0,1],[0,4],[3,1],[1,1],[4,4],[1,1,1],[4,4,4]]
WORDS += [[end]*power+[middle] for end in [0,3] for power in range(2,5) for middle in [4,1]]
WORDS += [[1,3,1],[0,4,4],[4,0,4],[3,1,1],[5,1,1,1],[5,3,1,1,1]]


def smooth(t):
    ip=index(product(t[0],t[1]));return ip in [index(t[0])+index(t[1]),index(t[2])+index(t[3])]


def check(folder,out):
    folder=Path(folder);info=rows(folder/'PRIMITIVE.tsv')[0];n=int(info['degree'])
    group=[parse_p(line.split('\t')[1],n) for line in (folder/'GROUP.tsv').read_text().splitlines()]
    if len(group)!=int(info['group_order']) or group!=sorted(set(group)):raise ValueError('group order/order/uniqueness')
    if info['group'] in ['m11','psl11']:
        a=parse_p('(2 10)(4 11)(5 7)(8 9)',n);b=parse_p('(1 4 3 8)(2 5 6 9)',n)
        if info['group']=='psl11':b=product(b,a,b,b,a,b)
        gens=[a,b];expected=7920 if info['group']=='m11' else 660
    else:
        # Independently identify a small generating subset of the exported
        # finite action. The declared affine/projective name is not certified
        # merely by matching its order; only its literal subgroup is covered.
        gens=[];generated={tuple(range(n))}
        for g in group:
            if g in generated:continue
            gens.append(g);queue=list(generated)
            for x in queue:
                for h in gens:
                    y=product(h,x)
                    if y not in generated:generated.add(y);queue.append(y)
            if len(generated)==len(group):break
        expected=len(group)
    generated={tuple(range(n))};queue=list(generated)
    for x in queue:
        for g in gens:
            y=product(g,x)
            if y not in generated:generated.add(y);queue.append(y)
    if len(generated)!=expected or generated!=set(group):raise ValueError('independent generator subgroup mismatch')
    frames=independent_frames(group,n)
    actual=[tuple(parse_p(s,n) for s in line.split('\t')[1:]) for line in (folder/'PAIR_FRAMES.tsv').read_text().splitlines()]
    if frames!=actual or len(frames)!=int(info['frames']):raise ValueError('independent class/centralizer pair coverage')
    start,end,total=int(info['start']),int(info['next']),int(info['total'])
    if not 0<=start<=end<=total or total!=len(frames)*len(group):raise ValueError('finite rank mass')
    journal={};statuses=Counter()
    for row in rows(folder/'RESOLUTIONS.tsv'):
        t=parse_seed(row['seed'],n);journal[state_key(t)]=row['status'];statuses[row['status']]+=1
    known=set()
    if (folder/'KNOWN.tsv').exists():
        for line in (folder/'KNOWN.tsv').read_text().splitlines():known.add(state_key(parse_seed(line.split('\t',1)[1],n)))
    for row in rows(folder/'RESOLUTIONS.tsv'):
        if row['status']=='closed_smooth':
            for s in rows(folder/row['folder']/'STATES.tsv'):known.add(state_key(parse_seed(s['tuple'],n)))
    count=Counter();proofs=Counter()
    for rank in range(start,end):
        a,b=frames[rank//len(group)];c=group[rank%len(group)];d=inverse(product(a,b,c));t=(a,b,c,d)
        count['occurrences']+=1
        if any(g==tuple(range(n)) for g in t):count['identity']+=1;continue
        genus=1-n+sum(map(index,t))//2
        if genus<=0:count['nonpositive_genus']+=1;continue
        if len(point_orbits(t))!=1:count['intransitive']+=1;continue
        count['qualifying']+=1
        if not smooth(t):raise ValueError('non-geodesic constructed starting tuple')
        proof=None
        for word in WORDS:
            u=t
            for move in word:u=hurwitz(u,move)
            if not smooth(u):proof=word;break
        if proof is not None:
            count['literal_negative']+=1;proofs[','.join(map(str,proof))]+=1;continue
        key=state_key(t)
        if key not in journal and key not in known:raise ValueError('unaccounted surviving occurrence')
        count['accounted_survivor']+=1
    result=dict(group=info['group'],degree=n,group_order=len(group),frames=len(frames),start=start,next=end,total=total,
        independent_complete_population=start==0 and end==total,
        population='all additive nonidentity first pairs in this literal subgroup modulo its simultaneous conjugation, all group-valued third entries; identity-free positive-genus transitive scope',
        counts=dict(count),literal_word_histogram=dict(proofs),journal_statuses=dict(statuses),
        no_unresolved=not any(statuses[s] for s in ['capped','timeout']),
        passed=True,sources={name:sha(folder/name) for name in ['GROUP.tsv','PAIR_FRAMES.tsv','PRIMITIVE.tsv','LEADS.tsv','RESOLUTIONS.tsv','COUNTS.tsv']},
        warning='This does not cover every primitive group or ambient A/S degree. Nonnegative genus scope/transitivity and tuple quotient remain separate from geometry.')
    Path(out).write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();print(json.dumps(check(args.folder,args.out),indent=2),flush=True)
