"""Incremental structural/marked-quotient audit of immutable smooth certificates.

Checks a deterministic marked quotient projection on the root pure-braid fiber.
A nonfaithful projection disproves that *unmarked* quotient data alone labels
the states. It does not disprove a richer marked elliptic translation model.
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

from profile_smooth_theory import profile
from profile_smooth_actions import action_profile
from verify_smooth_geometry import parse_tuple,cycles
from run_smooth_degree9 import sha
from verify_smooth_geometry import forced_partition,compose
from run_s18_linked_roots import state_key

_BASE_COMPONENTS={}

def install_bases(entries):
    for e in entries:
        if e['degree'] not in [9,12]:continue
        folder=Path(e['raw_folder'])
        if sha(folder/'STATES.tsv')!=e['table_sha256']:raise ValueError('base certificate changed')
        with (folder/'STATES.tsv').open() as stream:keys={state_key(parse_tuple(row['tuple'],e['degree'])) for row in csv.DictReader(stream,delimiter='\t')}
        _BASE_COMPONENTS.setdefault(e['degree'],[]).append((e,keys))

def double_quotient(t):
    n=len(t[0]);parts={forced_partition(t,[[0,x]]) for x in range(1,n)}
    parts=sorted(p for p in parts if len(p)*2==n and all(len(b)==2 for b in p))
    if len(parts)!=1:raise ValueError('double quotient partition missing/nonunique; no exclusion')
    blocks=parts[0];membership={x:j for j,b in enumerate(blocks) for x in b}
    q=tuple(tuple(membership[g[b[0]]] for b in blocks) for g in t)
    return q,blocks

def ramified_double_profile(item):
    folder=Path(item['raw_folder']);n=item['degree'];r=n//2
    if n%2 or r not in _BASE_COMPONENTS:raise ValueError('no certified double-lift base in this queue')
    if sha(folder/'STATES.tsv')!=item['table_sha256']:raise ValueError('double-lift source changed')
    with (folder/'STATES.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    keys=[];quotients=[];excess_types=Counter();part0=None;q0=None
    for row in rows:
        t=parse_tuple(row['tuple'],n);q,parts=double_quotient(t)
        if part0 is None:part0,q0=parts,q
        excess=tuple(n-len(cycles(g))-2*(r-len(cycles(h))) for g,h in zip(t,q))
        if sorted(excess)!=[0,0,1,1]:raise ValueError('not the two-simple-ramification stratum')
        genus=1-r+sum(r-len(cycles(h)) for h in q)//2
        if int(row['genus'])!=2*genus:raise ValueError('double relative genus mismatch')
        ip=r-len(cycles(compose(q[0],q[1])))
        if ip not in [(r-len(cycles(q[0])))+(r-len(cycles(q[1]))),
                       (r-len(cycles(q[2])))+(r-len(cycles(q[3])))]:raise ValueError('quotient not smooth')
        keys.append(state_key(q));quotients.append(q);excess_types[excess]+=1
    image=set(keys);matches=[(e,k) for e,k in _BASE_COMPONENTS[r] if image==k]
    if len(matches)!=1:raise ValueError('source quotient not exactly one known complete base component')
    base,_=matches[0]
    from verify_smooth_geometry import hurwitz
    for i,row in enumerate(rows):
        edges=list(map(int,row['edges'].split(',')))
        if len(edges)!=6:raise ValueError('double-lift six-edge table')
        for m,j in enumerate(edges):
            if state_key(hurwitz(quotients[i],m))!=keys[j]:raise ValueError('double-lift projection not braid equivariant')
    fibers=Counter(keys)
    if len(set(fibers.values()))!=1:raise ValueError('nonuniform transitive projection fibers')
    from sympy.combinatorics import Permutation,PermutationGroup
    t0=parse_tuple(rows[0]['tuple'],n)
    order=int(PermutationGroup([Permutation(list(p)) for p in t0]).order())
    qorder=int(PermutationGroup([Permutation(list(p)) for p in q0]).order())
    if order%qorder:raise ValueError('double quotient group order')
    kernel_code=None
    if n==18 and order==41472 and qorder==324:
        supports=[{x for x in range(r) if g[x]!=x} for g in q0];support=min((s for s in supports if s),key=len)
        coarse=forced_partition(q0,[sorted(support)])
        if len(coarse)==3 and all(len(b)==3 for b in coarse):
            membership={x:j for j,b in enumerate(part0) for x in b};kernel=set()
            group=PermutationGroup([Permutation(list(p)) for p in t0])
            for element in group.generate_schreier_sims():
                if any(membership[element(b[0])]!=j for j,b in enumerate(part0)):continue
                kernel.add(sum((element(b[0])==b[1])<<j for j,b in enumerate(part0)))
            expected={mask for mask in range(1<<r) if len({sum((mask>>x)&1 for x in b)%2 for b in coarse})==1}
            kernel_code=dict(order=len(kernel),rank=7,equal_parities_on_three_coarse_blocks=kernel==expected,
                coarse_blocks=[[x+1 for x in b] for b in coarse])
            if len(kernel)!=order//qorder:raise ValueError('enumerated double kernel mass')
    return dict(degree=n,states=len(rows),genus=int(rows[0]['genus']),mechanism='ramified-double-lift',
        primitive=False,imprimitive=True,blocks=r,block_size=2,
        block_partition=[[x+1 for x in b] for b in part0],group_order=order,quotient_order=qorder,
        kernel_order=order//qorder,base_degree=r,base_genus=base['genus'],base_states=base['states'],
        kernel_code=kernel_code,
        base_table_sha256=base['table_sha256'],quotient_projection_fiber=next(iter(fibers.values())),
        all_rows_and_six_quotient_edges_checked=True,
        ramification_excess_histogram={','.join(map(str,k)):v for k,v in excess_types.items()},
        raw_folder=str(folder),table_sha256=item['table_sha256'],representative=rows[0]['tuple'],
        constructor_seed=item['representative'],partition_for='STATES.tsv state0',
        action=action_profile(folder,item['table_sha256']),
        geometric_classification='not inferred; relative branch locations require separate analysis')


def primitive(t):
    n=len(t[0])
    for x in range(1,n):
        labels=list(range(n));labels[x]=0
        while True:
            changed=False
            for a in range(n):
                for b in range(a+1,n):
                    if labels[a]!=labels[b]:continue
                    for p in t[:3]:
                        i,j=labels[p[a]],labels[p[b]]
                        if i!=j:
                            lo,hi=min(i,j),max(i,j)
                            labels=[lo if v==hi else v for v in labels];changed=True
            if not changed:break
        size=labels.count(labels[0])
        if 1<size<n:return False
    return True


def marked_quotient(seed,n,r,k):
    p=profile(seed,n,r,k,compute_orders=False);t=parse_tuple(seed,n)
    blocks=[[x-1 for x in b] for b in p['block_partition']]
    membership={x:j for j,b in enumerate(blocks) for x in b}
    q=[tuple(membership[g[b[0]]] for b in blocks) for g in t]
    local=next(i for i,g in enumerate(q) if g==tuple(range(r)))
    mark=membership[next(x for x in range(n) if t[local][x]!=x)]
    triple=[g for i,g in enumerate(q) if i!=local];answers=[]
    for root in range(r):
        order=[root];labels={root:0}
        for x in order:
            for g in triple:
                if g[x] not in labels:labels[g[x]]=len(order);order.append(g[x])
        if len(order)!=r:raise ValueError('nontransitive quotient triple')
        answers.append((tuple(tuple(labels[g[x]] for x in order) for g in triple),labels[mark]))
    return min(answers)


def projection_profile(entry,max_fiber=256):
    folder=Path(entry['raw_folder']);n=entry['degree'];r=entry['blocks'];k=entry['block_size']
    if sha(folder/'STATES.tsv')!=entry['table_sha256']:raise ValueError('projection table changed')
    with (folder/'STATES.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    adj=[list(map(int,row['edges'].split(','))) for row in rows]
    words=[list(range(j-1,i,-1))+[i,i]+[m+3 for m in range(i+1,j)] for i in range(3) for j in range(i+1,4)]
    def walk(s,w):
        for m in w:s=adj[s][m]
        return s
    seen={0};queue=[0];edges={}
    for s in queue:
        edges[s]=[walk(s,w) for w in words]
        for h in edges[s]:
            if h not in seen:
                seen.add(h);queue.append(h)
                if len(queue)>max_fiber:return dict(status='fiber_cap',lower_bound=len(queue),not_an_exclusion=True)
    marks={s:marked_quotient(rows[s]['tuple'],n,r,k) for s in queue}
    triples={v[0] for v in marks.values()};image=set(marks.values());transitions={};collisions=0
    for s in queue:
        for j,h in enumerate(edges[s]):
            key=(marks[s],j)
            if key in transitions and transitions[key]!=marks[h]:collisions+=1
            else:transitions[key]=marks[h]
    return dict(status='checked',pure_fiber=len(queue),quotient_triple_classes=len(triples),
        marked_quotient_classes=len(image),faithful=len(image)==len(queue),
        deterministic_factor=collisions==0,transition_inconsistencies=collisions,
        projection_fiber_multiplicities=dict(Counter(Counter(marks.values()).values())),
        interpretation='Exact root-pure-fiber test of this canonical marked block projection; not geometric isotriviality or a universal factor theorem.')


def classify(item):
    folder=Path(item['raw_folder'])
    for name,field in [('STATES.tsv','table_sha256'),('SUMMARY.tsv','summary_sha256')]:
        if field in item and sha(folder/name)!=item[field]:raise ValueError('immutable profile input changed')
    try:
        p=profile(item['representative'],item['degree'],compute_orders=False)
        p.update(states=item['states'],raw_folder=str(folder),table_sha256=item['table_sha256'])
        p['point_group_order']=item.get('group_order')
        p['point_group_order_provenance']=item.get('group_order_provenance','producer only; not independently rechecked here')
        p['action']=item.get('action') or action_profile(folder,item['table_sha256'])
        p['marked_projection']=projection_profile(p)
        p['mechanism']='single-local/unramified'
    except ValueError as error:
        # Failure to recover this particular factor is not proof of primitivity.
        try:p=ramified_double_profile(item)
        except ValueError as lift_error:
            p=dict(degree=item['degree'],states=item['states'],raw_folder=str(folder),
                table_sha256=item['table_sha256'],representative=item['representative'],
                mechanism='unclassified',structural_error=str(error),lift_error=str(lift_error),
                primitive=primitive(parse_tuple(item['representative'],item['degree'])))
    return p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queue',type=Path,required=True);parser.add_argument('--start',type=int,default=0)
    parser.add_argument('--limit',type=int,default=10);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--all-block-systems',action='store_true')
    parser.add_argument('--batch-only',action='store_true')
    parser.add_argument('--lattice-only',action='store_true',help='all invariant partitions of bound representatives; no repeated table/action audit')
    args=parser.parse_args();args.out.mkdir();data=json.loads(args.queue.read_text())
    entries=data['components'];results=[]
    if not args.lattice_only:install_bases(entries)
    for i in range(args.start,min(len(entries),args.start+args.limit)):
        item=entries[i]
        result=(dict(degree=item['degree'],states=item['states'],representative=item['representative'],
                     raw_folder=item['raw_folder'],table_sha256=item['table_sha256'],mechanism='lattice-only')
                if args.lattice_only else classify(item))
        result['queue_index']=i;results.append(result)
        if args.all_block_systems or args.lattice_only:
            from smooth_block_lattice import block_lattice
            result['block_lattice']=block_lattice(parse_tuple(result['representative'],result['degree']))
        if not args.batch_only:(args.out/f'profile_{i:06}.json').write_text(json.dumps(result,indent=2)+'\n')
        if not args.lattice_only or i%100==0:
            print('PROFILE',i,result['degree'],result['states'],result['mechanism'],result.get('marked_projection'),flush=True)
    report=dict(queue_sha256=sha(args.queue),start=args.start,next=args.start+len(results),total=len(entries),
        complete=args.start+len(results)==len(entries),components=results)
    (args.out/'PROFILES.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
