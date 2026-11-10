"""Bounded, component-deduplicated replay of VERIFIED forest timeouts.

Overlapping partial quotient-state graphs prove two starts share a component.
Literal paths transport negative witnesses to every original start; the separate
native kernel checks all transported words. Partial/unknown results stay leads.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter, deque
from concurrent.futures import ThreadPoolExecutor
import csv
import fcntl
import json
from pathlib import Path
import subprocess
import sys

from audit_guided_native import audit_guided_native
from recover_forest_interruption import ensure_inactive, COUNT_HEADER, HEADER
from run_smooth_degree9 import sha, save
from run_s18_linked_roots import state_key
from odd_order_spin import parse_seed

PROJECT=Path(__file__).resolve().parents[2]


def collect(root):
    starts=[]; graph={}; owners={}; parent=[]; sources=[]
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def join(i,j):
        i,j=find(i),find(j)
        if i!=j:parent[max(i,j)]=min(i,j)
    for job in sorted(root.glob('guided_*forest')):
        if not (job/'VERIFIED').is_file():continue
        checks_path=job/'CHECKS.json'; checks=json.loads(checks_path.read_text())
        n=int(checks['construction']['degree'])
        partial={e['summary']['seed']:e for e in checks['components'] if e['summary']['status'] in ['timeout','capped']}
        with (job/'RESOLUTIONS.tsv').open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if row['status'] not in ['timeout','capped']:continue
                e=partial[row['seed']]; path=job/row['folder']/'STATES.tsv'
                if sha(path)!=e['table_sha256']:raise ValueError('partial source changed')
                with path.open() as table:rows=list(csv.DictReader(table,delimiter='\t'))
                keys=[(n,state_key(parse_seed(r['tuple'],n))) for r in rows]
                i=len(starts);parent.append(i)
                root_key=(n,state_key(parse_seed(row['seed'],n)))
                if root_key not in keys:raise ValueError('partial table misses root')
                starts.append(dict(degree=n,seed=row['seed'],root=root_key,states=len(rows),
                    source=str(job),candidate=row['candidate'],folder=row['folder'],table_sha256=e['table_sha256']))
                for key in keys:
                    if key in owners:join(i,owners[key])
                    else:owners[key]=i
                    graph.setdefault(key,{})
                for j,r in enumerate(rows):
                    edges=r['edges'].split(',')
                    if len(edges)!=6:raise ValueError('partial edge shape')
                    for m,target in enumerate(edges):
                        if target=='-':continue
                        q=keys[int(target)]
                        for a,b,mv in [(keys[j],q,m),(q,keys[j],(m+3)%6)]:
                            previous=graph[a].get(mv)
                            if previous is not None and previous!=b:raise ValueError('overlap graph is not equivariant')
                            graph[a][mv]=b
                sources.append(dict(path=str(path),sha256=e['table_sha256'],checks=str(checks_path),checks_sha256=sha(checks_path)))
    groups={}
    for i in range(len(starts)):groups.setdefault(find(i),[]).append(i)
    return starts,graph,list(groups.values()),sources


def paths_from(graph,root):
    paths={root:[]};queue=deque([root])
    while queue:
        x=queue.popleft()
        for m,y in graph[x].items():
            if y not in paths:paths[y]=paths[x]+[m];queue.append(y)
    return paths


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','out','binary','auditor']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--workers',type=int,default=6);p.add_argument('--seconds',type=int,default=20)
    args=p.parse_args();root=args.root.resolve()
    if not 1<=args.workers<=6 or args.seconds<1:raise ValueError('bounded replay budget')
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    ensure_inactive(root,int(json.loads((root/'LAUNCH.json').read_text())['pid']))
    starts,graph,groups,sources=collect(root);args.out.mkdir()
    save(args.out/'INPUTS.json',dict(starts=starts,groups=groups,sources=sources,
        binary_sha256=sha(args.binary),auditor_sha256=sha(args.auditor)))
    def work(item):
        group_id,ids=item;representative=max(ids,key=lambda i:starts[i]['states']);s=starts[representative]
        seed=args.out/f'seed_{group_id:05}.txt';seed.write_text(s['seed']+'\n')
        job=args.out/f'replay_{group_id:05}'
        subprocess.run([str(args.binary),'--d',str(s['degree']),'--replay',str(seed),'--seconds',str(args.seconds),'--out',str(job)],check=True,capture_output=True,text=True)
        checks=audit_guided_native(job,args.auditor,native_components=True)
        save(job/'CHECKS.json',checks);(job/'VERIFIED').write_text('Independent native replay\n')
        with (job/'RESOLUTIONS.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
        output=dict(representative=representative,original_indices=ids,status=row['status'],job=str(job),checks_sha256=sha(job/'CHECKS.json'))
        if row['status']=='nonsmooth':
            witness=list(map(int,row['word'].split(',')));paths=paths_from(graph,s['root']);proofs=[]
            for i in ids:
                path=paths[starts[i]['root']]
                word=[(m+3)%6 for m in reversed(path)]+witness
                proofs.append(dict(**starts[i],word=word))
            output['proofs']=proofs
        return output
    with ThreadPoolExecutor(max_workers=args.workers) as pool:results=list(pool.map(work,enumerate(groups)))
    proofs=[e for r in results for e in r.get('proofs',[])];native=[]
    for n in sorted({e['degree'] for e in proofs}):
        entries=[e for e in proofs if e['degree']==n];folder=args.out/f'original_witnesses_{n}';folder.mkdir()
        (folder/'LEADS.tsv').write_text(''.join(f"{i}\t{e['seed']}\n" for i,e in enumerate(entries,1)))
        (folder/'RESOLUTIONS.tsv').write_text(HEADER+''.join(f"{i}\tnonsmooth\t{','.join(map(str,e['word']))}\t\t{e['seed']}\n" for i,e in enumerate(entries,1)))
        count=len(entries);(folder/'COUNTS.tsv').write_text(COUNT_HEADER+f'{n}\treplay\t0\t{count}\t{count}\t{count}\t{count}\t0\t0\t0\t0\t0\n')
        checks=audit_guided_native(folder,args.auditor);save(folder/'CHECKS.json',checks)
        if checks['independent_native']['nonsmooth']!=count:raise ValueError('transported witness mass')
        native.append(dict(path=str(folder/'CHECKS.json'),sha256=sha(folder/'CHECKS.json'),records=count))
    report=dict(passed=True,original_occurrences=len(starts),partial_overlap_groups=len(groups),
        states_in_partial_union=len(graph),outcomes=dict(Counter(r['status'] for r in results)),
        original_starts_certified_negative=len(proofs),results=results,original_witness_audits=native,
        input_sha256=sha(args.out/'INPUTS.json'),scope='All current VERIFIED forest timeouts only. Unknown/partial cases remain unresolved. No whole-degree exclusion.')
    save(args.out/'RESULTS.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['results','original_witness_audits']}),flush=True)


if __name__=='__main__':main()
