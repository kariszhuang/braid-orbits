#!/usr/bin/env python3
"""Frozen, shared-lock exact smooth-orbit mission. No weak-search filters."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from analyze_fixed_braid_constraints import point_orbits, product
from odd_order_spin import cycle_lengths, hurwitz, parse_seed
from run_s18_linked_roots import certified_states, state_key

PROJECT=Path(__file__).resolve().parents[2]
REPORTED='[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]'
JOBS=[dict(m=m,fixed=f,shape=s) for m,f in [(2,1),(3,0),(3,1),(4,0),(4,1)]
      for s in ['333','331'] if (s=='333' and m>=3) or (s=='331' and f==1)]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,data):
    with Path(p).open('x') as f:json.dump(data,f,indent=2);f.write('\n')
def rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def values(t):
    n=len(t[0]);indices=[n-len(cycle_lengths(p)) for p in t]
    ip=n-len(cycle_lengths(product(*t[:2])))
    return dict(smooth=ip in [sum(indices[:2]),sum(indices[2:])],
                defect=n-ip-len(point_orbits(t[:2]))-len(point_orbits(t[2:]))+1,
                genus=1-n+sum(indices)//2,product_index=ip,
                left_orbits=len(point_orbits(t[:2])),right_orbits=len(point_orbits(t[2:])))


def audit_component(folder, collect_records=False):
    summary=rows(folder/'SUMMARY.tsv')[0];n=int(summary['degree'])
    states=rows(folder/'STATES.tsv');keys={};records=[]
    for row in states:
        t=parse_seed(row['tuple'],n)
        if product(*t)!=tuple(range(n)) or len(point_orbits(t))!=1:raise ValueError('invalid tuple')
        key=state_key(t)
        if key in keys:raise ValueError('duplicate conjugacy state')
        keys[key]=int(row['state']);v=values(t)
        if collect_records:records.append((key,row['tuple']))
        if str(v['smooth']).lower()!=row['smooth']:raise ValueError('smoothness mismatch')
        for name in ['defect','genus','product_index','left_orbits','right_orbits']:
            if v[name]!=int(row[name]):raise ValueError(name+' mismatch')
        if v['smooth'] and v['defect']!=0:raise ValueError('smooth-implies-zero violation')
    closed=summary['status']=='closed_smooth'
    if closed:
        if not all(r['smooth']=='true' for r in states):raise ValueError('false smooth certificate')
        adjacency=[]
        for row in states:
            t=parse_seed(row['tuple'],n);edges=list(map(int,row['edges'].split(',')))
            if edges!=[keys[state_key(hurwitz(t,m))] for m in range(6)]:raise ValueError('edge mismatch')
            adjacency.append(edges)
        reached={0};queue=[0]
        for i in queue:
            for j in adjacency[i]:
                if j not in reached:reached.add(j);queue.append(j)
        if len(reached)!=len(states):raise ValueError('disconnected certificate')
        if state_key(parse_seed(summary['seed'],n)) not in keys:raise ValueError('wrong root')
    result=dict(summary=summary,independent_six_move_closure=closed,
                state_keys=sorted(keys),table_sha256=sha(folder/'STATES.tsv'))
    if collect_records and closed:result['_closed_records']=records
    return result


def audit_job(folder):
    coverage=rows(folder/'COVERAGE.tsv')[0];n=int(coverage['degree'])
    leads={r['candidate']:parse_seed(r['seed'],n) for r in rows(folder/'LEADS.tsv')}
    resolutions=rows(folder/'RESOLUTIONS.tsv');status=Counter()
    for row in resolutions:
        t=leads[row['candidate']];status[row['status']]+=1
        if row['status']=='nonsmooth':
            for mv in row['word'].split(','):
                if mv:t=hurwitz(t,int(mv))
            if values(t)['smooth']:raise ValueError('false nonsmooth word')
    if len(resolutions)!=len(leads):raise ValueError('unresolved durable promoted records')
    components=[audit_component(p) for p in sorted(folder.glob('component_*'))]
    for c in components:c.pop('state_keys')
    if len(components)!=int(coverage['closed_components']):raise ValueError('component count mismatch')
    if status['nonsmooth']!=int(coverage['nonsmooth']):raise ValueError('rejection count mismatch')
    return dict(coverage=coverage,literal_resolution_counts=dict(status),components=components)


def weak_audit(directory):
    answer=[]
    for path in sorted(directory.glob('CERTIFIED_*weak.txt')):
        states=list(certified_states(path));v=[values(t) for t in states]
        answer.append(dict(certificate=path.name,states=len(states),smooth_states=sum(x['smooth'] for x in v),
                           defect_histogram=dict(Counter(x['defect'] for x in v)),sha256=sha(path)))
    if len(answer)!=9 or any(r['smooth_states'] for r in answer):raise ValueError('weak catalog control failed')
    return answer


def worker(root,lock_fd):
    # Inherited open file description retains the campaign lock across detach.
    lock=os.fdopen(lock_fd,'a');plan=json.loads((root/'PLAN.json').read_text())
    for name,digest in plan['source_hashes'].items():
        if sha(root/'source'/name)!=digest:raise ValueError('frozen source mismatch')
    binary=root/'smooth_orbits'
    if sha(binary)!=plan['binary_sha256']:raise ValueError('binary mismatch')
    save(root/'STATE.json',dict(status='running',pid=os.getpid(),started=dt.datetime.now().isoformat()))
    source=root/'source';reported=root/'reported_A9'
    subprocess.run([str(binary),'--mode','profile','--d','9','--tuple',REPORTED,'--cap','50000','--seconds','60','--out',str(reported)],check=True)
    check=audit_component(reported);check.pop('state_keys')
    if check['summary']['states']!='12' or check['summary']['genus']!='2' or check['summary']['group_order']!='324' or check['summary']['full_group']!='false' or not check['independent_six_move_closure']:raise ValueError('reported positive control failed')
    # Independent literal enumeration, feasible for this 324-element control.
    t=parse_seed(REPORTED,9);identity=tuple(range(9));group={identity};queue=[identity]
    for h in queue:
        for g in t:
            q=product(g,h)
            if q not in group:group.add(q);queue.append(q)
    if len(group)!=324:raise ValueError('independent reported group-order mismatch')
    check['independent_literal_group_order']=len(group)
    save(root/'CONTROLS.json',dict(reported=check,weak_catalog=weak_audit(source/'evidence/certificates')))
    def run_job(j):
        name=f"A{3*j['m']+j['fixed']}_{j['shape']}";folder=root/name
        cmd=[str(binary),'--mode','enumerate','--m',str(j['m']),'--fixed',str(j['fixed']),'--shape',j['shape'],
             '--out',str(folder),'--seconds',str(plan['seconds']),'--cap','50000']
        with (root/(name+'.log')).open('x') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        audit=audit_job(folder);save(root/(name+'_CHECKS.json'),audit)
        print(name,json.dumps(audit['coverage']),flush=True);return audit
    def block_job(k):
        folder=root/f'block_family_{3*k}'
        subprocess.run([str(binary),'--mode','block-family','--block-size',str(k),'--seconds','30','--cap','50000','--out',str(folder)],check=True)
        answer=audit_component(folder);answer.pop('state_keys')
        if answer['summary']['genus']!='2':raise ValueError('block family genus')
        from math import factorial
        expected=3*factorial(k)**3//2
        if int(answer['summary']['group_order'])!=expected:raise ValueError('block family order prediction failed')
        if answer['summary']['status']=='nonsmooth':raise ValueError('block family smooth proof regression')
        save(root/f'block_family_{3*k}_CHECKS.json',answer);return answer
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        results=list(pool.map(run_job,plan['jobs']))
        blocks=list(pool.map(block_job,range(3,9)))
    save(root/'CHECKS.json',dict(controls=check,jobs=results,
        block_family=blocks,
        fully_resolved=all(r['coverage']['fully_resolved']=='true' for r in results),
        warning='Specific 3-cycle/3^m 1^f passports; not all alternating tuples. All transitive groups eligible, full A/S generation separately marked.'))
    save(root/'FINISHED.json',dict(status='complete',ended=dt.datetime.now().isoformat()))
    lock.close()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true')
    p.add_argument('--seconds',type=int,default=600);p.add_argument('--workers',type=int,default=5)
    p.add_argument('--worker',type=Path);p.add_argument('--lock-fd',type=int);args=p.parse_args()
    if args.worker:
        try:worker(args.worker,args.lock_fd)
        except Exception as error:
            save(args.worker/'FAILURE.json',dict(error=repr(error)));raise
        return
    if not args.run:print(json.dumps(dict(jobs=JOBS,workers=args.workers,seconds=args.seconds),indent=2));return
    if not 1<=args.workers<=7 or args.seconds<=0:p.error('1..7 workers, positive seconds')
    dest=Path('data')
    if not Path('data').is_dir() or shutil.disk_usage(dest).free<30*2**30:raise RuntimeError('T7 mounted and 30 GiB free required')
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    root=dest/('smooth_mission_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'));root.mkdir()
    source=root/'source';source.mkdir();hashes={}
    paths=[PROJECT/'Cargo.toml',PROJECT/'Cargo.lock',*sorted((PROJECT/'src').rglob('*.rs')),
           *[_scripts_path.find(s) for s in ['run_smooth_mission.py','run_s18_linked_roots.py','analyze_fixed_braid_constraints.py','odd_order_spin.py']],
           *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
           PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv',
           PROJECT/'evidence/smooth_mission_20261011/PLAN.md',PROJECT/'evidence/smooth_mission_20261011/MATH.md']
    for path in paths:
        relative=path.relative_to(PROJECT);target=source/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target);hashes[str(relative)]=sha(target)
    env=dict(os.environ,CARGO_BUILD_JOBS='1',CARGO_TARGET_DIR=str(PROJECT/'target/smooth-mission'))
    print('OUTPUT='+str(root),flush=True)
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo','test','--release','--bin','smooth_orbits','smooth_tests','--','--test-threads=1'],cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run(['cargo','build','--release','--bin','smooth_orbits'],cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    binary=root/'smooth_orbits';shutil.copy2(PROJECT/'target/smooth-mission/release/smooth_orbits',binary)
    save(root/'PLAN.json',dict(jobs=JOBS,workers=args.workers,seconds=args.seconds,source_hashes=hashes,binary_sha256=sha(binary),
         quotient='simultaneous S_n conjugation',objective='smooth everywhere, not weak pair failure',
         previous_campaign='data/campaign_20261010_212021_242369'))
    with (root/'WORKER.log').open('x') as log:
        child=subprocess.Popen([sys.executable,str(source/'scripts/smooth/run_smooth_mission.py'),'--worker',str(root),'--lock-fd',str(lock.fileno())],
            cwd=source,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json',dict(pid=child.pid,started=dt.datetime.now().isoformat()))
    lock.close();print('PID='+str(child.pid),flush=True)


if __name__=='__main__':main()
