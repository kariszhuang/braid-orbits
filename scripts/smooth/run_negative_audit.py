#!/usr/bin/env python3
"""Resumable independent negative-only census audit. Plan-only by default.

Does not enumerate a new search population or modify predecessor evidence.
Positive/partial outcomes require the existing general Python verifier instead.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse, csv, datetime as dt, fcntl, json, os, shutil, subprocess, sys
import threading, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from zoneinfo import ZoneInfo
from run_smooth_degree9 import COLS, rows, sha, state, now
from run_smooth_frontier import audit_split_frames, audit_guided

PROJECT=Path(__file__).resolve().parents[2]
DEFAULT=Path('data/smooth_frontier_20261012_123427_295724')

def save(path,data):
    # Interrupted JSON writes cannot become trusted resumed frame reports.
    path=Path(path);temp=path.with_name(path.name+f'.tmp-{os.getpid()}')
    with temp.open('w') as stream:json.dump(data,stream,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    temp.replace(path)

class Native:
    def __init__(self,binary,n):
        self.p=subprocess.Popen([str(binary),'--server',str(n)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
        if self.p.stdout.readline().strip()!='READY':raise ValueError('native verifier did not start')
    def frame(self,folder,a,b,positive,limit=0):
        fields=[str(folder),','.join(map(str,a)),','.join(map(str,b)),str(int(positive)),str(limit)]
        self.p.stdin.write('\t'.join(fields)+'\n');self.p.stdin.flush()
        line=self.p.stdout.readline()
        if not line:raise RuntimeError('native verifier failed; no exclusion certified')
        return json.loads(line)
    def close(self):
        self.p.stdin.close()
        self.p.wait(timeout=30)
        self.p.stdout.close()
        if self.p.returncode:raise RuntimeError(f'native exit {self.p.returncode}')

def checked_frame(result,folder,fi,expected,positive,n):
    if result['count']!=expected:raise ValueError(f'frame {fi} coverage mismatch')
    local=rows(folder/'COUNTS.tsv')
    if len(local)!=1 or int(local[0]['frame'])!=fi or [int(local[0][c]) for c in COLS]!=expected:raise ValueError('local/global mismatch')
    from math import factorial
    if (folder/'COMPLETE').read_text()!=f'{fi}\t{factorial(n)}\n':raise ValueError('completion marker')
    if positive:
        row=rows(folder/'OMITTED.tsv')[0]
        if result['omitted']!=[int(row['identity']),int(row['nonpositive_genus'])]:raise ValueError('omission mismatch')
    if not result['complete']:raise ValueError('prefix cannot certify a frame')

def hashes(folder,positive):
    return {name:sha(folder/name) for name in ['OUTCOMES.bin','WITNESSES.tsv','COUNTS.tsv','COMPLETE']+(['OMITTED.tsv'] if positive else [])}

def audit(root,folder,binary,workers):
    coverage=rows(folder/'COVERAGE.tsv')[0];n=int(coverage['degree'])
    if coverage['enumeration_complete']!='true' or coverage['fully_resolved']!='true' or any(int(coverage[c]) for c in ['closed_occurrences','closed_components','unresolved']):
        raise ValueError('native negative-only scope not satisfied; use general independent Python audit')
    positive=(folder/'POSITIVE_ONLY').exists()
    state(root,status='verifying_frame_population',pid=os.getpid(),workers=1,gpu=False)
    frames=audit_split_frames(folder,n,None);counts=rows(folder/'FRAME_COUNTS.tsv')
    if len(counts)!=len(frames) or [int(r['frame']) for r in counts]!=list(range(len(frames))):raise ValueError('frame gaps/duplicates')
    out=root/'audit_frames';out.mkdir(exist_ok=True);total=[0]*8;weighted=[0,0];completed=0
    start=time.monotonic();queue=[]
    for fi,(a,b,mass) in enumerate(frames):
        # Resume only our atomic, input-bound completed outputs. Historical
        # Python outputs are left untouched and the native replay checks every
        # frame again, providing an additional complete differential comparison.
        record=out/f'frame_{fi:05}.json'
        if record.exists():
            r=json.loads(record.read_text());expected=[int(counts[fi][c]) for c in COLS]
            checked_frame(r,folder/'frames'/f'frame_{fi:05}',fi,expected,positive,n)
            if r['binary_sha256']!=sha(binary) or r['sources']!=hashes(folder/'frames'/f'frame_{fi:05}',positive):raise ValueError('changed resume inputs or verifier')
            total=[x+y for x,y in zip(total,r['count'])];weighted[0]+=mass*r['count'][1];weighted[1]+=mass*r['count'][2];completed+=1
        else:queue.append(fi)
    save(root/'FRAME_POPULATION.json',dict(frames=len(frames),sources={p.name:sha(p) for p in folder.glob('*.tsv') if not p.name.startswith('._')},scope_positive_only=positive))
    binary_hash=sha(binary);mutex=threading.Lock();cursor=0;stop=threading.Event()
    def lane(worker_id):
        nonlocal cursor,completed,total,weighted
        server=Native(binary,n)
        try:
            while not stop.is_set():
                with mutex:
                    if (root/'STOP').exists():stop.set();break
                    if cursor>=len(queue):break
                    fi=queue[cursor];cursor+=1
                if not Path('data').is_dir() or shutil.disk_usage(root).free<30*2**30:raise RuntimeError('T7 guard')
                a,b,mass=frames[fi];local=folder/'frames'/f'frame_{fi:05}';expected=[int(counts[fi][c]) for c in COLS]
                result=server.frame(local,a,b,positive);checked_frame(result,local,fi,expected,positive,n)
                sources=hashes(local,positive)
                previous=folder/'audit_frames'/f'frame_{fi:05}.json'
                if previous.exists():
                    old=json.loads(previous.read_text())
                    if old['count']!=result['count'] or old['sources']!=sources:raise ValueError('Python/native completed-frame discrepancy')
                result.update(frame=fi,sources=sources,binary_sha256=binary_hash,python_frame_sha256=sha(previous) if previous.exists() else None)
                save(record:=out/f'frame_{fi:05}.json',result)
                with mutex:
                    total=[x+y for x,y in zip(total,result['count'])];weighted[0]+=mass*result['count'][1];weighted[1]+=mass*result['count'][2];completed+=1
                    elapsed=time.monotonic()-start
                    state(root,status='auditing',pid=os.getpid(),workers=workers,gpu=False,completed_frames=completed,total_frames=len(frames),elapsed_seconds=elapsed,
                          eta_seconds=elapsed*(len(frames)-completed)/completed if completed else None)
                    if completed%25==0 or completed==len(frames):print(f'NATIVE_AUDIT frames={completed}/{len(frames)} elapsed={elapsed:.1f}',flush=True)
        except BaseException:stop.set();raise
        finally:server.close()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(lane,i) for i in range(workers)]):future.result()
    if completed!=len(frames):return dict(exhaustive=False,completed_frames=completed,total_frames=len(frames),stopped=True)
    names=['enumerated','transitive_smooth_occurrences','both_geodesic_occurrences','nonidentity_occurrences','identity_occurrences','nonsmooth_orbit_occurrences','closed_occurrences','unresolved']
    if total!=[int(coverage[c]) for c in names]:raise ValueError('global coverage mismatch')
    if [*weighted,2*weighted[0]-weighted[1]]!=[int(coverage[c]) for c in ['labeled_left_smooth','labeled_both_smooth','labeled_smooth_tuples']]:raise ValueError('labeled mass mismatch')
    return dict(exhaustive=True,coverage=coverage,components={},independently_reenumerated_third_occurrences=total[0],literal_orbit_outcomes_checked=total[3],audit_workers=workers,
                sources={p.name:sha(p) for p in folder.glob('*.tsv') if not p.name.startswith('._')},frame_population_sha256=sha(root/'FRAME_POPULATION.json'),
                frame_reports={p.name:sha(p) for p in out.glob('*.json') if not p.name.startswith('._')},elapsed_seconds=time.monotonic()-start)

def worker(root,lock_fd):
    lock=os.fdopen(lock_fd,'a');plan=json.loads((root/'PLAN.json').read_text());helper=None
    stop=threading.Event();children=set();mutex=threading.Lock()
    def guided(lane_id):
        epoch=0
        while not stop.is_set() and not (root/'STOP').exists():
            while list(root.glob(f'guided_{lane_id}_{epoch:07}_*')):epoch+=1
            n=[11,12,13,14][epoch%4] if lane_id==0 else [12,14][epoch%2]
            lane='forest' if lane_id==0 else 'blocks';name=f'guided_{lane_id}_{epoch:07}_{n}_{lane}';folder=root/name
            with (root/f'{name}.log').open('x') as log:
                child=subprocess.Popen([str(root/'smooth_forest_fast'),'--d',str(n),'--lane',lane,'--frame-reuse','256',
                    '--seed',str(300000003+lane_id*1000003+epoch),'--seconds','60','--out',str(folder)],stdout=log,stderr=subprocess.STDOUT)
                with mutex:children.add(child)
                try:
                    while child.poll() is None:
                        if shutil.disk_usage(root).free<30*2**30:child.terminate();child.wait();raise RuntimeError('guided T7 guard')
                        if stop.is_set():break
                        time.sleep(.5)
                    # A requested orderly audit handoff lets the bounded job
                    # finish and independently verifies it, preserving its leads.
                    if stop.is_set():child.wait(timeout=120)
                    if child.returncode:raise RuntimeError('guided producer failed')
                finally:
                    with mutex:children.discard(child)
            report=audit_guided(folder);save(folder/'CHECKS.json',report);(folder/'VERIFIED').write_text(now()+'\n')
            print('GUIDED_VERIFIED',name,json.dumps(report['construction']),flush=True);epoch+=1
    try:
        for name,digest in plan['source_hashes'].items():
            if sha(root/'source'/name)!=digest:raise ValueError('frozen source changed')
        if sha(root/'smooth_negative_audit')!=plan['binary_sha256']:raise ValueError('binary changed')
        if sha(root/'smooth_forest_fast')!=plan['guided_binary_sha256']:raise ValueError('guided binary changed')
        predecessor=Path(plan['predecessor']);saved=json.loads((predecessor/'STATE.json').read_text())
        if saved['status'] not in ['stopped','complete']:raise ValueError('predecessor must be safely stopped or complete')
        if plan['producer_table_hashes']!={p.name:sha(p) for p in (predecessor/'degree_10').glob('*.tsv') if not p.name.startswith('._')}:raise ValueError('producer tables changed')
        if shutil.which('caffeinate'):helper=subprocess.Popen(['caffeinate','-i','-w',str(os.getpid())])
        with ThreadPoolExecutor(max_workers=2) as pool:
            lanes=[pool.submit(guided,i) for i in range(2)]
            try:checks=audit(root,predecessor/'degree_10',root/'smooth_negative_audit',plan['workers'])
            finally:stop.set()
            for future in lanes:future.result()
        save(root/'EXACT_CHECKS.json',checks)
        status='complete' if checks['exhaustive'] else 'stopped'
        save(root/'FINISHED.json',dict(status=status,ended=now(),audit_sha256=sha(root/'EXACT_CHECKS.json'),predecessor=plan['predecessor']))
        state(root,status=status,workers=0,gpu=False)
    except BaseException as e:
        stop.set()
        with mutex:
            for child in children:
                if child.poll() is None:child.terminate()
        save(root/'FAILURE.json',dict(error=repr(e),ended=now()));state(root,status='failed',workers=0,error=repr(e));raise
    finally:
        if helper is not None and helper.poll() is None:helper.terminate()
        lock.close()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--predecessor',type=Path,default=DEFAULT)
    p.add_argument('--binary',type=Path,default=PROJECT/'target/audit-fast/release/smooth_negative_audit');p.add_argument('--checks',type=Path)
    p.add_argument('--guided-checks',type=Path,default=Path('data/forest_fast_benchmark_20261012_final/RESULTS.json'))
    p.add_argument('--workers',type=int,default=4);p.add_argument('--worker',type=Path);p.add_argument('--lock-fd',type=int);p.add_argument('--resume',type=Path)
    args=p.parse_args()
    if args.worker:worker(args.worker,args.lock_fd);return
    if args.workers not in range(1,5):p.error('audit workers 1..4, plus two guided workers')
    if not args.run:print(json.dumps(dict(predecessor=str(args.predecessor),workers=args.workers,scope='complete independent negative-only degree10 audit; every rank and literal word',gpu=False),indent=2));return
    control=PROJECT/'.run-control';control.mkdir(exist_ok=True);lock=(control/'.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if args.resume:
        root=args.resume
        if (root/'FINISHED.json').exists() and json.loads((root/'FINISHED.json').read_text())['status']=='complete':raise ValueError('already complete')
        # STOP is deliberately not removed automatically: resume is explicit.
        if (root/'STOP').exists():raise ValueError('remove STOP explicitly before resuming')
    else:
        if args.checks is None:p.error('validated --checks required')
        checks=json.loads(args.checks.read_text())
        if checks['binary_sha256']!=sha(args.binary) or not checks['passed']:raise ValueError('unbound/failed validation')
        guided_binary=PROJECT/'target/forest-fast/release/smooth_forest_fast';guided_checks=json.loads(args.guided_checks.read_text())
        if guided_checks['binary_sha256']!=sha(guided_binary) or guided_checks['controls']!=66 or guided_checks['control_states']!=25192:raise ValueError('unbound guided controls')
        dest=Path('data')
        if not Path('data').is_dir() or shutil.disk_usage(dest).free<30*2**30:raise RuntimeError('T7 guard')
        root=dest/('smooth_audit_fast_'+dt.datetime.now(ZoneInfo('America/Chicago')).strftime('%Y%m%d_%H%M%S_%f'));root.mkdir();source=root/'source';source.mkdir();hashes_={}
        paths=[PROJECT/'Cargo.toml',PROJECT/'Cargo.lock',*sorted((PROJECT/'src').rglob('*.rs')),
               *[s for s in sorted((PROJECT/'scripts').rglob('*.py')) if not s.name.startswith('._')],
               PROJECT/'evidence/negative_audit_acceleration_20261012/PLAN.md',PROJECT/'evidence/negative_audit_acceleration_20261012/MATH.md']
        for path in paths:
            rel=path.relative_to(PROJECT);target=source/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target);hashes_[str(rel)]=sha(target)
        shutil.copy2(args.binary,root/'smooth_negative_audit');shutil.copy2(args.checks,root/'VALIDATION.json')
        shutil.copy2(guided_binary,root/'smooth_forest_fast');shutil.copy2(args.guided_checks,root/'GUIDED_VALIDATION.json')
        save(root/'PLAN.json',dict(started=now(),predecessor=str(args.predecessor),guided_predecessors=[str(args.predecessor)],workers=args.workers,gpu=False,
            scope='positive-genus identity-free all-passport degree10 negative-only audit, complete rank/witness replay; not new discovery',
            stopping='STOP after current frames; partial outputs never exclusions; input-bound frame resume',source_hashes=hashes_,binary_sha256=sha(root/'smooth_negative_audit'),validation_sha256=sha(root/'VALIDATION.json'),
            guided_binary_sha256=sha(root/'smooth_forest_fast'),producer_table_hashes={p.name:sha(p) for p in (args.predecessor/'degree_10').glob('*.tsv') if not p.name.startswith('._')}))
    with (root/f'WORKER_{time.time_ns()}.log').open('x') as log:
        child=subprocess.Popen([sys.executable,str(root/'source/scripts/smooth/run_negative_audit.py'),'--worker',str(root),'--lock-fd',str(lock.fileno())],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json',dict(pid=child.pid,launched=now()));lock.close();print('OUTPUT='+str(root));print('PID='+str(child.pid))

if __name__=='__main__':main()
