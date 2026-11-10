#!/usr/bin/env python3
"""Frozen optimized guided successor. Plan-only unless --queue/--run is given."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import plistlib
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from zoneinfo import ZoneInfo

from run_smooth_frontier import audit_guided, values
from run_smooth_degree9 import sha, save, state, now
from run_s18_linked_roots import state_key
from odd_order_spin import cycle_lengths, parse_seed
from audit_guided_native import audit_guided_native,known_keys,prime_known_cache
from recover_forest_interruption import excluded_ids
from verified_smooth_bank import VerifiedSmoothBank

PROJECT=Path(__file__).resolve().parents[2]
DEFAULT_PREDECESSOR=Path('data/smooth_frontier_20261012_123427_295724')
MIN_FREE_BYTES=10*2**30


def storage_available(root):
    return Path(root).is_dir() and shutil.disk_usage(root).free>=MIN_FREE_BYTES

def compiled_source_files(folder, expected):
    """Bind a clean, isolated build without stashing unrelated workspace edits.

    All Rust files, not merely the three changed kernels, must match the
    completed validation receipt. The frozen run then contains the actual
    compilation inputs alongside its independently bound binary.
    """
    folder=Path(folder).resolve()
    if not expected:raise ValueError('compiled source validation required')
    for rel in expected:
        path=Path(rel)
        if path.is_absolute() or '..' in path.parts or not rel.startswith('src/') or path.suffix!='.rs':
            raise ValueError('invalid compiled source path')
    actual={str(p.relative_to(folder)):p for p in (folder/'src').rglob('*.rs')}
    if set(actual)!=set(expected):raise ValueError('compiled source file set differs')
    if any(sha(actual[rel])!=digest for rel,digest in expected.items()):
        raise ValueError('compiled source changed after validation')
    return actual

def resolved_native_history(folder, checks):
    """An audited zero-unresolved guided job cannot add a pending replay.

    Replay journals still need reading to retire inherited leads. Old reports
    without both agreeing counters also keep the literal recovery path.
    """
    return (Path(folder).name.startswith('guided_')
            and (Path(folder)/'VERIFIED').is_file()
            and checks.get('construction',{}).get('unresolved') in ('0',0)
            and checks.get('independent_native',{}).get('unresolved')==0)

def replace_queue_guard(path,data):
    # Lifecycle pointer, not immutable evidence. A dead historical waiter must
    # not prevent registering its already-spawned authorized successor.
    temp=path.with_name(path.name+'.tmp')
    with temp.open('w') as stream:json.dump(data,stream,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    temp.replace(path)

def predecessor_finished(root, inherit_stopped=False):
    finished=json.loads((root/'FINISHED.json').read_text())
    if finished['status']=='complete':return finished
    if inherit_stopped and finished['status']=='stopped' and json.loads((root/'STATE.json').read_text())['status']=='stopped':
        # STOP is finalized only after bounded children and their audits exit.
        return finished
    raise RuntimeError('predecessor audit did not complete; successor not started')

def profile_runtime(predecessor):
    runtime=shutil.which('uv')
    if runtime:return runtime
    # launchd's minimal PATH may omit the already checked profiler executable.
    plan=predecessor/'PLAN.json'
    if plan.is_file():
        saved=json.loads(plan.read_text()).get('profile_uv')
        if saved and Path(saved).is_absolute() and Path(saved).is_file() and os.access(saved,os.X_OK):return saved
    return None

def design_job(worker_id,epoch):
    # Pilot: joint construction gave worse promotion yield, not an improvement.
    # Keep it exploratory one job in eight on the second discovery lane.
    construction='forest' if worker_id<2 else 'solved' if worker_id==2 else 'quotient'
    if worker_id==1 and epoch%8==0:construction='joint'
    degrees=[11,12,13,14,16] if construction=='forest' else [11,12,13,14] if construction=='joint' else [15,18,20,21,24,25] if construction=='solved' else [12,15,16,18,20,21,24,25]
    return degrees[(epoch+worker_id)%len(degrees)],construction,1 if construction in ['solved','quotient'] else 256

PRIMITIVE_TASKS={0:[('psl11',11),('psl2_16',17)],1:[('m11',11),('psl3_3',13)]}

def mechanism_job(worker_id,epoch,cursors,double_cursors=None):
    if worker_id<2:
        for group,n in PRIMITIVE_TASKS[worker_id]:
            if not cursors.get(f'{group}_{n}',{}).get('producer_complete',False):
                return n,'primitive',1,group
        return [11,13,17][epoch%3],'forest',256,None
    if worker_id==2:return [12,15,16,18,20,21,24,25][epoch%8],'mutate',1,None
    if worker_id==3:return [11,12,13,14,16][epoch%5],'forest',256,None
    if worker_id==4:
        if epoch%2==0 and not (double_cursors or {}).get('24',{}).get('producer_complete',False):return 24,'double_lift',1,None
        return [12,15,16,18,20,21,24,25][epoch%8],'quotient',1,None
    raise ValueError('sixth mechanism slot is the structural audit, not a producer')

def launchd_pid(label):
    result=subprocess.run(['/bin/launchctl','list',label],capture_output=True,text=True)
    match=re.search(r'"PID"\s*=\s*(\d+)\s*;',result.stdout)
    if result.returncode or match is None:raise RuntimeError('launchd service has no running PID')
    return int(match.group(1))

def novelty_job(worker_id,epoch,cursors):
    if worker_id==4:
        for group in ['asp16','agl16']:
            if not cursors.get(group+'_16',{}).get('producer_complete',False):return 16,'primitive',1,group
        return [11,13,16,17,18][epoch%5],'forest',256,None
    if worker_id in [0,1]:return [11,12,13,14,15,16,17,18][(epoch+worker_id)%8],'triangle',1,None
    if worker_id==3 and epoch%16==0:return [12,15,16,18,20,21,24,25][(epoch//16)%8],'quotient',1,None
    if worker_id in [2,3]:return [12,14,15,16,18,20,21,24,25][(epoch+worker_id)%9],'returns',1,None
    if worker_id==5:return [12,15,16,18,20,21,24,25][epoch%8],('quotient' if epoch%16==0 else 'returns'),1,None
    raise ValueError('six-slot novelty mission')

def triangle_target_job(worker_id,epoch,spec):
    """Exact d16 target plus exact-genus-three higher-degree index shards."""
    shards=spec['passport_minimum_shards']
    if worker_id not in range(3) or len(shards)!=3:raise ValueError('three disjoint triangle shards required')
    band=shards[worker_id]
    low=int(band['passport_min'])
    first=max(int(spec['minimum_degree']),int(band.get('minimum_degree',spec['minimum_degree'])),2*low+1)
    last=min(int(spec['maximum_degree']),int(band.get('maximum_degree',spec['maximum_degree'])))
    degrees=[d for d in range(first,last+1) if (d-1)//2>=low]
    if not degrees:raise ValueError(f'empty triangle passport shard {worker_id}')
    n=degrees[epoch%len(degrees)]
    band=dict(band)
    # If the final minimum index is at least lo, both first-pair indices
    # must also be at least lo. Sampling smaller pair indices can never enter
    # this shard, so reject them at construction rather than after the tuple.
    band.setdefault('pair_index_min',low)
    if band.get('exact_passport'):
        band['third_method']='middle_repair'
    else:
        band['third_method']=['random','joint','mixed','strong','centralizer'][(epoch+worker_id)%5]
    return n,'triangle',1,None,band

def affine_interval_cover(proofs,total):
    """Aggregate already independently checked, bound intervals; no gap claims."""
    cursor=0;qualifying=negative=survivors=0;unresolved=False
    for p in sorted(proofs,key=lambda x:x['start']):
        if not p['passed'] or p['start']!=cursor or p['next']<cursor or p['total']!=total:raise ValueError('affine coverage gap/overlap/binding')
        cursor=p['next'];c=p['counts'];qualifying+=c['qualifying'];negative+=c['literal_negative'];survivors+=c['accounted_survivor']
        unresolved|=not p['no_unresolved']
    if cursor>total:raise ValueError('affine coverage past population')
    return dict(next=cursor,total=total,independent_complete_population=cursor==total,
                no_unresolved=not unresolved,negative_census=cursor==total and not unresolved and qualifying==negative,
                counts=dict(qualifying=qualifying,literal_negative=negative,accounted_survivor=survivors))

def spawn_launchd(root,command,label,timeout=30):
    # One-off user service: no login installation, recurrence or auto-restart.
    # launchd owns the actual controller, not a short-lived launcher process.
    # Internal startup log; the Python worker opens its own T7 log. Use an
    # explicit non-KeepAlive plist: legacy `submit` implicitly restarts failures.
    folder=Path(tempfile.mkdtemp(prefix='hurwitz-launchd-'));log=folder/'bootstrap.log'
    config=folder/'service.plist';domain=f'gui/{os.getuid()}'
    config.write_bytes(plistlib.dumps(dict(Label=label,ProgramArguments=list(map(str,command)),
                     RunAtLoad=True,KeepAlive=False,StandardOutPath=str(log),StandardErrorPath=str(log))))
    shutil.copyfile(config,root/'LAUNCHD_CONFIG.plist')
    save(root/'LAUNCHD_BOOTSTRAP.json',dict(label=label,domain=domain,log=str(log),
         config=str(config),config_sha256=sha(config),keepalive=False))
    subprocess.run(['/bin/launchctl','bootstrap',domain,str(config)],check=True,capture_output=True,text=True)
    deadline=time.monotonic()+timeout
    while time.monotonic()<deadline:
        try:
            pid=launchd_pid(label);record=json.loads((root/'STATE.json').read_text())
            if record['pid']==pid and record['status']=='queued':return pid
        except (RuntimeError,OSError,ValueError,KeyError):pass
        time.sleep(.02)
    # Parent still owns the campaign lock: no producer has been admitted.
    (root/'STOP').write_text('Failed to bind one-off launchd service to queued controller.\n')
    subprocess.run(['/bin/launchctl','bootout',domain+'/'+label],capture_output=True,text=True)
    raise RuntimeError('launchd/controller registration did not agree; new service stopped')

def activate_prepared(root):
    root=root.resolve();plan=json.loads((root/'PLAN.json').read_text())
    if plan.get('supervisor')!='launchd' or not plan.get('design_search'):raise ValueError('prepared managed design search required')
    if (root/'LAUNCH.json').exists():
        # A delayed task replay cannot restart a previously activated run.
        print('ALREADY_ACTIVATED='+str(root));return
    if json.loads((root/'STATE.json').read_text())['status']!='prepared':raise ValueError('not prepared')
    predecessor_finished(Path(plan['predecessor']),plan.get('inherit_stopped',False))
    for name,digest in plan['source_hashes'].items():
        if sha(root/'source'/name)!=digest:raise ValueError('prepared source changed')
    for name,field in [('smooth_forest_fast','binary_sha256'),('smooth_negative_audit','audit_binary_sha256'),
                       ('KNOWN_STATES.tsv','known_bank_sha256'),('VALIDATION.json','validation_sha256'),
                       ('CONTROLS.tsv','controls_sha256')]:
        if sha(root/name)!=plan[field]:raise ValueError('prepared artifact changed: '+name)
    for name,field in [('MUTATION_SEEDS.tsv','mutation_seeds_sha256'),('PROFILE_INPUT.json','profile_input_sha256'),('PRIMITIVE_EXCLUSIONS.json','primitive_exclusions_sha256'),('DOUBLE_BASES.tsv','double_bases_sha256')]:
        if plan.get('mechanism_search') and sha(root/name)!=plan[field]:raise ValueError('prepared mechanism input changed: '+name)
    if plan.get('extra_leads_sha256') and sha(root/'EXTRA_LEADS.json')!=plan['extra_leads_sha256']:raise ValueError('prepared saved leads changed')
    if plan.get('novelty_search') and sha(root/'NOVELTY_INPUT.json')!=plan['novelty_input_sha256']:raise ValueError('prepared novelty inputs changed')
    control=Path(plan['campaign_lock']).parent
    queue_lock=(control/'.forest-fast-queue.lock').open('a');fcntl.flock(queue_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    queue_guard=control/'FOREST_FAST_QUEUED.json'
    if queue_guard.exists():
        previous=json.loads(queue_guard.read_text())
        try:os.kill(previous['pid'],0)
        except ProcessLookupError:pass
        else:raise RuntimeError('another successor is queued/running')
    lock=Path(plan['campaign_lock']).open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if not storage_available(root):raise RuntimeError('T7 guard')
    if (root/'STOP').exists():raise RuntimeError('prepared run cancelled by STOP')
    label='local.hurwitz.'+root.name
    command=[sys.executable,str(root/'source/scripts/smooth/run_forest_fast.py'),'--worker',str(root)]
    pid=spawn_launchd(root,command,label)
    replace_queue_guard(queue_guard,dict(pid=pid,output=str(root),queued=now(),workers=plan['workers'],launchd_label=label))
    save(root/'LAUNCH.json',dict(pid=pid,queued=now(),supervisor='launchd',launchd_label=label))
    lock.close();queue_lock.close();print('OUTPUT='+str(root));print('PID='+str(pid))

def worker(root,lock_fd):
    plan=json.loads((root/'PLAN.json').read_text())
    lock=os.fdopen(lock_fd,'a') if lock_fd is not None else Path(plan['campaign_lock']).open('a')
    if (root/'source').is_dir():os.chdir(root/'source')
    if plan.get('supervisor')=='launchd':
        with (root/'WORKER.log').open('x') as log:
            os.dup2(log.fileno(),sys.stdout.fileno());os.dup2(log.fileno(),sys.stderr.fileno())
    state(root,status='queued',pid=os.getpid(),workers=0,waiting_for=plan['predecessor'])
    # Block in the OS, consuming no search core. The predecessor keeps its lock,
    # audit and frozen configuration. No competing pool is created.
    fcntl.flock(lock,fcntl.LOCK_EX)
    stop=threading.Event();children=set();mutex=threading.Lock();helper=None
    def run(command,log):
        with log.open('x')as f:
            child=subprocess.Popen(list(map(str,command)),stdout=f,stderr=subprocess.STDOUT)
            with mutex:children.add(child)
            try:
                while child.poll()is None:
                    if (root/'STOP').exists():stop.set()
                    if not storage_available(root):
                        child.terminate();child.wait();raise RuntimeError('T7 space/mount guard; partial seeds preserved')
                    time.sleep(.5)
                if child.returncode:raise RuntimeError(f'child failed exit={child.returncode}')
            finally:
                with mutex:children.discard(child)
    try:
        predecessor=Path(plan['predecessor'])
        finished=predecessor_finished(predecessor,plan.get('inherit_stopped',False))
        if not storage_available(root):raise RuntimeError('T7 unavailable/space guard')
        for name,digest in plan['source_hashes'].items():
            if sha(root/'source'/name)!=digest:raise ValueError('frozen source changed')
        if sha(root/'smooth_forest_fast')!=plan['binary_sha256']:raise ValueError('frozen binary changed')
        native=root/'smooth_negative_audit' if plan.get('audit_binary_sha256') else None
        if native is not None and sha(native)!=plan['audit_binary_sha256']:raise ValueError('frozen auditor changed')
        if plan.get('design_search') and sha(root/'KNOWN_STATES.tsv')!=plan['known_bank_sha256']:raise ValueError('frozen known bank changed')
        if plan.get('extra_leads_sha256') and sha(root/'EXTRA_LEADS.json')!=plan['extra_leads_sha256']:raise ValueError('frozen extra leads changed')
        if plan.get('novelty_search') and sha(root/'NOVELTY_INPUT.json')!=plan['novelty_input_sha256']:raise ValueError('frozen novelty inputs changed')
        if plan.get('mechanism_search'):
            for name,field in [('MUTATION_SEEDS.tsv','mutation_seeds_sha256'),('PROFILE_INPUT.json','profile_input_sha256'),('PRIMITIVE_EXCLUSIONS.json','primitive_exclusions_sha256'),('DOUBLE_BASES.tsv','double_bases_sha256')]:
                if sha(root/name)!=plan[field]:raise ValueError('mechanism input changed')
        if (root/'STOP').exists():state(root,status='stopped',workers=0);return
        if shutil.which('caffeinate'):helper=subprocess.Popen(['caffeinate','-i','-w',str(os.getpid())])
        # Read only immutable VERIFIED jobs after the predecessor releases its
        # lock. Deduplicate saved capped seeds by exact S_n canonical state.
        pending=[];seen=set();resolved_replays=set();historical_components=set();skipped_history=0
        history=[predecessor]
        parent_plan=json.loads((predecessor/'PLAN.json').read_text())
        history.extend(Path(p) for p in parent_plan.get('guided_predecessors',[]))
        for old_root in history:
            if json.loads((old_root/'STATE.json').read_text())['status'] not in ['complete','stopped']:raise ValueError('guided history still active')
        for folder in sorted(f for old_root in history for pattern in ['guided_*','replay_*'] for f in old_root.glob(pattern)):
            if not folder.is_dir():continue
            verified=(folder/'VERIFIED').exists()
            if verified:
                checks=json.loads((folder/'CHECKS.json').read_text());n=int(checks['construction']['degree'])
                # These identities affect the discovery counter only, never a
                # producer's positive bank or acceptance without full closure.
                for item in checks.get('closed_identities',[]):
                    historical_components.add((int(item['degree']),tuple(tuple(p) for p in item['key'])))
                if resolved_native_history(folder,checks):
                    skipped_history+=1;continue
            else:
                # Safely stopped predecessor may have a lead flushed before
                # an interrupted BFS. Replay such tuples, never call them done.
                raw=(folder/'LEADS.tsv')
                if not raw.exists():continue
                n=None
                log=folder.with_name(folder.name+'.log')
                # Degree is recovered from the persisted tuple support AND
                # the old job name, which explicitly includes n for this root.
                parts=folder.name.split('_')
                if len(parts)>=5 and parts[-2].isdigit():n=int(parts[-2])
                elif len(parts)==3 and parts[0]=='replay' and parts[-1].isdigit():n=int(parts[-1])
                elif len(parts)==2 and parts[0]=='replay' and parts[-1].isdigit():
                    inherited=json.loads((folder.parent/'INHERITED_LEADS.json').read_text())
                    n=int(inherited[int(parts[-1])]['degree'])
                if n is None:raise ValueError('partial guided job lacks recoverable degree')
            import csv
            if verified:
                stream=(folder/'RESOLUTIONS.tsv').open()
                lead_rows=csv.DictReader(stream,delimiter='\t')
            else:
                stream=None
                exclusions=excluded_ids(folder.parent,folder)
                stream=raw.open()
                def interrupted_rows(source,retired):
                    for line in source:
                        if not line.endswith('\n'):
                            if plan.get('novelty_search'):break
                            raise ValueError('incomplete durable lead retained; manual recovery needed')
                        i,t=line.rstrip('\n').split('\t',1)
                        if i not in retired:yield dict(candidate=i,seed=t,status='interrupted')
                lead_rows=interrupted_rows(stream,exclusions)
            for row in lead_rows:
                if verified and folder.name.startswith('replay_') and row['status'] in ['nonsmooth','closed_smooth']:
                    resolved_replays.add((n,state_key(parse_seed(row['seed'],n))))
                if row['status']not in ['timeout','capped','interrupted']:continue
                t=parse_seed(row['seed'],n);v=values(t)
                if not v['smooth']or v['genus']<=0 or any(p==tuple(range(n)) for p in t):raise ValueError('out-of-scope recovered lead')
                key=(n,state_key(parse_seed(row['seed'],n)))
                if key in seen:continue
                seen.add(key);pending.append(dict(degree=n,seed=row['seed'],source=str(folder),candidate=row['candidate']))
            if stream is not None:stream.close()
        # A STOP can happen before an inherited queue is fully issued. Carry
        # its unissued rows too, retiring only independently VERIFIED replays.
        inherited_path=predecessor/'INHERITED_LEADS.json'
        if inherited_path.exists():
            for item in json.loads(inherited_path.read_text()):
                n=int(item['degree']);key=(n,state_key(parse_seed(item['seed'],n)))
                if key not in seen and key not in resolved_replays:
                    seen.add(key);pending.append(item)
        if plan.get('extra_leads_sha256'):
            for item in json.loads((root/'EXTRA_LEADS.json').read_text()):
                n=int(item['degree']);t=parse_seed(item['seed'],n);v=values(t)
                if n not in range(9,26) or not v['smooth'] or v['genus']<=0 or any(p==tuple(range(n)) for p in t):raise ValueError('extra lead scope')
                key=(n,state_key(t))
                if key in seen or key in resolved_replays:continue
                seen.add(key);pending.append(item)
        triangle_target=plan.get('triangle_target')
        if triangle_target:
            before=len(pending);kept=[];rejected=dict(degree=0,genus=0,passport=0,high_fixed_point=0)
            for item in pending:
                n=int(item['degree'])
                if not int(triangle_target['minimum_degree'])<=n<=int(triangle_target['maximum_degree']):
                    rejected['degree']+=1;continue
                t=parse_seed(item['seed'],n);ix=sorted(n-len(cycle_lengths(p)) for p in t)
                if not values(t)['smooth'] or values(t)['genus']<int(triangle_target['minimum_genus']):
                    rejected['genus']+=1;continue
                if ix[1]+ix[2]>n-1 or ix[0]+ix[3]<n:
                    rejected['passport']+=1;continue
                high=ix[3]
                if not any(n-len(cycle_lengths(p))==high and 1 not in cycle_lengths(p) for p in t):
                    rejected['high_fixed_point']+=1;continue
                kept.append(item)
            pending=kept
            save(root/'TRIANGLE_REPLAY_FILTER.json',dict(inherited_leads=before,
                triangle_replays_retained=len(kept),retained_degrees=sorted({int(x['degree']) for x in kept}),
                rejected=rejected,scope='Only saved smooth unresolved seeds in the fresh degree window that meet the checked triangle budget, genus floor, and fixed-point-free high-entry condition are replayed. All other raw sources remain preserved in predecessor history.'))
        novelty=json.loads((root/'NOVELTY_INPUT.json').read_text()) if plan.get('novelty_search') else None
        if novelty:
            for item in novelty['negative_frontier']:
                folder=Path(item['folder']);checks_path=folder/'CHECKS.json'
                if sha(checks_path)!=item['checks_sha256']:raise ValueError('frontier proof changed')
                report=json.loads(checks_path.read_text())
                if report['independent_native']['closed'] or report['independent_native']['unresolved']:raise ValueError('frontier is not negative-only')
                for name,digest in report['sources'].items():
                    if sha(folder/name)!=digest:raise ValueError('frontier journal changed')
                with (folder/'LEADS.tsv').open() as stream:
                    for line in stream:
                        _,seed=line.rstrip('\n').split('\t',1);n=item['degree']
                        resolved_replays.add((n,state_key(parse_seed(seed,n))))
            pending=[e for e in pending if (e['degree'],state_key(parse_seed(e['seed'],e['degree']))) not in resolved_replays]
        save(root/'INHERITED_LEADS.json',pending);save(root/'PREDECESSOR_COMPLETE.json',finished)
        save(root/'HISTORY_RECOVERY.json',dict(zero_unresolved_native_journals_skipped=skipped_history,
            inherited_component_identities=len(historical_components),pending_leads=len(pending),
            scope='Metadata optimization only; all nonzero/legacy/replay/partial journals retain literal recovery.'))
        cursor=0
        prior_epochs=finished.get('epochs',[]) if plan.get('inherit_stopped') else []
        if any(not isinstance(e,int) or e<0 for e in prior_epochs):raise ValueError('invalid predecessor epoch checkpoint')
        if len(prior_epochs)==plan['workers']:
            epoch=prior_epochs[:]
            epoch_handoff=dict(previous_workers=len(prior_epochs),new_workers=plan['workers'],
                previous_epochs=prior_epochs,epochs=epoch,policy='same worker count; continue each lane after its finalized predecessor epoch')
        elif prior_epochs:
            next_epoch=max(prior_epochs)+1
            epoch=[(prior_epochs[i]+1 if i<len(prior_epochs) else next_epoch) for i in range(plan['workers'])]
            epoch_handoff=dict(previous_workers=len(prior_epochs),new_workers=plan['workers'],
                previous_epochs=prior_epochs,epochs=epoch,policy='worker topology changed; continue corresponding old lane ids after their own final epochs and seed new lanes beyond every old epoch')
        else:
            epoch=[0]*plan['workers']
            epoch_handoff=dict(previous_workers=0,new_workers=plan['workers'],previous_epochs=[],epochs=epoch,
                policy='no predecessor epoch checkpoint; initialize lanes at zero')
        save(root/'EPOCH_HANDOFF.json',epoch_handoff)
        mutex_queue=threading.Lock()
        progress=dict(completed_jobs=0,closed=0,unresolved=0)
        discoveries={}
        if plan.get('mechanism_search'):save(root/'KNOWN_CACHE.json',prime_known_cache(root/'KNOWN_STATES.tsv',native))
        all_known=set(known_keys(root/'KNOWN_STATES.tsv')) if plan.get('design_search') else set()
        all_known.update(historical_components)
        profile_queue=json.loads((root/'PROFILE_INPUT.json').read_text()) if plan.get('mechanism_search') else None
        if profile_queue:
            all_known.update((e['degree'],tuple(tuple(p) for p in e['component_identity_key']))
                             for e in profile_queue['components'] if 'component_identity_key' in e)
        primitive_cursors={};double_cursors={};profiles_done=novelty['profile_start'] if novelty else 0;profile_summary={};primitive_audited=set();primitive_jobs={};affine_proofs={}
        if plan.get('mechanism_search'):
            exclusions=json.loads((root/'PRIMITIVE_EXCLUSIONS.json').read_text())
            for key,e in exclusions.items():
                proof=Path(e['path'])
                if sha(proof)!=e['sha256']:raise ValueError('primitive exclusion proof changed')
                report=json.loads(proof.read_text())
                if not report['passed'] or not report['independent_complete_population'] or not report['no_unresolved'] or report['counts'].get('qualifying',0)!=report['counts'].get('literal_negative',0):raise ValueError('not an exact negative family')
                if any(sha(proof.parent/name)!=digest for name,digest in report['sources'].items()):raise ValueError('primitive exclusion population source changed')
                primitive_cursors[key]=dict(next=report['total'],total=report['total'],producer_complete=True)
                primitive_audited.add(key)
        if novelty:
            for key,item in novelty.get('primitive_prefixes',{}).items():
                proof=Path(item['path']);report=json.loads(proof.read_text())
                if sha(proof)!=item['sha256'] or not report['passed'] or report['start']!=0 or not report['no_unresolved']:raise ValueError('affine prefix binding')
                if any(sha(proof.parent/name)!=digest for name,digest in report['sources'].items()):raise ValueError('affine prefix source changed')
                if key in primitive_cursors:raise ValueError('prefix overlaps retired population')
                primitive_cursors[key]=dict(next=report['next'],total=report['total'],producer_complete=report['next']==report['total'],
                    group_sha256=report['sources']['GROUP.tsv'],frames_sha256=report['sources']['PAIR_FRAMES.tsv'])
                affine_proofs[key]=[(proof,report)]
        shared_bank=VerifiedSmoothBank(root/'VERIFIED_BANKS',root/'KNOWN_STATES.tsv',plan['known_bank_sha256']) if plan.get('design_search') else None
        end=time.monotonic()+plan['hours']*3600 if plan['hours']>0 else float('inf')
        state(root,status='searching',pid=os.getpid(),workers=plan['workers'],gpu=False,pending_leads=len(pending))
        def lane(worker_id):
            nonlocal cursor,profiles_done
            banks={}
            while not stop.is_set()and not(root/'STOP').exists()and time.monotonic()<end:
                replay=None;bank_binding=None
                with mutex_queue:
                    can_replay=worker_id in [2,3] if novelty else not(plan.get('mechanism_search') and worker_id<2)
                    if cursor<len(pending) and can_replay:replay=(cursor,pending[cursor]);cursor+=1
                if replay:
                    i,item=replay;name=f'replay_{i:06}';n=item['degree'];seed_file=root/f'{name}.txt'
                    seed_file.write_text(item['seed']+'\n')
                    command=[root/'smooth_forest_fast','--d',n,'--seconds',plan['replay_seconds'],'--replay',seed_file,'--out',root/name]
                else:
                    e=epoch[worker_id];n=[11,12,13,14][(e+worker_id)%4]
                    if plan.get('mechanism_search') and worker_id==5 and (not novelty or profiles_done<len(profile_queue['components'])):
                        with mutex_queue:
                            audit_task=next(((key,folders[:]) for key,folders in primitive_jobs.items()
                                if key not in primitive_audited and primitive_cursors[key]['producer_complete']),None)
                        if audit_task:
                            key,folders=audit_task
                            for job in folders:
                                name='family_'+job.name;run([sys.executable,root/'source/scripts/smooth/audit_primitive_family.py',
                                    '--folder',job,'--out',job/'EXACT_POPULATION.json'],root/f'{name}.log')
                            with mutex_queue:primitive_audited.add(key)
                            print('PRIMITIVE_FAMILY_AUDITED',key,flush=True);continue
                        if profiles_done>=len(profile_queue['components']):
                            state(root/f'lane_{worker_id}',status='profile_complete',components=profiles_done)
                            stop.wait(1);continue
                        name=f'profile_{profiles_done:06}';folder=root/name
                        batch=min(8,len(profile_queue['components'])-profiles_done) if novelty else 1
                        command=[plan['profile_uv'],'run','--isolated','--no-project','--python','3.13','--with','sympy==1.14.0','python',
                            root/'source/scripts/smooth/profile_mechanism_breaker.py','--queue',root/'PROFILE_INPUT.json',
                            '--start',profiles_done,'--limit',batch,'--out',folder]
                        if novelty:command+=['--all-block-systems','--batch-only']
                        run(command,root/f'{name}.log')
                        result=json.loads((folder/'PROFILES.json').read_text())
                        if result['queue_sha256']!=plan['profile_input_sha256'] or result['start']!=profiles_done or result['next']!=profiles_done+batch:raise ValueError('profile cursor/input binding')
                        for value in result['components']:
                            kind=value['mechanism'];profile_summary[kind]=profile_summary.get(kind,0)+1
                            projection=value.get('marked_projection',{})
                            if projection.get('status')=='checked':
                                key='faithful_projection' if projection['faithful'] else 'nonfaithful_projection';profile_summary[key]=profile_summary.get(key,0)+1
                                if not projection['deterministic_factor']:profile_summary['projection_not_factor']=profile_summary.get('projection_not_factor',0)+1
                        profiles_done+=batch;epoch[worker_id]+=1
                        temp=root/'PROFILE_SUMMARY.tmp';temp.write_text(json.dumps(dict(checked=profiles_done,total=len(profile_queue['components']),counts=profile_summary,inherited=novelty['profile_start'] if novelty else 0),indent=2)+'\n');temp.replace(root/'PROFILE_SUMMARY.json')
                        state(root/f'lane_{worker_id}',status='profiling',latest=name,components=profiles_done)
                        print('PROFILE_CHECKED',profiles_done,flush=True);continue
                    group=None
                    triangle_band=None
                    if triangle_target:
                        n,construction,reuse,group,triangle_band=triangle_target_job(worker_id,e,triangle_target)
                    elif novelty:
                        n,construction,reuse,group=novelty_job(worker_id,e,primitive_cursors)
                    elif plan.get('mechanism_search'):
                        n,construction,reuse,group=mechanism_job(worker_id,e,primitive_cursors,double_cursors)
                    elif plan.get('design_search'):
                        n,construction,reuse=design_job(worker_id,e)
                    else:
                        construction='blocks'if worker_id==plan['workers']-1 else'forest'
                        if construction=='blocks':n=[12,14][e%2]
                        reuse=[1,64,256,1024][worker_id%4]
                    name=f'guided_{worker_id}_{e:07}_{n}_{construction.replace("_","-")}'
                    depth=0 if novelty and construction=='primitive' else plan.get('block_deep_level',0) if construction=='blocks' else plan.get('deep_level',3)
                    command=[root/'smooth_forest_fast','--d',n,'--lane',construction,'--frame-reuse',reuse,'--deep-level',depth,
                             '--seed',200000003+1000003*worker_id+e,'--seconds',plan['job_seconds'],'--out',root/name]
                    if novelty or triangle_target:command+=['--full-closure-budget']
                    if triangle_target:
                        command+=['--triangle-min-genus',triangle_target['minimum_genus'],
                                  '--triangle-low-min',triangle_band.get('pair_index_min',2),
                                  '--triangle-low-max',triangle_band.get('pair_index_max',63),
                                  '--triangle-passport-min',triangle_band['passport_min'],
                                  '--triangle-passport-max',triangle_band['passport_max'],
                                  '--triangle-third-method',triangle_band['third_method'],
                                  '--triangle-strength-divisor',triangle_band['strong_divisor'],
                                  '--closure-seconds',triangle_target['closure_seconds'],
                                  '--cap',triangle_target['closure_cap']]
                        if triangle_band.get('exact_genus'):command+=['--target-genus',triangle_band['exact_genus']]
                        if triangle_band.get('exact_passport'):
                            command+=['--triangle-exact-passport',','.join(map(str,triangle_band['exact_passport']))]
                    if construction=='blocks':command+=['--no-orbit-cache']
                    if group:
                        key=f'{group}_{n}';command+=['--primitive-group',group,'--primitive-cursor',primitive_cursors.get(key,{}).get('next',0)]
                    if construction=='mutate':command+=['--mutation-seeds',root/'MUTATION_SEEDS.tsv']
                    if construction=='double_lift':command+=['--double-bases',root/'DOUBLE_BASES.tsv','--primitive-cursor',double_cursors.get(str(n),{}).get('next',0)]
                    if plan.get('mechanism_search') and worker_id<2 and construction=='forest':command+=['--target-genus',1]
                    if plan.get('design_search'):
                        bank_binding=shared_bank.snapshot(n)
                        command+=['--compact-frames','--export-negative-cache','--known-bank',bank_binding[0]]
                        if n in banks:command+=['--negative-bank',banks[n]]
                    epoch[worker_id]+=1
                if plan.get('design_search') and sha(root/'KNOWN_STATES.tsv')!=plan['known_bank_sha256']:raise ValueError('frozen known bank changed')
                if bank_binding and sha(bank_binding[0])!=bank_binding[1]:raise ValueError('immutable shared bank changed')
                run(command,root/f'{name}.log')
                checks=audit_guided_native(root/name,native,collect_records=bool(shared_bank),
                    native_components=plan.get('native_components',False)) if native is not None else audit_guided(root/name)
                if (root/name/'PRIMITIVE.tsv').exists():
                    from run_smooth_degree9 import rows
                    info=rows(root/name/'PRIMITIVE.tsv')[0];key=info['group']+'_'+info['degree']
                    start=int(info['start']);next_cursor=int(info['next']);total=int(info['total'])
                    if start!=primitive_cursors.get(key,{}).get('next',0) or next_cursor-start!=int(checks['construction']['constructed']) or not 0<=start<=next_cursor<=total:raise ValueError('primitive finite cursor mass')
                    binding=dict(group_sha256=sha(root/name/'GROUP.tsv'),frames_sha256=sha(root/name/'PAIR_FRAMES.tsv'))
                    previous=primitive_cursors.get(key,{})
                    if previous and any(previous[field]!=value for field,value in binding.items()):raise ValueError('primitive population changed between slices')
                    with mutex_queue:primitive_cursors[key]=dict(next=next_cursor,total=total,producer_complete=info['producer_complete']=='true',**binding)
                    checks['primitive_population']=dict(info,**binding,independent_full_population_coverage=False)
                    if novelty:
                        proof=json.loads(subprocess.check_output([str(native),'--primitive-population',str(root/name)],text=True))
                        proof['sources']={f:sha(root/name/f) for f in ['GROUP.tsv','PAIR_FRAMES.tsv','PRIMITIVE.tsv','LEADS.tsv','RESOLUTIONS.tsv','COUNTS.tsv']}
                        save(root/name/'EXACT_POPULATION.json',proof)
                        checks['primitive_population']['independent_interval_coverage']=True
                        with mutex_queue:
                            primitive_audited.add(key)
                            affine_proofs.setdefault(key,[]).append((root/name/'EXACT_POPULATION.json',proof))
                            aggregate=affine_interval_cover([p for _,p in affine_proofs[key]],total)
                            aggregate.update(group=info['group'],degree=n,passed=True,
                                proofs=[dict(path=str(path),sha256=sha(path)) for path,_ in affine_proofs[key]],
                                scope='Complete geodesic first-pair frames in this named affine action; not a whole-degree census.')
                            save(root/f'AFFINE_COVERAGE_{key}.json',aggregate)
                published_records=checks.pop('_closed_records',[])
                if bank_binding and checks.get('input_banks',{}).get('--known-bank',{}).get('sha256')!=bank_binding[1]:raise ValueError('known bank audit binding')
                save(root/name/'CHECKS.json',checks)
                (root/name/'VERIFIED').write_text(now()+'\n')
                if (root/name/'DOUBLE_LIFT.tsv').exists():
                    from run_smooth_degree9 import rows
                    info=rows(root/name/'DOUBLE_LIFT.tsv')[0]
                    if int(info['start'])!=double_cursors.get(str(n),{}).get('next',0) or int(info['next'])-int(info['start'])!=int(checks['construction']['constructed']):raise ValueError('double-lift finite cursor mass')
                    double_cursors[str(n)]=dict(next=int(info['next']),total=int(info['total']),producer_complete=info['producer_complete']=='true')
                    temp=root/'DOUBLE_PROGRESS.tmp';temp.write_text(json.dumps(double_cursors,indent=2)+'\n');temp.replace(root/'DOUBLE_PROGRESS.json')
                if checks.get('primitive_population'):
                    with mutex_queue:
                        primitive_jobs.setdefault(key,[]).append(root/name)
                        temp=root/'PRIMITIVE_PROGRESS.tmp';temp.write_text(json.dumps(primitive_cursors,indent=2)+'\n');temp.replace(root/'PRIMITIVE_PROGRESS.json')
                if shared_bank and published_records:
                    publication=shared_bank.publish(root/name,checks,published_records)
                    save(root/name/'BANK_PUBLICATION.json',publication)
                if plan.get('design_search') and checks.get('exported_negative_bank'):
                    banks[n]=root/name/'NEGATIVE_CACHE.tsv'
                state(root/f'lane_{worker_id}',status='searching',latest=name,jobs=epoch[worker_id])
                with mutex_queue:
                    progress['completed_jobs']+=1
                    progress['closed']+=int(checks['construction']['closed'])
                    progress['unresolved']+=int(checks['construction']['unresolved'])
                    for item in checks.get('closed_identities',[]):
                        exact_key=tuple(tuple(p) for p in item['key'])
                        if (item['degree'],exact_key) in all_known:continue
                        identity=json.dumps([item['degree'],item['key']],separators=(',',':'))
                        discoveries.setdefault(identity,dict(**item,job=name))
                    progress['distinct_closed_components']=len(discoveries)
                    # Atomic lifecycle index, not the immutable certificates.
                    temp=root/'NEW_COMPONENTS.tmp';temp.write_text(json.dumps(list(discoveries.values()),indent=2)+'\n')
                    temp.replace(root/'NEW_COMPONENTS.json')
                    state(root,status='searching',pid=os.getpid(),workers=plan['workers'],gpu=False,
                          pending_leads=max(0,len(pending)-cursor),epochs=epoch[:],**progress)
                print('GUIDED_VERIFIED',name,json.dumps(checks['construction']),flush=True)
        for i in range(plan['workers']):(root/f'lane_{i}').mkdir()
        with ThreadPoolExecutor(max_workers=plan['workers'])as pool:
            futures=[pool.submit(lane,i)for i in range(plan['workers'])]
            for future in as_completed(futures):
                try:future.result()
                except BaseException:
                    stop.set();raise
        save(root/'FINISHED.json',dict(status='stopped'if(root/'STOP').exists()else'budget_complete',ended=now(),epochs=epoch))
        state(root,status='stopped'if(root/'STOP').exists()else'budget_complete',workers=0)
    except BaseException as error:
        stop.set()
        with mutex:
            for child in children:
                if child.poll()is None:child.terminate()
        save(root/'FAILURE.json',dict(error=repr(error),ended=now()));state(root,status='failed',workers=0,error=repr(error))
        raise
    finally:
        if helper is not None and helper.poll()is None:helper.terminate()
        lock.close()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--queue',action='store_true');p.add_argument('--run',action='store_true')
    p.add_argument('--prepare',action='store_true',help='freeze checked inputs without activating; for an explicit later launch')
    p.add_argument('--activate',type=Path,help='idempotently activate a prepared frozen managed run')
    p.add_argument('--predecessor',type=Path,default=DEFAULT_PREDECESSOR)
    p.add_argument('--binary',type=Path,default=PROJECT/'target/forest-fast/release/smooth_forest_fast')
    p.add_argument('--compiled-source',type=Path,help='isolated Rust source tree bound by completed validation; preserves unrelated dirty workspace files')
    p.add_argument('--checks',type=Path);p.add_argument('--workers',type=int,default=6)
    p.add_argument('--extra-leads',type=Path,help='independently checked saved pilot caps; replay before new proposals')
    p.add_argument('--audit-binary',type=Path)
    p.add_argument('--inherit-stopped',action='store_true')
    p.add_argument('--launchd',action='store_true',help='one-off macOS-owned controller; survives launcher exit')
    p.add_argument('--deep-level',type=int,default=3,choices=range(4))
    p.add_argument('--block-deep-level',type=int,default=0,choices=range(4))
    p.add_argument('--design-search',action='store_true',help='solved composite-degree lifts, joint forests, verified banks, compact frames')
    p.add_argument('--native-components',action='store_true',help='separately validated complete independent Rust component audit')
    p.add_argument('--mechanism-search',action='store_true',help='two primitive lanes, mutation, baseline, quotient and one profile slot')
    p.add_argument('--novelty-input',type=Path,help='validated six-slot structural novelty mission')
    p.add_argument('--triangle-target',action='store_true',help='dedicated three-worker, disjoint TRIANGLE index-passport search')
    p.add_argument('--mutation-seeds',type=Path);p.add_argument('--profile-input',type=Path)
    p.add_argument('--hours',type=float,default=0);p.add_argument('--job-seconds',type=int,default=60)
    p.add_argument('--worker',type=Path);p.add_argument('--lock-fd',type=int)
    args=p.parse_args()
    if args.worker:worker(args.worker,args.lock_fd);return
    if args.activate:activate_prepared(args.activate);return
    if args.prepare and (args.run or args.queue or not args.launchd or not args.design_search):p.error('--prepare requires --launchd --design-search without --run/--queue')
    if args.launchd and sys.platform!='darwin':p.error('--launchd requires macOS')
    if args.launchd and args.queue and not args.run:
        p.error('--launchd requires --run admission; use a finalized predecessor, not an unlocked queue')
    if args.workers not in range(1,7)or args.hours<0 or args.job_seconds<1:p.error('workers1..6, nonnegative hours, positive job seconds')
    if args.triangle_target:
        if args.workers!=3 or not args.design_search or not args.native_components or args.mechanism_search or args.novelty_input:
            p.error('--triangle-target requires exactly 3 workers, --design-search, --native-components, and no other lane bundle')
    proposal=dict(predecessor=str(args.predecessor),workers=args.workers,gpu=False,
                  design_search=args.design_search,
                  native_components=args.native_components,
                  mechanism_search=args.mechanism_search,
                  novelty_search=args.novelty_input is not None,
                  supervisor='launchd' if args.launchd else 'detached_subprocess',
                  scope='identity-free positive-genus biased forest/joint11..16, solved15,18,20,21,24,25 and mixed-quotient12,15,16,18,20,21,24,25; saved unresolved replay first'
                        if args.design_search else 'identity-free positive-genus degrees11..14, biased forest/block lanes; saved unresolved replay first',
                  stopping='STOP sentinel, disk guard, or optional hours; no census claim',
                  activation='explicit finalized STOP handoff and shared lock release' if args.inherit_stopped else
                             'wait for predecessor complete independent audit and shared lock release')
    if args.triangle_target:
        if args.inherit_stopped is False:p.error('--triangle-target requires --inherit-stopped to retain completed prior evidence and safely hand off the lock')
        proposal['triangle_target']=dict(
            minimum_degree=16,maximum_degree=63,minimum_genus=3,
            closure_seconds=120,closure_cap=200000,
            passport_minimum_shards=[
                dict(passport_min=7,passport_max=7,pair_index_min=7,pair_index_max=7,
                     minimum_degree=16,maximum_degree=16,strong_divisor=8,exact_genus=3,
                     exact_passport=[7,7,7,15]),
                dict(passport_min=6,passport_max=12,minimum_degree=17,maximum_degree=63,
                     pair_index_min=6,pair_index_max=63,strong_divisor=8,exact_genus=3),
                dict(passport_min=13,passport_max=31,minimum_degree=27,maximum_degree=63,
                     pair_index_min=13,pair_index_max=63,strong_divisor=8,exact_genus=3)],
            previous_coverage='Reviewed genus-3 TRIANGLE work excludes degrees through 15. The user-reported degree-16 table leaves (7,7,7,15) as the search target; (7,7,8,14) is also ruled out by the reviewed high-index bound, while the other reported degree-16 computational exclusions remain certificate-audit items. The preserved six-worker constructor run made 660 jobs in degrees 11-18 and recorded zero full closures.',
            mathematical_filters='Exact source genus 3; TRIANGLE r2+r3<=d-1; reviewed h>=floor(d/2)+2g+1 bound; maximum-index entry fixed-point-free (ties allowed); lane0 exact passport (7,7,7,15); reviewed C12 exclusion of maximum-index type 15*2^((d-15)/2); reviewed equal-low theorems exclude any low involution and two special 3*2^5*1^3 low entries. Pair-index lower bounds match each disjoint r1 shard.',
            search_scope='Biased discovery. Each positive still requires full six-move closure and independent native audit.')
        proposal['scope']='Dedicated three-worker genus-3 TRIANGLE search: worker0 builds additive degree-16 low pairs from varied non-involution cycle types by an explicit two-component incidence forest, then tries up to64 strong/joint/random third entries per pair and logs each local rejection stage; workers1/2 cover exact-genus-3 passports from degree17 onward in disjoint r1=6..12 and r1=13..31 bands. Reviewed index/type exclusions run before closure. All positives require full six-move closure and independent native audit.'
    if args.mechanism_search:proposal['scope']='two primitive/prime-degree genus-one lanes, controlled mutations, forest baseline, quotient/ramified-double lift lane and one structural audit slot; four slots recover saved leads, two explore immediately'
    if not(args.run or args.queue or args.prepare):print(json.dumps(proposal,indent=2));return
    if args.checks is None:p.error('completed --checks RESULTS.json required')
    checks=json.loads(args.checks.read_text())
    compiled=compiled_source_files(args.compiled_source,checks.get('compiled_source_hashes')) if args.compiled_source else None
    for rel,digest in checks.get('launch_script_hashes',{}).items():
        path=PROJECT/rel
        if not rel.startswith('scripts/') or sha(path)!=digest:raise ValueError('launch script changed after scoped validation: '+rel)
    if checks['binary_sha256']!=sha(args.binary)or checks['controls']!=66 or checks['control_states']!=25192:raise ValueError('checks do not bind validated binary/catalog')
    if args.audit_binary and (not checks.get('passed') or checks.get('audit_binary_sha256')!=sha(args.audit_binary)):
        raise ValueError('checks do not bind validated native auditor')
    known_bank=args.checks.parent/'known_states.tsv'
    if args.design_search:
        if not args.audit_binary:raise ValueError('design search requires independent native auditor')
        if not checks.get('design_passed') or checks.get('known_bank_sha256')!=sha(known_bank):raise ValueError('design/bank checks required')
        proposal['scope']+='; exact known states skipped with immutable shared verified degree banks'
    if args.native_components and (not args.design_search or not checks.get('native_components_passed')):
        raise ValueError('native component audit requires its completed bound validation')
    profile_uv=profile_runtime(args.predecessor) if args.mechanism_search else None
    if args.mechanism_search:
        if not args.design_search or args.workers!=6 or not checks.get('mechanism_passed'):raise ValueError('mechanism search needs validated design and six total slots')
        if args.mutation_seeds is None or args.profile_input is None or not profile_uv:raise ValueError('bound mutation/profile inputs and isolated SymPy runtime required')
        if sha(args.mutation_seeds)!=checks['mutation_seeds_sha256'] or sha(args.profile_input)!=checks['profile_queue_sha256']:raise ValueError('mechanism input validation binding')
        if sha(args.checks.parent/'primitive_exclusions.json')!=checks['primitive_exclusions_sha256']:raise ValueError('retired family input binding')
        if sha(args.checks.parent/'double_bases.tsv')!=checks['double_bases_sha256']:raise ValueError('double-lift input binding')
    if args.novelty_input:
        if not args.mechanism_search or not checks.get('novelty_passed') or sha(args.novelty_input)!=checks['novelty_input_sha256']:raise ValueError('validated novelty input required')
        proposal['scope']='triangle passports11--18, multi-active-block return forests12--25, exact affine actions16, saved lead recovery and batched complete block-lattice profiles; six CPU slots and occasional known-family controls'
    if args.inherit_stopped:predecessor_finished(args.predecessor,True)
    control=PROJECT/'.run-control';control.mkdir(exist_ok=True)
    queue_lock=(control/'.forest-fast-queue.lock').open('a')
    fcntl.flock(queue_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    queue_guard=control/'FOREST_FAST_QUEUED.json'
    if queue_guard.exists():
        previous=json.loads(queue_guard.read_text())
        try:os.kill(previous['pid'],0)
        except ProcessLookupError:pass
        else:raise RuntimeError('a queued/running optimized successor already exists')
    lock=(control/'.campaign.lock').open('a')
    if args.run or args.prepare:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    dest=Path('data')
    if not storage_available(dest):raise RuntimeError('T7 mount/free-space guard')
    root=dest/('smooth_fast_'+dt.datetime.now(ZoneInfo('America/Chicago')).strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir();source=root/'source';source.mkdir();hashes={}
    paths=[PROJECT/'Cargo.toml',PROJECT/'Cargo.lock',*sorted((PROJECT/'src').rglob('*.rs')),
           *[s for s in sorted((PROJECT/'scripts').rglob('*.py'))if not s.name.startswith('._')],
           PROJECT/'scripts/gpu/forest_screen.swift',PROJECT/'scripts/gpu/forest_screen.metal',
           PROJECT/'scripts/gpu/forest_middle_screen.swift',PROJECT/'scripts/gpu/forest_middle_screen.metal',
           PROJECT/'evidence/forest_acceleration_20261012/PLAN.md',
           PROJECT/'evidence/forest_acceleration_20261012/MATH.md',
           PROJECT/'evidence/forest_second_pass_20261012/PLAN.md',
           PROJECT/'evidence/forest_third_pass_20261012/PLAN.md',
           *sorted((PROJECT/'evidence/smooth_discovery_design_20261012').glob('*')),
           *sorted((PROJECT/'evidence/smooth_theory_frontier_20261013').glob('*')),
           *sorted((PROJECT/'evidence/smooth_mechanism_breaker_20261013').glob('*')),
           *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
           PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv']
    for path in paths:
        rel=path.relative_to(PROJECT);target=source/rel;target.parent.mkdir(parents=True,exist_ok=True)
        if compiled is not None and str(rel).startswith('src/'):
            if str(rel) not in compiled:raise ValueError('workspace/build Rust file set differs')
            path=compiled[str(rel)]
        shutil.copy2(path,target);hashes[str(rel)]=sha(target)
    if compiled is not None:
        for rel,digest in checks['compiled_source_hashes'].items():
            if hashes.get(rel)!=digest:raise ValueError('frozen compiled source differs')
    shutil.copy2(args.binary,root/'smooth_forest_fast');shutil.copy2(args.checks,root/'VALIDATION.json')
    shutil.copy2(args.checks.parent/'controls.tsv',root/'CONTROLS.tsv')
    if args.design_search:shutil.copy2(known_bank,root/'KNOWN_STATES.tsv')
    if args.mechanism_search:
        shutil.copy2(args.mutation_seeds,root/'MUTATION_SEEDS.tsv');shutil.copy2(args.profile_input,root/'PROFILE_INPUT.json')
        shutil.copy2(args.checks.parent/'primitive_exclusions.json',root/'PRIMITIVE_EXCLUSIONS.json')
        shutil.copy2(args.checks.parent/'double_bases.tsv',root/'DOUBLE_BASES.tsv')
    if args.extra_leads:
        if not args.design_search:raise ValueError('extra leads require design search')
        shutil.copy2(args.extra_leads,root/'EXTRA_LEADS.json')
    if args.novelty_input:shutil.copy2(args.novelty_input,root/'NOVELTY_INPUT.json')
    if args.audit_binary:shutil.copy2(args.audit_binary,root/'smooth_negative_audit')
    save(root/'PLAN.json',dict(**proposal,started=now(),source_hashes=hashes,minimum_free_bytes=MIN_FREE_BYTES,
         binary_sha256=sha(root/'smooth_forest_fast'),validation_sha256=sha(root/'VALIDATION.json'),
         controls_sha256=sha(root/'CONTROLS.tsv'),
         known_bank_sha256=sha(root/'KNOWN_STATES.tsv') if args.design_search else None,
         extra_leads_sha256=sha(root/'EXTRA_LEADS.json') if args.extra_leads else None,
         novelty_input_sha256=sha(root/'NOVELTY_INPUT.json') if args.novelty_input else None,
         mutation_seeds_sha256=sha(root/'MUTATION_SEEDS.tsv') if args.mechanism_search else None,
         profile_input_sha256=sha(root/'PROFILE_INPUT.json') if args.mechanism_search else None,
         primitive_exclusions_sha256=sha(root/'PRIMITIVE_EXCLUSIONS.json') if args.mechanism_search else None,
         double_bases_sha256=sha(root/'DOUBLE_BASES.tsv') if args.mechanism_search else None,
         profile_uv=profile_uv,
         audit_binary_sha256=sha(root/'smooth_negative_audit') if args.audit_binary else None,
         inherit_stopped=args.inherit_stopped,
         deep_level=args.deep_level,
         block_deep_level=args.block_deep_level,
         campaign_lock=str(control/'.campaign.lock'),
         hours=args.hours,job_seconds=args.job_seconds,replay_seconds=120))
    command=[sys.executable,str(source/'scripts/smooth/run_forest_fast.py'),'--worker',str(root)]
    if args.prepare:
        state(root,status='prepared',workers=0,gpu=False,activation='explicit scheduled task; not running')
        lock.close();queue_lock.close();print('PREPARED='+str(root));return
    label='local.hurwitz.'+root.name if args.launchd else None
    if args.launchd:
        pid=spawn_launchd(root,command,label)
    else:
        with (root/'WORKER.log').open('x')as log:
            child=subprocess.Popen([*command,'--lock-fd',str(lock.fileno())],
                                   cwd=source,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,pass_fds=(lock.fileno(),))
        pid=child.pid
    replace_queue_guard(queue_guard,dict(pid=pid,output=str(root),queued=now(),workers=args.workers,launchd_label=label))
    save(root/'LAUNCH.json',dict(pid=pid,queued=now(),supervisor=proposal['supervisor'],launchd_label=label));lock.close()
    print('OUTPUT='+str(root));print('PID='+str(pid))
    queue_lock.close()

if __name__=='__main__':main()
