"""Finalize a dead smooth controller, retaining all evidence and pending leads.

Only independently verified literal negatives from interrupted journals may be
retired. Adapter COUNTS are verifier bounds, NOT original proposal statistics.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import fcntl
import json
import os
import shlex
from pathlib import Path
import subprocess

from audit_guided_native import audit_guided_native
from run_smooth_degree9 import sha, save, state, now

PROJECT=Path(__file__).resolve().parents[2]
HEADER='candidate\tstatus\tword\tfolder\tseed\n'
COUNT_HEADER='degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'


def negative_prefix(folder,dest,n,binary):
    """Copy matching complete negative rows; original journals remain untouched."""
    dest.mkdir()
    count=0;maximum=0;last=0;partial=False
    with (folder/'LEADS.tsv').open() as leads,(folder/'RESOLUTIONS.tsv').open() as resolutions,\
         (dest/'LEADS.tsv').open('x') as out_leads,(dest/'RESOLUTIONS.tsv').open('x') as out_res:
        header=resolutions.readline()
        if header!=HEADER:
            if not header.endswith('\n') and HEADER.startswith(header):partial=True
            else:raise ValueError('recovery resolution header')
        out_res.write(HEADER)
        for line in resolutions:
            if not line.endswith('\n'):partial=True;break
            fields=line.rstrip('\n').split('\t')
            lead=leads.readline()
            if not lead.endswith('\n') or len(fields)!=5 or lead!=fields[0]+'\t'+fields[4]+'\n':
                raise ValueError('recovery durable lead mismatch')
            candidate=int(fields[0])
            if candidate<=last:raise ValueError('recovery candidate ordering')
            last=candidate
            if fields[1]=='nonsmooth':
                out_leads.write(lead);out_res.write(line);count+=1;maximum=candidate
            elif fields[1] not in ['closed_smooth','capped','timeout']:
                raise ValueError('unsupported recovery status')
    # Match the native verifier's original format. The only needed quantitative
    # fields are candidate upper bound and zero exceptional outcomes. Do not
    # report these synthetic rows as completed construction work.
    (dest/'COUNTS.tsv').write_text(COUNT_HEADER+f'{n}\treplay\t0\t{maximum}\t{count}\t{count}\t{count}\t0\t0\t0\t0\t0\n')
    checks=audit_guided_native(dest,binary)
    if checks['independent_native']['nonsmooth']!=count:raise ValueError('recovery negative mass')
    save(dest/'CHECKS.json',checks)
    return dict(prefix=str(dest),records=count,trailing_incomplete_resolution=partial,
                original_sources={name:sha(folder/name) for name in ['LEADS.tsv','RESOLUTIONS.tsv']},
                checks_sha256=sha(dest/'CHECKS.json'),
                quantitative_scope='literal-negative journal recovery only; adapter COUNTS are not original proposal totals')


def excluded_ids(root,folder):
    """Exact, hash-bound skip set, empty if this job has no audited recovery."""
    path=root/'RECOVERY.json'
    if not path.exists():return set()
    record=json.loads(path.read_text()).get('negative_prefixes',{}).get(folder.name)
    if record is None:return set()
    for name,digest in record['original_sources'].items():
        if sha(folder/name)!=digest:raise ValueError('interrupted journal changed after recovery')
    prefix=Path(record['prefix']);checks_path=prefix/'CHECKS.json'
    if sha(checks_path)!=record['checks_sha256']:raise ValueError('recovery checks changed')
    checks=json.loads(checks_path.read_text())
    for name,digest in checks['sources'].items():
        if sha(prefix/name)!=digest:raise ValueError('audited recovery prefix changed')
    native=checks['independent_native']
    if native['records']!=record['records'] or native['nonsmooth']!=record['records'] or native['closed'] or native['unresolved']:
        raise ValueError('recovery is not a complete negative-only proof')
    with (prefix/'LEADS.tsv').open() as stream:
        ids={line.split('\t',1)[0] for line in stream}
    if len(ids)!=record['records']:raise ValueError('recovery IDs not unique')
    return ids


def ensure_inactive(root,pid):
    try:os.kill(pid,0)
    except ProcessLookupError:pass
    else:raise RuntimeError('controller PID still exists; do not finalize or kill it')
    rows=subprocess.check_output(['ps','-axo','pid,command'],text=True)
    for line in rows.splitlines()[1:]:
        parts=line.strip().split(None,1)
        if len(parts)!=2:continue
        try:argv=shlex.split(parts[1])
        except ValueError:continue
        if not argv:continue
        executable=Path(argv[0]).name
        if executable.startswith('python') and str(root/'source/scripts/smooth/run_forest_fast.py') in argv:
            raise RuntimeError('bound controller still alive')
        if argv[0]==str(root/'smooth_forest_fast') and '--out' in argv:
            target=Path(argv[argv.index('--out')+1]).resolve()
            if target.is_relative_to(root):raise RuntimeError('bound producer still alive')
        if argv[0]==str(root/'smooth_negative_audit') and '--guided' in argv:
            target=Path(argv[argv.index('--guided')+1]).resolve()
            if target.is_relative_to(root):raise RuntimeError('bound auditor still alive')


def recovered_epochs(root, plan, old_state):
    workers=int(plan['workers'])
    epochs=list(old_state.get('epochs',[0]*workers))
    if len(epochs)!=workers:raise ValueError('recovery worker/epoch shape')
    for folder in root.glob('guided_*'):
        parts=folder.name.split('_')
        if folder.is_dir() and len(parts)==5:
            worker,epoch=int(parts[1]),int(parts[2])
            if not 0<=worker<workers:raise ValueError('recovery worker outside plan')
            epochs[worker]=max(epochs[worker],epoch+1)
    profile=root/'PROFILE_SUMMARY.json'
    if plan.get('mechanism_search') and profile.exists():
        epochs[-1]=max(epochs[-1],int(json.loads(profile.read_text())['checked']))
    return epochs


def preserve_original_state(root):
    # A repeated recovery must never replace the original failed snapshot.
    with (root/'STATE_BEFORE_RECOVERY.json').open('xb') as stream:
        stream.write((root/'STATE.json').read_bytes())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--finalize',action='store_true')
    args=p.parse_args();root=args.root.resolve()
    plan=json.loads((root/'PLAN.json').read_text());launch=json.loads((root/'LAUNCH.json').read_text())
    lock=(PROJECT/'.run-control/.campaign.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    ensure_inactive(root,int(launch['pid']))
    if (root/'FINISHED.json').exists():raise ValueError('already finalized; do not overwrite')
    if not args.finalize:
        print(json.dumps(dict(root=str(root),inactive=True,action='audit complete jobs and negative journal prefixes; preserve all raw files and recover all remaining leads')));return
    if Path(json.loads((PROJECT/'.run-control/FOREST_FAST_QUEUED.json').read_text())['output']).resolve()!=root:
        raise ValueError('root is not the current stopped-by-interruption pointer')
    if sha(root/'smooth_negative_audit')!=plan['audit_binary_sha256']:raise ValueError('frozen auditor changed')
    old_state=json.loads((root/'STATE.json').read_text());epochs=recovered_epochs(root,plan,old_state)
    original=root/'STATE_BEFORE_RECOVERY.json'
    if original.exists():
        if original.read_bytes()!=(root/'STATE.json').read_bytes():raise ValueError('recovery retry no longer matches original failed state')
    else:preserve_original_state(root)
    completed=[];prefixes={};partial=[]
    for folder in sorted(root.glob('guided_*')):
        if not folder.is_dir():continue
        _,worker,epoch,n,lane=folder.name.split('_');worker,epoch,n=map(int,(worker,epoch,n))
        epochs[worker]=max(epochs[worker],epoch+1)
        if (folder/'VERIFIED').exists():continue
        if (folder/'COUNTS.tsv').exists() and (folder/'COUNTS.tsv').stat().st_size:
            checks=audit_guided_native(folder,root/'smooth_negative_audit')
            save(folder/'CHECKS.json',checks);(folder/'VERIFIED').write_text(now()+'\n');completed.append(folder.name)
        else:
            partial.append(folder.name)
            if (folder/'LEADS.tsv').exists() and (folder/'RESOLUTIONS.tsv').exists():
                dest=root/('recovered_negative_'+folder.name)
                if dest.exists():
                    import time
                    # Preserve incomplete products of a failed recovery; never raw journals.
                    dest.rename(dest.with_name(dest.name+'_prior_'+str(time.time_ns())))
                prefixes[folder.name]=negative_prefix(folder,dest,n,root/'smooth_negative_audit')
    record=dict(ended=now(),cause=old_state.get('error','controller and bound children absent; termination cause unknown'),
                original_state=old_state,epochs=epochs,newly_verified_jobs=completed,
                interrupted_jobs=partial,negative_prefixes=prefixes,
                invariant='raw journals/frozen source unchanged; all other durable starts remain unresolved')
    save(root/'RECOVERY.json',record)
    (root/'STOP').write_text('Finalized inactive controller; raw partial leads preserved.\n')
    save(root/'FINISHED.json',dict(status='stopped',ended=now(),epochs=epochs,interrupted=True,
         recovery_sha256=sha(root/'RECOVERY.json')))
    state(root,status='stopped',workers=0,reason='dead controller recovered; pending raw starts inherited')
    print(json.dumps(dict(root=str(root),epochs=epochs,newly_verified_jobs=len(completed),
                         checked_negative_prefix_records=sum(r['records'] for r in prefixes.values()),interrupted_jobs=len(partial))))


if __name__=='__main__':main()
