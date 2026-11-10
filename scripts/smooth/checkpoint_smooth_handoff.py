#!/usr/bin/env python3
"""Explicit validated handoff; no deletion or modification of frozen inputs."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse, json, os, signal, subprocess, time
from pathlib import Path
from run_negative_audit import save, sha, state, now

def command(pid):
    return subprocess.check_output(['ps','-p',str(pid),'-o','command='],text=True).strip()

def main():
    p=argparse.ArgumentParser();p.add_argument('--stop-old',action='store_true');p.add_argument('--predecessor',type=Path,required=True);p.add_argument('--queued',type=Path,required=True)
    p.add_argument('--finalize-stopped',action='store_true')
    args=p.parse_args();old=args.predecessor;queued=args.queued
    old_launch=json.loads((old/'LAUNCH.json').read_text());old_pid=old_launch['pid'];q_pid=json.loads((queued/'LAUNCH.json').read_text())['pid']
    if args.finalize_stopped:
        if json.loads((queued/'STATE.json').read_text())['status']!='stopped' or not (queued/'CANCELLED.json').exists():raise ValueError('no cancelled queue checkpoint')
        if subprocess.run(['pgrep','-g',str(old_pid)],capture_output=True).returncode!=1:raise ValueError('old group still active')
    else:
        if json.loads((queued/'STATE.json').read_text())['status']!='queued':raise ValueError('successor not idle queued')
        if str(queued/'source/scripts/smooth/run_forest_fast.py') not in command(q_pid):raise ValueError('queued PID identity mismatch')
        if str(old/'source/scripts/smooth/run_smooth_frontier.py') not in command(old_pid):raise ValueError('controller PID identity mismatch')
        if os.getpgid(old_pid)!=old_pid or os.getpgid(old_pid)==os.getpgrp():raise ValueError('unsafe process group')
    if not args.stop_old:print(json.dumps(dict(predecessor=str(old),pid=old_pid,queued=str(queued),queued_pid=q_pid,action='cancel idle queue; stop validated old process group; preserve all files')));return
    # Cancel the idle waiter FIRST, so it cannot accidentally acquire the
    # predecessor lock before the new independent audit is launched.
    if not args.finalize_stopped:
        os.kill(q_pid,signal.SIGTERM)
        for _ in range(100):
            try:os.kill(q_pid,0)
            except ProcessLookupError:break
            time.sleep(.05)
        else:raise RuntimeError('queued waiter did not stop')
        save(queued/'CANCELLED.json',dict(ended=now(),reason='superseded by degree10 native-audit handoff',pid=q_pid))
        state(queued,status='stopped',workers=0,reason='superseded before activation')
    process_rows=subprocess.check_output(['ps','-axo','pid,ppid,pgid,command'],text=True)
    # Freeze all producer/auditor children before capturing the recovery list;
    # SIGTERM+SIGCONT then retires precisely this validated controller group.
    if not args.finalize_stopped:os.killpg(old_pid,signal.SIGSTOP)
    try:
        completed={};invalid=[]
        for file in sorted((old/'degree_10/audit_frames').glob('*.json')):
            if file.name.startswith('._'):continue
            try:json.loads(file.read_text())
            except json.JSONDecodeError:invalid.append(str(file))
            else:completed[file.name]=sha(file)
        partial=[str(f) for f in old.glob('guided_*') if f.is_dir() and not(f/'VERIFIED').exists()]
        save(old/'STOPPED_CHECKPOINT.json',dict(ended=now(),reason='user-authorized native audit acceleration',pid=old_pid,
            python_completed_frame_reports=completed,unfinished_python_reports=invalid,partial_guided_jobs=partial,
            producer_table_hashes={f.name:sha(f) for f in (old/'degree_10').glob('*.tsv')},
            processes=[line for line in process_rows.splitlines() if line.split() and len(line.split())>=3 and line.split()[2]==str(old_pid)],
            scope='original complete producer unchanged; Python audit incomplete, never marked complete'))
    finally:
        if not args.finalize_stopped:
            os.killpg(old_pid,signal.SIGTERM);os.killpg(old_pid,signal.SIGCONT)
    for _ in range(200):
        output=subprocess.run(['pgrep','-g',str(old_pid)],capture_output=True,text=True)
        if output.returncode==1:break
        time.sleep(.05)
    else:raise RuntimeError('old process group still present; no new pool authorized until gone')
    state(old,status='stopped',workers=0,reason='native audit successor; producer and raw leads preserved')
    print('CHECKPOINT='+str(old/'STOPPED_CHECKPOINT.json'));print('PYTHON_FRAMES='+str(len(completed)));print('PARTIAL_GUIDED='+str(len(partial)))

if __name__=='__main__':main()
