#!/usr/bin/env python3
"""Freeze an independent re-audit of completed small-census outputs; no rerun."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import sys

from run_smooth_mission import rows, save, sha

PROJECT = Path(__file__).resolve().parents[2]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(); root = args.root.resolve(); out = args.out.resolve()
    lock = (PROJECT/'.run-control/.campaign.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    plan = json.loads((root/'PLAN.json').read_text())
    for name, digest in plan['source_hashes'].items():
        if sha(root/'source'/name) != digest:
            raise ValueError('changed original frozen source')
    if sha(root/'smooth_small_census') != plan['binary_sha256']:
        raise ValueError('changed original binary')
    if not all(rows(root/f'degree_{n}/COVERAGE.tsv')[0]['enumeration_complete']=='true' for n in plan['degrees']):
        raise ValueError('incomplete producer population')
    out.mkdir(); source = out/'source'; source.mkdir(); hashes = {}
    for name in ['verify_smooth_small_run.py', 'run_smooth_small_census.py', 'run_smooth_mission.py', 'run_s18_linked_roots.py', 'analyze_fixed_braid_constraints.py', 'odd_order_spin.py', 'test_smooth_small_census.py']:
        shutil.copy2(_scripts_path.find(name), source/name); hashes[name] = sha(source/name)
    inputs = {str(path.relative_to(root)): sha(path) for path in root.glob('degree_*/*.tsv')}
    inputs.update({str(path.relative_to(root)): sha(path) for path in root.glob('degree_*/component_*/*.tsv')})
    save(out/'PLAN.json', dict(input_root=str(root), original_plan_sha256=sha(root/'PLAN.json'),
                              binary_sha256=plan['binary_sha256'], verifier_hashes=hashes, input_hashes=inputs,
                              started=dt.datetime.now().isoformat(), correction='Auditor frame parser now uses comma-space separators. Producer data unchanged.'))
    def verify(n):
        with (out/f'degree_{n}.log').open('x') as log:
            subprocess.run([sys.executable, str(source/'run_smooth_small_census.py'), '--audit', str(root/f'degree_{n}'), '--audit-out', str(out/f'degree_{n}_CHECKS.json')], cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
        answer = json.loads((out/f'degree_{n}_CHECKS.json').read_text()); print('AUDITED', n, flush=True); return answer
    print('AUDIT_OUTPUT='+str(out), flush=True)
    try:
        with ThreadPoolExecutor(max_workers=min(4, len(plan['degrees']))) as pool:
            jobs = list(pool.map(verify, plan['degrees']))
        if any(sha(root/name) != digest for name, digest in inputs.items()):
            raise ValueError('input changed during audit')
        complete = all(j['exhaustive'] and j['coverage']['fully_resolved']=='true' for j in jobs)
        save(out/'CHECKS.json', dict(jobs=jobs, fully_resolved=complete, input_root=str(root), input_plan_sha256=sha(root/'PLAN.json')))
        save(out/'FINISHED.json', dict(status='complete' if complete else 'incomplete', ended=dt.datetime.now().isoformat()))
    except Exception as e:
        save(out/'FAILURE.json', dict(error=repr(e))); raise
    finally:
        lock.close()


if __name__=='__main__':
    main()
