#!/usr/bin/env python3
"""4+2 smooth frontier: exact positive-genus degree10, guided degree11--14."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time
from zoneinfo import ZoneInfo

from run_smooth_degree9 import audit_parallel, rows, save, sha, state, now
from run_smooth_mission import audit_component, values
from odd_order_spin import parse_seed, hurwitz, cycle_lengths, inverse
from analyze_fixed_braid_constraints import product
from run_smooth_small_census import partition, pair_key_and_centralizer
from math import factorial

PROJECT = Path(__file__).resolve().parents[2]
PREVIOUS = Path('data/smooth_degree9_20261012_114225_888930')
# Each worker runs one single-threaded producer, then audits its saved leads.
# Degree10 is exclusively exact; block shapes require >=3 blocks and size>=2.
GUIDED = [[(11, 'forest'), (12, 'forest'), (13, 'forest'), (14, 'forest')],
          [(12, 'blocks'), (14, 'blocks')]]


def split_geodesics(p):
    """All shortest transposition paths to identity only split cycles."""
    seen = {p}; queue = [p]
    for q in queue:
        yield q
        remaining = set(range(len(q)))
        while remaining:
            root = min(remaining); cycle = []; x = root
            while x in remaining:
                remaining.remove(x); cycle.append(x); x = q[x]
            for i, x in enumerate(cycle):
                for y in cycle[i+1:]:
                    images = list(q); images[x], images[y] = images[y], images[x]
                    a = tuple(images)
                    if a not in seen: seen.add(a); queue.append(a)


def audit_split_frames(folder, n, _perms):
    frames = rows(folder/'FRAMES.tsv'); keys = set(); result = []
    for i, row in enumerate(frames):
        a, b, _, _ = parse_seed(f"[{row['a']}, {row['b']}, (), ()]", n)
        if int(row['frame']) != i or 2*n-len(cycle_lengths(a))-len(cycle_lengths(b)) != n-len(cycle_lengths(product(a,b))):
            raise ValueError('invalid geodesic frame')
        key, z = pair_key_and_centralizer(a,b)
        if key in keys or factorial(n)//z != int(row['pair_orbit']): raise ValueError('duplicate/wrong-mass frame')
        keys.add(key); result.append((a,b,factorial(n)//z))
    expected = set(); mass = 0
    from collections import Counter
    for lens in partition(n):
        p = list(range(n)); offset = 0
        for length in lens:
            for j in range(length): p[offset+j] = offset+(j+1)%length
            offset += length
        p = tuple(p); ip = n-len(lens); raw = 0
        for a in split_geodesics(p):
            b = product(inverse(a),p)
            if 2*n-len(cycle_lengths(a))-len(cycle_lengths(b)) != ip: raise ValueError('false split geodesic')
            key, _ = pair_key_and_centralizer(a,b); expected.add(key); raw += 1
        z = 1
        for length, multiplicity in Counter(lens).items(): z *= length**multiplicity*factorial(multiplicity)
        mass += raw*factorial(n)//z
    if keys != expected or mass != sum(w for _,_,w in result): raise ValueError('incomplete split population')
    return result


def audit_guided(folder):
    report = dict(sampled=True, construction=rows(folder/'COUNTS.tsv')[0], components=[])
    resolutions = rows(folder/'RESOLUTIONS.tsv')
    leads = []
    if (folder/'LEADS.tsv').exists():
        leads = [line.split('\t', 1) for line in (folder/'LEADS.tsv').read_text().splitlines()]
    if [(r['candidate'], r['seed']) for r in resolutions] != [(i, t) for i, t in leads]:
        raise ValueError('missing, reordered or unclassified durable lead')
    n = int(report['construction']['degree'])
    for row in resolutions:
        t = parse_seed(row['seed'], n)
        v = values(t)
        if not v['smooth'] or v['genus'] <= 0 or any(p == tuple(range(n)) for p in t):
            raise ValueError('out-of-scope guided lead')
        if row['status'] == 'nonsmooth':
            word = tuple(map(int, row['word'].split(',')))
            if not word or any(m not in range(6) for m in word):
                raise ValueError('invalid literal word')
            for m in word:
                t = hurwitz(t, m)
            if values(t)['smooth']:
                raise ValueError('false rejection')
        else:
            result = audit_component(folder/row['folder'])
            result.pop('state_keys')
            if result['summary']['status'] != row['status']:
                raise ValueError('resolution/certificate mismatch')
            report['components'].append(result)
    return report


def validation(exact, guided, root):
    # Independently audit both full and focused populations and compare the
    # focused closed components against the positive-genus part of the full run.
    reports = {}
    for focus in [False, True]:
        folder = root/('focused' if focus else 'full')
        cmd = [str(exact), '--d', '5', '--workers', '4', '--seconds', '120', '--out', str(folder)]
        if focus:
            cmd.append('--positive-only')
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        reports[str(focus)] = audit_parallel(folder, 4, audit_split_frames)
    full = reports['False']; focused = reports['True']
    full_positive = sorted(c['table_sha256'] for c in full['components'].values() if int(c['summary']['genus']) > 0)
    focus_positive = sorted(c['table_sha256'] for c in focused['components'].values())
    if full_positive != focus_positive or not all(r['coverage']['fully_resolved'] == 'true' for r in reports.values()):
        raise ValueError('focus altered positive-genus components')
    controls = []
    for n, lane in [(11, 'forest'), (12, 'blocks'), (14, 'blocks')]:
        folder = root/f'guided_{n}_{lane}'
        subprocess.run([str(guided), '--d', str(n), '--lane', lane, '--seed', '41', '--seconds', '2',
                        '--out', str(folder)], check=True, stdout=subprocess.DEVNULL)
        controls.append(audit_guided(folder))
    result = dict(full_and_focused_degree5_independently_audited=True,
                  focused_positive_components_equal_full=True, guided_smokes=controls)
    save(root/'VALIDATION.json', result)
    return result


def worker(root, lock_fd):
    lock = os.fdopen(lock_fd, 'a')
    plan = json.loads((root/'PLAN.json').read_text())
    children = set(); mutex = threading.Lock(); stop = threading.Event()
    helper = None
    def run(cmd, logfile):
        with logfile.open('x') as log:
            child = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
            with mutex: children.add(child)
            try:
                while child.poll() is None:
                    if shutil.disk_usage(root).free < 30*2**30:
                        child.terminate(); child.wait()
                        raise RuntimeError('T7 free-space guard; partial evidence preserved')
                    time.sleep(1)
                if child.returncode:
                    raise RuntimeError(f'child failed: {cmd[0]} exit {child.returncode}')
            finally:
                with mutex: children.discard(child)
    def guided(lane_id):
        epoch = 0; completed = []
        while not stop.is_set():
            n, lane = GUIDED[lane_id][epoch % len(GUIDED[lane_id])]
            name = f'guided_{lane_id}_{epoch:06}_{n}_{lane}'
            folder = root/name
            run([str(root/'smooth_guided'), '--d', str(n), '--lane', lane, '--seed', str(1000003*lane_id+epoch+81),
                 '--seconds', str(plan['guided_job_seconds']), '--out', str(folder)], root/f'{name}.log')
            checks = audit_guided(folder)
            save(folder/'CHECKS.json', checks)
            (folder/'VERIFIED').write_text(now()+'\n')
            completed.append(dict(folder=name, construction=checks['construction']))
            state(root/f'lane_{lane_id}', status='guided', jobs=len(completed), latest=name, worker=1)
            print('GUIDED_VERIFIED', name, json.dumps(checks['construction']), flush=True)
            epoch += 1
        return completed
    try:
        for name, digest in plan['source_hashes'].items():
            if sha(root/'source'/name) != digest: raise ValueError('frozen source changed')
        for name, digest in plan['binary_hashes'].items():
            if sha(root/name) != digest: raise ValueError('frozen binary changed')
        if shutil.which('caffeinate'):
            helper = subprocess.Popen(['caffeinate', '-i', '-w', str(os.getpid())])
        for i in range(2): (root/f'lane_{i}').mkdir()
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(guided, i) for i in range(2)]
            try:
                state(root, status='enumerating', pid=os.getpid(), exact_workers=4, guided_workers=2)
                folder = root/'degree_10'
                run([str(root/'smooth_parallel_census'), '--d', '10', '--workers', '4', '--positive-only',
                     '--seconds', '31536000', '--out', str(folder), '--cap', '100000',
                     '--closure-seconds', '30'], root/'degree_10.log')
                print('ENUMERATED', json.dumps(rows(folder/'COVERAGE.tsv')[0]), flush=True)
                state(root, status='auditing', pid=os.getpid(), exact_workers=4, guided_workers=2)
                checks = audit_parallel(folder, 4, audit_split_frames)
                save(root/'EXACT_CHECKS.json', checks)
            finally:
                stop.set()
            sampled = [f.result() for f in futures]
        complete = checks['exhaustive'] and checks['coverage']['fully_resolved'] == 'true'
        save(root/'FINISHED.json', dict(ended=now(), status='complete' if complete else 'incomplete',
                                     guided_batches=sampled, scope='positive-genus identity-free degree10'))
        state(root, status='complete' if complete else 'incomplete', workers=0)
    except BaseException as error:
        stop.set()
        with mutex:
            for child in children:
                if child.poll() is None: child.terminate()
        save(root/'FAILURE.json', dict(error=repr(error), ended=now()))
        state(root, status='failed', workers=0, error=repr(error))
        raise
    finally:
        if helper is not None and helper.poll() is None: helper.terminate()
        lock.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true')
    p.add_argument('--worker', type=Path); p.add_argument('--lock-fd', type=int)
    p.add_argument('--guided-job-seconds', type=int, default=60)
    args = p.parse_args()
    if args.worker: worker(args.worker, args.lock_fd); return
    if args.guided_job_seconds <= 0: p.error('positive guided-job-seconds required')
    if not args.run:
        print(json.dumps(dict(exact='degree10, identity-free, positive-genus, all passports/transitive groups',
                             workers='4 exact +2 guided', guided=GUIDED, gpu=False,
                             stop='guided stops after exact independent audit, preserving caps'), indent=2)); return
    if json.loads((PREVIOUS/'FINISHED.json').read_text())['status'] != 'complete':
        raise RuntimeError('degree9 predecessor has not passed its independent audit')
    dest = Path('data')
    if not Path('data').is_dir() or shutil.disk_usage(dest).free < 30*2**30:
        raise RuntimeError('mounted T7 with 30GiB free required')
    control = PROJECT/'.run-control'; control.mkdir(exist_ok=True)
    lock = (control/'.campaign.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = dest/('smooth_frontier_'+dt.datetime.now(ZoneInfo('America/Chicago')).strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir(); source = root/'source'; source.mkdir(); hashes = {}
    paths = [PROJECT/'Cargo.toml', PROJECT/'Cargo.lock', *sorted((PROJECT/'src').rglob('*.rs')),
             *[p for p in sorted((PROJECT/'scripts').rglob('*.py')) if not p.name.startswith('._')],
             *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
             PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv',
             PROJECT/'evidence/smooth_frontier_20261012/PLAN.md']
    for path in paths:
        rel = path.relative_to(PROJECT); target = source/rel; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target); hashes[str(rel)] = sha(target)
    env = dict(os.environ, CARGO_BUILD_JOBS='1', RUSTFLAGS='-C target-cpu=native',
               CARGO_TARGET_DIR=str(PROJECT/'target/smooth-parallel'))
    print('OUTPUT='+str(root), flush=True)
    with (root/'BUILD.log').open('x') as log:
        for binary in ['smooth_parallel_census', 'smooth_guided']:
            subprocess.run(['cargo', 'test', '--release', '--bin', binary, '_tests', '--', '--test-threads=1'],
                           cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cargo', 'build', '--release', '--bin', 'smooth_parallel_census', '--bin', 'smooth_guided'],
                       cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, str(source/'scripts/smooth/test_smooth_degree9.py')],
                       cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, str(source/'scripts/smooth/test_smooth_frontier.py')],
                       cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
    binaries = {}
    for name in ['smooth_parallel_census', 'smooth_guided']:
        shutil.copy2(PROJECT/'target/smooth-parallel/release'/name, root/name)
        binaries[name] = sha(root/name)
    v = root/'validation'; v.mkdir()
    validation(root/'smooth_parallel_census', root/'smooth_guided', v)
    save(root/'PLAN.json', dict(started=now(), predecessor=str(PREVIOUS), source_hashes=hashes,
             binary_hashes=binaries, exact_workers=4, guided_workers=2, gpu=False,
             guided_job_seconds=args.guided_job_seconds, guided=GUIDED,
             exact_scope='all degree10 positive-genus identity-free everywhere-smooth components modulo S10'))
    with (root/'WORKER.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(source/'scripts/smooth/run_smooth_frontier.py'), '--worker', str(root),
                                  '--lock-fd', str(lock.fileno())], cwd=source, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True, pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json', dict(pid=child.pid, started=now()))
    save(control/f'SMOOTH_FRONTIER_{root.name}.json', dict(output=str(root), pid=child.pid, workers=6))
    lock.close(); print('PID='+str(child.pid), flush=True)


if __name__ == '__main__': main()
