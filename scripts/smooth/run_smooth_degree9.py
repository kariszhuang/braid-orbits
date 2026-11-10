#!/usr/bin/env python3
"""Frozen degree-nine all-passport census; six Rust threads, then six audit processes."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import datetime as dt
import fcntl
import itertools
import json
from math import factorial
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time
from zoneinfo import ZoneInfo

from analyze_fixed_braid_constraints import product, point_orbits
from odd_order_spin import inverse, hurwitz, cycle_lengths
from run_smooth_mission import audit_component, rows, save, sha, values
from run_s18_linked_roots import state_key
from run_smooth_small_census import audit_frames, audit as serial_audit

PROJECT = Path(__file__).resolve().parents[2]
COLS = ['enumerated', 'transitive', 'both_geodesic', 'nonidentity', 'identity',
        'nonsmooth_orbit', 'closed_occurrences', 'unresolved']
_PERMS = None
_CLOSED = None


def now():
    return dt.datetime.now(ZoneInfo('America/Chicago')).isoformat()


def state(root, **data):
    # Mutable lifecycle pointer, unlike immutable proof/plan artifacts.
    tmp = root/'STATE.tmp'
    tmp.write_text(json.dumps(dict(updated=now(), **data), indent=2)+'\n')
    tmp.replace(root/'STATE.json')


def init_audit(n, closed):
    global _PERMS, _CLOSED
    _PERMS = [(p, n-len(cycle_lengths(p))) for p in itertools.permutations(range(n))]
    _CLOSED = closed


def binary_records(path):
    with path.open('rb') as stream:
        while block := stream.read(8*131072):
            if len(block) % 8:
                raise ValueError('truncated binary record')
            yield from struct.iter_unpack('<II', block)


def audit_frame(task):
    folder, fi, a, b, expected = task
    folder = Path(folder)/'frames'/f'frame_{fi:05}'
    n = len(a)
    if not (folder/'COMPLETE').exists():
        raise ValueError('missing complete frame marker')
    dictionary = rows(folder/'WITNESSES.tsv')
    if [int(row['token']) for row in dictionary] != list(range(len(dictionary))):
        raise ValueError('dictionary numbering')
    witnesses = []
    for row in dictionary:
        status, data = row['status'], row['word_or_component']
        if status == 'nonsmooth':
            data = tuple(map(int, data.split(',')))
            if not data or any(not 0 <= move < 6 for move in data):
                raise ValueError('bad witness word')
        elif status != 'closed_smooth':
            # Partial closure stays unresolved, never an exclusion.
            audit_component(folder.parents[1]/data)
        witnesses.append((status, data))
    records = iter(binary_records(folder/'OUTCOMES.bin'))
    p = product(a, b); ip = n-len(cycle_lengths(p))
    left = len(point_orbits((a, b))); identity = tuple(range(n))
    count = [0]*8
    positive_only = (folder.parents[1]/'POSITIVE_ONLY').exists()
    omitted = [0, 0]
    for rank, (c, ic) in enumerate(_PERMS):
        count[0] += 1
        d = inverse(product(p, c)); t = (a, b, c, d)
        if positive_only:
            if any(q == identity for q in t):
                omitted[0] += 1; continue
            if ip+ic+n-len(cycle_lengths(d)) <= 2*n-2:
                omitted[1] += 1; continue
        if len(point_orbits(t[:3])) != 1:
            continue
        count[1] += 1
        right = len(point_orbits(t[2:]))
        if n-ip-left-right+1 != 0 or right != 1:
            raise ValueError('independent implication violation')
        if ic+n-len(cycle_lengths(d)) == ip:
            count[2] += 1
        if any(q == identity for q in t):
            count[4] += 1; continue
        count[3] += 1
        record_rank, token = next(records)
        if record_rank != rank or token >= len(witnesses):
            raise ValueError('missing/reordered/out-of-range outcome')
        status, data = witnesses[token]
        if status == 'nonsmooth':
            for move in data:
                t = hurwitz(t, move)
            if values(t)['smooth']:
                raise ValueError('false literal negative witness')
            count[5] += 1
        elif status == 'closed_smooth':
            if state_key(t) not in _CLOSED[data]:
                raise ValueError('false component membership')
            count[6] += 1
        else:
            count[7] += 1
    if next(records, None) is not None or count != expected:
        raise ValueError('extra outcomes or per-frame coverage mismatch')
    local = rows(folder/'COUNTS.tsv')
    if len(local) != 1 or int(local[0]['frame']) != fi or [int(local[0][c]) for c in COLS] != count:
        raise ValueError('local/global count mismatch')
    if (folder/'COMPLETE').read_text() != f'{fi}\t{factorial(n)}\n':
        raise ValueError('wrong completion marker')
    if positive_only:
        row = rows(folder/'OMITTED.tsv')[0]
        if omitted != [int(row['identity']), int(row['nonpositive_genus'])]:
            raise ValueError('scope omission mismatch')
    return dict(frame=fi, count=count, sources={name: sha(folder/name) for name in
                ['OUTCOMES.bin', 'WITNESSES.tsv', 'COUNTS.tsv', 'COMPLETE'] +
                (['OMITTED.tsv'] if positive_only else [])})


def audit_parallel(folder, workers=6, frame_verifier=None):
    folder = Path(folder)
    coverage = rows(folder/'COVERAGE.tsv')[0]; n = int(coverage['degree'])
    if coverage['enumeration_complete'] != 'true':
        return dict(exhaustive=False, coverage=coverage, warning='Unfinished construction is not a census.')
    perms = [(p, n-len(cycle_lengths(p))) for p in itertools.permutations(range(n))]
    frames = (frame_verifier or audit_frames)(folder, n, perms)
    del perms
    counts = rows(folder/'FRAME_COUNTS.tsv')
    if len(counts) != len(frames) or [int(r['frame']) for r in counts] != list(range(len(frames))):
        raise ValueError('frame gaps/duplicates')
    components = {p.name: audit_component(p) for p in sorted(folder.glob('component_*'))}
    closed = {name: set(c['state_keys']) for name, c in components.items()}
    if sum(map(len, closed.values())) != len(set().union(*closed.values())):
        raise ValueError('duplicated closed components')
    audit_dir = folder/'audit_frames'; audit_dir.mkdir()
    tasks = [(str(folder), fi, a, b, [int(counts[fi][c]) for c in COLS])
             for fi, (a, b, _) in enumerate(frames)]
    total = [0]*8; weighted = [0, 0]
    with ProcessPoolExecutor(max_workers=workers, initializer=init_audit, initargs=(n, closed)) as pool:
        for result in pool.map(audit_frame, tasks, chunksize=1):
            fi = result['frame']; count = result['count']
            save(audit_dir/f'frame_{fi:05}.json', result)
            total = [x+y for x, y in zip(total, count)]
            weighted[0] += frames[fi][2]*count[1]; weighted[1] += frames[fi][2]*count[2]
            if (fi+1) % 10 == 0 or fi+1 == len(frames):
                print(f'AUDIT frames={fi+1}/{len(frames)}', flush=True)
    names = ['enumerated', 'transitive_smooth_occurrences', 'both_geodesic_occurrences',
             'nonidentity_occurrences', 'identity_occurrences', 'nonsmooth_orbit_occurrences',
             'closed_occurrences', 'unresolved']
    if total != [int(coverage[x]) for x in names] or len(components) != int(coverage['closed_components']):
        raise ValueError('global coverage mismatch')
    if [weighted[0], weighted[1], 2*weighted[0]-weighted[1]] != [int(coverage[x]) for x in
            ['labeled_left_smooth', 'labeled_both_smooth', 'labeled_smooth_tuples']]:
        raise ValueError('labeled mass mismatch')
    for c in components.values():
        c.pop('state_keys')
    return dict(exhaustive=True, coverage=coverage, components=components,
                independently_reenumerated_third_occurrences=total[0],
                literal_orbit_outcomes_checked=total[3], audit_workers=workers,
                sources={p.name: sha(p) for p in folder.glob('*.tsv')})


def verify_small(binary, serial, root):
    reports = {}
    for kind, executable, workers in [('serial', serial, 1), ('parallel1', binary, 1), ('parallel6', binary, 6)]:
        folder = root/kind
        cmd = [str(executable), '--d', '4', '--out', str(folder), '--seconds', '120']
        if kind != 'serial':
            cmd += ['--workers', str(workers)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        check = serial_audit(folder) if kind == 'serial' else audit_parallel(folder, workers)
        reports[kind] = check
    fields = [k for k in reports['serial']['coverage'] if k != 'elapsed_seconds']
    baseline = [reports['serial']['coverage'][k] for k in fields]
    for check in reports.values():
        if [check['coverage'][k] for k in fields] != baseline:
            raise ValueError('serial/parallel census discrepancy')
    # Matched complete n7 runs: actual medians, not an assumed 6x gain.
    timings = {1: [], 6: []}
    for repeat in range(3):
        for workers in [1, 6] if repeat % 2 == 0 else [6, 1]:
            folder = root/f'benchmark_{workers}_{repeat}'
            subprocess.run([str(binary), '--d', '7', '--workers', str(workers), '--out', str(folder),
                            '--seconds', '120'], check=True, stdout=subprocess.DEVNULL)
            cov = rows(folder/'COVERAGE.tsv')[0]
            if cov['fully_resolved'] != 'true':
                raise ValueError('benchmark incomplete')
            timings[workers].append(float(cov['elapsed_seconds']))
    import statistics
    result = dict(s4_serial_parallel_independently_audited=True,
                  timings_seconds=timings, degree7_median_speedup=
                  statistics.median(timings[1])/statistics.median(timings[6]),
                  limit='Small regression panel; not a degree9 or whole-pipeline speedup.')
    save(root/'VALIDATION.json', result)
    return result


def worker(root, lock_fd):
    lock = os.fdopen(lock_fd, 'a')
    plan = json.loads((root/'PLAN.json').read_text()); caffeinate = None
    child = None
    try:
        for name, digest in plan['source_hashes'].items():
            if sha(root/'source'/name) != digest:
                raise ValueError('changed frozen source')
        binary = root/'smooth_parallel_census'
        if sha(binary) != plan['binary_sha256']:
            raise ValueError('changed frozen binary')
        if shutil.which('caffeinate'):
            caffeinate = subprocess.Popen(['caffeinate', '-i', '-w', str(os.getpid())])
        state(root, status='enumerating', pid=os.getpid(), workers=plan['workers'])
        folder = root/'degree_9'
        with (root/'degree_9.log').open('x') as log:
            child = subprocess.Popen([str(binary), '--d', '9', '--workers', str(plan['workers']),
                    '--out', str(folder), '--seconds', str(plan['seconds']), '--cache-entries', '1000000',
                    '--cap', '100000', '--closure-seconds', '30'], stdout=log, stderr=subprocess.STDOUT)
            while child.poll() is None:
                if shutil.disk_usage(root).free < 30*2**30:
                    child.terminate(); child.wait()
                    raise RuntimeError('disk guard: stopped at 30GiB; partial frame files preserved')
                time.sleep(2)
            if child.returncode:
                raise RuntimeError(f'producer exit {child.returncode}')
        print('ENUMERATED', json.dumps(rows(folder/'COVERAGE.tsv')[0]), flush=True)
        state(root, status='auditing', pid=os.getpid(), workers=plan['workers'])
        checks = audit_parallel(folder, plan['workers'])
        save(root/'CHECKS.json', checks)
        complete = checks['exhaustive'] and checks['coverage']['fully_resolved'] == 'true'
        save(root/'FINISHED.json', dict(status='complete' if complete else 'incomplete', ended=now()))
        state(root, status='complete' if complete else 'incomplete', workers=0)
    except BaseException as error:
        if child is not None and child.poll() is None:
            child.terminate(); child.wait()
        save(root/'FAILURE.json', dict(error=repr(error), ended=now()))
        state(root, status='failed', workers=0, error=repr(error))
        raise
    finally:
        if caffeinate is not None and caffeinate.poll() is None:
            caffeinate.terminate()
        lock.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true'); p.add_argument('--workers', type=int, default=6)
    p.add_argument('--seconds', type=int, default=31536000)
    p.add_argument('--worker', type=Path); p.add_argument('--lock-fd', type=int)
    p.add_argument('--audit', type=Path)
    args = p.parse_args()
    if args.worker:
        worker(args.worker, args.lock_fd); return
    if args.audit:
        print(json.dumps(audit_parallel(args.audit, args.workers), indent=2)); return
    if not 1 <= args.workers <= 10 or args.seconds <= 0:
        p.error('positive seconds and 1..10 workers required')
    if not args.run:
        print(json.dumps(dict(degree=9, workers=args.workers, population='all passports/all transitive groups',
                    quotient='simultaneous S9 conjugation', gpu=False,
                    ram='shared permutation/inverse tables, bounded million-entry proof cache, buffered binary outcomes',
                    completion='independent full reenumeration and every literal rejection replay; unresolved stays unresolved'), indent=2)); return
    dest = Path('data')
    if not Path('data').is_dir() or shutil.disk_usage(dest).free < 30*2**30:
        raise RuntimeError('mounted T7 with 30GiB free required')
    control = PROJECT/'.run-control'; control.mkdir(exist_ok=True)
    lock = (control/'.campaign.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = dest/('smooth_degree9_'+dt.datetime.now(ZoneInfo('America/Chicago')).strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir(); source = root/'source'; source.mkdir(); hashes = {}
    paths = [PROJECT/'Cargo.toml', PROJECT/'Cargo.lock', *sorted((PROJECT/'src').rglob('*.rs')),
        *[_scripts_path.find(s) for s in ['run_smooth_degree9.py', 'test_smooth_degree9.py',
          'run_smooth_small_census.py', 'test_smooth_small_census.py', 'run_smooth_mission.py',
          'run_s18_linked_roots.py', 'analyze_fixed_braid_constraints.py', 'odd_order_spin.py']],
        *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
        PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv',
        PROJECT/'evidence/smooth_small_census_20261011/MATH.md',
        PROJECT/'evidence/smooth_degree9_20261012/PLAN.md']
    for path in paths:
        rel = path.relative_to(PROJECT); target = source/rel; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target); hashes[str(rel)] = sha(target)
    env = dict(os.environ, CARGO_BUILD_JOBS='1', RUSTFLAGS='-C target-cpu=native',
               CARGO_TARGET_DIR=str(PROJECT/'target/smooth-parallel'))
    print('OUTPUT='+str(root), flush=True)
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo', 'test', '--release', '--bin', 'smooth_parallel_census', '_tests',
                        '--', '--test-threads=1'], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cargo', 'build', '--release', '--bin', 'smooth_parallel_census', '--bin', 'smooth_small_census'],
                       cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, str(source/'scripts/smooth/test_smooth_degree9.py')],
                       cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = root/'smooth_parallel_census'
    shutil.copy2(PROJECT/'target/smooth-parallel/release/smooth_parallel_census', binary)
    validation = root/'validation'; validation.mkdir()
    result = verify_small(binary, PROJECT/'target/smooth-parallel/release/smooth_small_census', validation)
    print('VALIDATION='+json.dumps(result), flush=True)
    save(root/'PLAN.json', dict(degree=9, workers=args.workers, seconds=args.seconds, source_hashes=hashes,
          binary_sha256=sha(binary), quotient='simultaneous S9', identity_entries='included for implication, excluded only for nontrivial existence',
          population='all geodesic first-pair frames/all 9! thirds; all passports/genera/transitive groups',
          record_format='per-frame <II little endian: third rank, witness token; literal dictionary; COMPLETE after full frame',
          cache_entries=1000000, cap=100000, closure_seconds=30, started=now()))
    with (root/'WORKER.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(source/'scripts/smooth/run_smooth_degree9.py'), '--worker', str(root),
                '--lock-fd', str(lock.fileno())], cwd=source, stdout=log, stderr=subprocess.STDOUT,
                start_new_session=True, pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json', dict(pid=child.pid, started=now()))
    save(control/f'SMOOTH_DEGREE9_{root.name}.json', dict(output=str(root), pid=child.pid, workers=args.workers))
    lock.close(); print('PID='+str(child.pid), flush=True)


if __name__ == '__main__':
    main()
