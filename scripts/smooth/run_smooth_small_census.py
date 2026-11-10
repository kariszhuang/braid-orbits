#!/usr/bin/env python3
"""Shared-lock, frozen exhaustive smooth population in degrees 5--8."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import fcntl
import itertools
import json
from math import factorial
import os
from pathlib import Path
import shutil
import subprocess
import sys

from analyze_fixed_braid_constraints import product, point_orbits
from odd_order_spin import inverse, hurwitz, parse_seed, cycle_lengths
from run_smooth_mission import audit_component, rows, save, sha, values
from run_s18_linked_roots import state_key

PROJECT = Path(__file__).resolve().parents[2]


def partition(n, minimum=1):
    if not n:
        yield ()
    for k in range(minimum, n+1):
        for tail in partition(n-k, k):
            yield (k, *tail)


def pair_key_and_centralizer(a, b):
    # A separate tuple-image implementation, independent of Rust's nibble keys.
    from collections import Counter
    components = []
    for points in point_orbits((a, b)):
        rooted = []
        for root in points:
            order = [root]; labels = {root: 0}
            for x in order:
                for p in (a, b):
                    y = p[x]
                    if y not in labels:
                        labels[y] = len(order); order.append(y)
            rooted.append(tuple(tuple(labels[p[x]] for x in order) for p in (a, b)))
        best = min(rooted)
        components.append((len(points), best, rooted.count(best)))
    components.sort()
    z = 1
    for _, _, aut in components:
        z *= aut
    for m in Counter((size, key) for size, key, _ in components).values():
        z *= factorial(m)
    return tuple((size, key) for size, key, _ in components), z


def audit_frames(folder, n, perms):
    frames = rows(folder/'FRAMES.tsv'); keys = set(); result = []
    for i, row in enumerate(frames):
        if int(row['frame']) != i:
            raise ValueError('frame numbering')
        a, b, _, _ = parse_seed(f"[{row['a']}, {row['b']}, (), ()]", n)
        ip = n-len(cycle_lengths(product(a, b)))
        if 2*n-len(cycle_lengths(a))-len(cycle_lengths(b)) != ip:
            raise ValueError('nongeodesic frame')
        key, z = pair_key_and_centralizer(a, b)
        if key in keys or factorial(n)//z != int(row['pair_orbit']):
            raise ValueError('duplicate or wrong-mass frame')
        keys.add(key); result.append((a, b, factorial(n)//z))
    expected = set(); mass = 0
    for lens in partition(n):
        p = list(range(n)); off = 0
        for length in lens:
            for x in range(length):
                p[off+x] = off+(x+1)%length
            off += length
        p = tuple(p); ip = n-len(lens); raw = 0
        for a, ia in perms:
            b = product(inverse(a), p)
            if ia+n-len(cycle_lengths(b)) == ip:
                key, z = pair_key_and_centralizer(a, b)
                expected.add(key); raw += 1
        from collections import Counter
        z = 1
        for l, multiplicity in Counter(lens).items():
            z *= l**multiplicity*factorial(multiplicity)
        mass += raw*factorial(n)//z
    if keys != expected or mass != sum(w for _, _, w in result):
        raise ValueError('incomplete frame population')
    return result


def audit(folder):
    coverage = rows(folder/'COVERAGE.tsv')[0]; n = int(coverage['degree'])
    if coverage['enumeration_complete'] != 'true':
        return dict(coverage=coverage, exhaustive=False, warning='Incomplete finite prefix, not a census.')
    perms = [(p, n-len(cycle_lengths(p))) for p in itertools.permutations(range(n))]
    frames = audit_frames(folder, n, perms)
    counts = rows(folder/'FRAME_COUNTS.tsv')
    if len(counts) != len(frames):
        raise ValueError('missing frame counts')
    components = {p.name: audit_component(p) for p in sorted(folder.glob('component_*'))}
    closed_keys = {name: set(c['state_keys']) for name, c in components.items()}
    if sum(map(len, closed_keys.values())) != len(set().union(*closed_keys.values())):
        raise ValueError('duplicate smooth components')
    for c in components.values():
        c.pop('state_keys')
    totals = [0]*8; weighted = [0, 0]; outcomes = 0
    with (folder/'OUTCOMES.tsv').open() as stream:
        import csv
        records = iter(csv.DictReader(stream, delimiter='\t'))
        for fi, (a, b, weight) in enumerate(frames):
            p = product(a, b); ip = n-len(cycle_lengths(p)); count = [0]*8
            left = len(point_orbits((a, b)))
            for rank, (c, ic) in enumerate(perms):
                count[0] += 1
                d = inverse(product(p, c)); t = (a, b, c, d)
                if len(point_orbits(t[:3])) != 1:
                    continue
                count[1] += 1
                right = len(point_orbits(t[2:]))
                if n-ip-left-right+1 != 0 or right != 1:
                    raise ValueError('independent implication violation')
                if ic+n-len(cycle_lengths(d)) == ip:
                    count[2] += 1
                if any(q == tuple(range(n)) for q in t):
                    count[4] += 1; continue
                count[3] += 1; row = next(records); outcomes += 1
                if (int(row['frame']), int(row['third_rank'])) != (fi, rank):
                    raise ValueError('missing, extra or reordered outcome')
                status = row['status']; data = row['word_or_component']
                if status == 'nonsmooth':
                    for move in data.split(','):
                        if move:
                            t = hurwitz(t, int(move))
                    if values(t)['smooth']:
                        raise ValueError('false literal nonsmooth witness')
                    count[5] += 1
                elif status == 'closed_smooth':
                    if state_key(t) not in closed_keys[data]:
                        raise ValueError('false closed membership')
                    count[6] += 1
                else:
                    count[7] += 1
                    audit_component(folder/data)
            expected = [int(counts[fi][name]) for name in ['enumerated', 'transitive', 'both_geodesic', 'nonidentity', 'identity', 'nonsmooth_orbit', 'closed_occurrences', 'unresolved']]
            if count != expected or int(counts[fi]['frame']) != fi:
                raise ValueError('per-frame census mismatch')
            totals = [x+y for x, y in zip(totals, count)]
            weighted[0] += weight*count[1]; weighted[1] += weight*count[2]
        if next(records, None) is not None:
            raise ValueError('extra outcome rows')
    names = ['enumerated', 'transitive_smooth_occurrences', 'both_geodesic_occurrences', 'nonidentity_occurrences', 'identity_occurrences', 'nonsmooth_orbit_occurrences', 'closed_occurrences', 'unresolved']
    if totals != [int(coverage[x]) for x in names] or len(components) != int(coverage['closed_components']):
        raise ValueError('coverage totals mismatch')
    if int(coverage['labeled_left_smooth']) != weighted[0] or int(coverage['labeled_both_smooth']) != weighted[1] or int(coverage['labeled_smooth_tuples']) != 2*weighted[0]-weighted[1]:
        raise ValueError('labeled mass mismatch')
    return dict(coverage=coverage, exhaustive=True, independently_reenumerated_third_occurrences=totals[0],
                independently_tested_smooth_occurrences=totals[1], literal_orbit_outcomes_checked=outcomes,
                components=components, sources={p.name: sha(p) for p in folder.glob('*.tsv')})


def worker(root, lock_fd):
    lock = os.fdopen(lock_fd, 'a'); plan = json.loads((root/'PLAN.json').read_text())
    try:
        for name, digest in plan['source_hashes'].items():
            if sha(root/'source'/name) != digest:
                raise ValueError('changed frozen source')
        binary = root/'smooth_small_census'
        if sha(binary) != plan['binary_sha256']:
            raise ValueError('changed frozen binary')
        save(root/'STATE.json', dict(status='running', pid=os.getpid(), started=dt.datetime.now().isoformat()))
        def job(n):
            folder = root/f'degree_{n}'
            with (root/f'degree_{n}.log').open('x') as log:
                subprocess.run([str(binary), '--d', str(n), '--out', str(folder), '--seconds', str(plan['seconds']), '--cap', '50000'], stdout=log, stderr=subprocess.STDOUT, check=True)
            print('ENUMERATED', n, json.dumps(rows(folder/'COVERAGE.tsv')[0]), flush=True)
            result = audit(folder)
            save(root/f'degree_{n}_CHECKS.json', result)
            print('INDEPENDENTLY_AUDITED', n, flush=True)
            return result
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            results = list(pool.map(job, plan['degrees']))
        complete = all(r['exhaustive'] and r['coverage']['fully_resolved'] == 'true' for r in results)
        save(root/'CHECKS.json', dict(jobs=results, fully_resolved=complete))
        save(root/'FINISHED.json', dict(status='complete' if complete else 'incomplete', ended=dt.datetime.now().isoformat()))
    except Exception as error:
        save(root/'FAILURE.json', dict(error=repr(error))); raise
    finally:
        lock.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true'); p.add_argument('--seconds', type=int, default=3600)
    p.add_argument('--workers', type=int, default=4); p.add_argument('--degrees', type=int, nargs='+', default=[5, 6, 7, 8])
    p.add_argument('--worker', type=Path); p.add_argument('--lock-fd', type=int)
    p.add_argument('--audit', type=Path); p.add_argument('--audit-out', type=Path)
    args = p.parse_args()
    if args.worker:
        worker(args.worker, args.lock_fd); return
    if args.audit:
        if not args.audit_out:
            p.error('--audit-out required')
        save(args.audit_out, audit(args.audit)); return
    if not args.run:
        print(json.dumps(dict(degrees=args.degrees, workers=args.workers, seconds=args.seconds, exhaustive_smooth_population=True), indent=2)); return
    if not 1 <= args.workers <= 7 or args.seconds <= 0 or any(not 3 <= n <= 9 for n in args.degrees) or len(set(args.degrees)) != len(args.degrees):
        p.error('valid distinct degrees, 1..7 workers and positive seconds required')
    dest = Path('data')
    if not Path('data').is_dir() or shutil.disk_usage(dest).free < 30*2**30:
        raise RuntimeError('mounted T7 with 30 GiB free required')
    lock = (PROJECT/'.run-control/.campaign.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = dest/('smooth_small_census_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')); root.mkdir()
    source = root/'source'; source.mkdir(); hashes = {}
    paths = [PROJECT/'Cargo.toml', PROJECT/'Cargo.lock', *sorted((PROJECT/'src').rglob('*.rs')),
             *[_scripts_path.find(s) for s in ['run_smooth_small_census.py', 'test_smooth_small_census.py', 'run_smooth_mission.py', 'run_s18_linked_roots.py', 'analyze_fixed_braid_constraints.py', 'odd_order_spin.py']],
             *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
             PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv',
             PROJECT/'evidence/smooth_small_census_20261011/PLAN.md',
             PROJECT/'evidence/smooth_small_census_20261011/MATH.md']
    for path in paths:
        rel = path.relative_to(PROJECT); target = source/rel; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target); hashes[str(rel)] = sha(target)
    env = dict(os.environ, CARGO_BUILD_JOBS='1', CARGO_TARGET_DIR=str(PROJECT/'target/smooth-small'))
    print('OUTPUT='+str(root), flush=True)
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo', 'test', '--release', '--bin', 'smooth_small_census', '_tests', '--', '--test-threads=1'], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cargo', 'build', '--release', '--bin', 'smooth_small_census'], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = root/'smooth_small_census'; shutil.copy2(PROJECT/'target/smooth-small/release/smooth_small_census', binary)
    save(root/'PLAN.json', dict(degrees=args.degrees, workers=args.workers, seconds=args.seconds, source_hashes=hashes, binary_sha256=sha(binary), quotient='exact simultaneous S_n first-pair reduction; full n! third enumeration', identity_entries='included for implication, separated from nontrivial orbit search'))
    with (root/'WORKER.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(source/'scripts/smooth/run_smooth_small_census.py'), '--worker', str(root), '--lock-fd', str(lock.fileno())], cwd=source, stdout=log, stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json', dict(pid=child.pid, started=dt.datetime.now().isoformat()))
    lock.close(); print('PID='+str(child.pid), flush=True)


if __name__ == '__main__':
    main()
