#!/usr/bin/env python3
"""Frozen, shared-lock proof-directed special smoothness mission."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import fcntl
import hashlib
import itertools
import json
from math import factorial, gcd
import os
from pathlib import Path
import shutil
import subprocess
import sys

from analyze_fixed_braid_constraints import point_orbits, product
from odd_order_spin import cycle_lengths, inverse, parse_seed
from run_smooth_mission import audit_component, rows, save, sha, values

PROJECT = Path(__file__).resolve().parents[2]
KNOWN_ROOT = Path('data/smooth_mission_20261011_231953_552478')
JOBS = sorted([dict(blocks=r, size=k) for r in (3, 5, 7) for k in range(3, 9)
               if r*k <= 25 and (r, k) != (3, 3)], key=lambda j: (j['blocks']*j['size'], j['blocks']))


def known_folders(root):
    candidates = [root/'reported_A9', *sorted(root.glob('A*/component_*')), *sorted(root.glob('block_family_*'))]
    return [path for path in candidates if path.is_dir() and (path/'SUMMARY.tsv').is_file()]


def construct(r, k, q, y):
    """Independent tuple-image construction, not Rust's packed Perm engine."""
    x = tuple((i+1) % k for i in range(k))
    voltage = [tuple(range(k)) for _ in range(r)]
    voltage[0], voltage[q], voltage[-q] = x, y, inverse(product(y, x))
    positions = [(j*(q+1)) % r for j in range(r)]
    # product acts rightmost first: traverse returns in reversed point order.
    tau = inverse(product(*(voltage[j] for j in reversed(positions))))
    a = tuple(tau) + tuple(range(k, r*k))
    b = tuple(((j+1) % r)*k+i for j in range(r) for i in range(k))
    c = tuple(((j+q) % r)*k+voltage[j][i] for j in range(r) for i in range(k))
    return a, b, c, inverse(product(a, b, c))


def local_type(p):
    return tuple(sorted(l for l in cycle_lengths(p) if l > 1))


def recover_blocks(t, r, k):
    """A relabeling-invariant certificate of the transversal hypothesis."""
    types = [local_type(p) for p in t]
    small = [i for i, typ in enumerate(types) if typ != (r,)*k]
    if len(small) != 1:
        return None
    small = small[0]; large = [i for i in range(4) if i != small]
    n = r*k
    for q in range(1, r):
        if gcd(q, r) != 1 or gcd(q+1, r) != 1:
            continue
        for offsets in set(itertools.permutations((1, q, (-q-1) % r))):
            shifts = [0]*4
            for j, offset in zip(large, offsets):
                shifts[j] = offset
            colors = {0: 0}; queue = [0]; valid = True
            for x in queue:
                for p, shift in zip(t, shifts):
                    y = p[x]; expected = (colors[x]+shift) % r
                    if y in colors and colors[y] != expected:
                        valid = False; break
                    if y not in colors:
                        colors[y] = expected; queue.append(y)
                if not valid:
                    break
            if not valid or len(colors) != n or sorted(Counter(colors.values()).values()) != [k]*r:
                continue
            support = {colors[x] for x in range(n) if t[small][x] != x}
            if len(support) == 1:
                return dict(small_position=small+1, shifts=shifts, colors=[colors[x] for x in range(n)])
    return None


def key_digest(keys):
    return hashlib.sha256(json.dumps(sorted(keys), separators=(',', ':')).encode()).hexdigest()


def audit_job(folder, job, representatives):
    r, k = job['blocks'], job['size']; n = r*k
    coverage = rows(folder/'COVERAGE.tsv')[0]
    inputs = rows(folder/'CANDIDATES.tsv'); result_rows = rows(folder/'RESOLUTIONS.tsv')
    if len(inputs) != len(result_rows) or len(inputs) != int(coverage['selected']):
        raise ValueError('missing selected-candidate resolution')
    independent = Counter(); identity = old = examined = 0
    qs = [q for q in range(1, r) if gcd(q, r) == gcd(q+1, r) == 1]
    # Count the finite voltage population independently, on k points rather
    # than canonicalizing rk-point tuples. No closure/generation per occurrence.
    x = tuple((i+1) % k for i in range(k))
    for q in qs:
        for y in itertools.permutations(range(k)):
            voltage = [tuple(range(k)) for _ in range(r)]
            voltage[0], voltage[q], voltage[-q] = x, y, inverse(product(y, x))
            positions = [(j*(q+1)) % r for j in range(r)]
            tau = inverse(product(*(voltage[j] for j in reversed(positions))))
            typ = local_type(tau); examined += 1
            if not typ:
                identity += 1
            elif r == 3 and typ == (3,):
                old += 1
            else:
                independent[q, str(list(typ))] += 1
    census_complete = coverage['enumeration_complete'] == 'true'
    if census_complete:
        recorded = {(int(row['q']), row['small_type']): int(row['occurrences']) for row in rows(folder/'TYPES.tsv')}
        if dict(independent) != recorded or examined != int(coverage['enumerated']) or identity != int(coverage['identity']) or old != int(coverage['old_template']):
            raise ValueError('independent construction mass mismatch')
    checks = []; folders = {}; proof_count = 0
    permutations = list(itertools.permutations(range(k)))
    for row, resolution in zip(inputs, result_rows):
        if row['candidate'] != resolution['candidate']:
            raise ValueError('misordered candidate resolution')
        q, rank = int(row['q']), int(row['y_rank'])
        t = parse_seed(row['seed'], n)
        if t != construct(r, k, q, permutations[rank]):
            raise ValueError('wrong exported literal construction')
        witness = recover_blocks(t, r, k); v = values(t)
        if witness is None or product(*t) != tuple(range(n)) or len(point_orbits(t)) != 1 or not v['smooth'] or v['defect'] != 0:
            raise ValueError('false theorem input')
        predicted = 1+(k*(r-3)+(n-len(cycle_lengths(t[0]))))//2
        if v['genus'] != predicted or predicted != int(row['genus']):
            raise ValueError('genus prediction mismatch')
        proof_count += 1
        checks.append(dict(candidate=row['candidate'], theorem_hypotheses_checked=True,
                           genus=predicted, seed=row['seed'], coloring=witness, closure_status=resolution['status']))
    for path in sorted([*folder.glob('component_*'), *folder.glob('partial_*')]):
        audit = audit_component(path); summary = audit['summary']
        for row in rows(path/'STATES.tsv'):
            if recover_blocks(parse_seed(row['tuple'], n), r, k) is None:
                raise ValueError('saved state lost theorem hypothesis')
        from sympy.combinatorics import Permutation, PermutationGroup
        t = parse_seed(summary['seed'], n)
        order = int(PermutationGroup([Permutation(list(p)) for p in t]).order())
        if order != int(summary['group_order']):
            raise ValueError('independent exact group order mismatch')
        audit['component_id'] = key_digest(audit['state_keys']) if audit['independent_six_move_closure'] else None
        audit['independent_group_order'] = order
        folders[path.name] = audit
    from run_s18_linked_roots import state_key
    for row, resolution in zip(inputs, result_rows):
        if resolution['status'] == 'alias':
            audit = folders[resolution['folder']]
            if not audit['independent_six_move_closure'] or state_key(parse_seed(row['seed'], n)) not in set(audit['state_keys']):
                raise ValueError('false closed-state alias')
    return dict(coverage=coverage, independent_mass_complete=census_complete,
                independent_enumerated=examined, selected_theorem_inputs_checked=proof_count,
                selected=checks, folders=folders,
                warning='Only selected representatives closed; counts are voltage occurrences, not all braid components. Isotriviality untested.')


def worker(root, lock_fd):
    lock = os.fdopen(lock_fd, 'a'); plan = json.loads((root/'PLAN.json').read_text())
    try:
        for name, digest in plan['source_hashes'].items():
            if sha(root/'source'/name) != digest:
                raise ValueError('changed frozen source: '+name)
        if sha(root/'smooth_special') != plan['binary_sha256'] or sha(root/'KNOWN_COMPONENTS.json') != plan['known_sha256']:
            raise ValueError('changed binary or known component catalog')
        save(root/'STATE.json', dict(status='running', pid=os.getpid(), started=dt.datetime.now().isoformat()))
        def job(j):
            name = f"r{j['blocks']}_k{j['size']}"; folder = root/name
            cmd = [str(root/'smooth_special'), '--blocks', str(j['blocks']), '--size', str(j['size']),
                   '--out', str(folder), '--seconds', str(plan['seconds']),
                   '--closure-seconds', str(plan['closure_seconds']), '--cap', str(plan['cap']),
                   '--representatives', str(plan['representatives'])]
            with (root/(name+'.log')).open('x') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            print('CONSTRUCTED', name, json.dumps(rows(folder/'COVERAGE.tsv')[0]), flush=True)
            answer = audit_job(folder, j, plan['representatives'])
            save(root/(name+'_CHECKS.json'), answer)
            print('AUDITED', name, flush=True)
            return name, answer
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            answers = list(pool.map(job, plan['jobs']))
        known = {x['component_id'] for x in json.loads((root/'KNOWN_COMPONENTS.json').read_text())}
        unique = {}; partial = []
        for name, answer in answers:
            for folder, check in answer['folders'].items():
                item = dict(folder=name+'/'+folder, summary=check['summary'],
                            independent_group_order=check['independent_group_order'], table_sha256=check['table_sha256'])
                if check['independent_six_move_closure']:
                    ident = check['component_id']; item.update(component_id=ident, known_component=ident in known)
                    unique.setdefault(ident, item)
                else:
                    partial.append(item)
        complete = all(a['independent_mass_complete'] for _, a in answers)
        save(root/'CATALOG.json', dict(distinct_closed_components=list(unique.values()),
             additional_closed_components=sum(not c['known_component'] for c in unique.values()),
             partial_state_counts=partial, known_component_count=len(known),
             warning='Additional means distinct from the prior computed catalog, not new to literature or non-isotrivial.'))
        save(root/'CHECKS.json', dict(jobs=[dict(name=name, coverage=a['coverage'], independent_mass_complete=a['independent_mass_complete']) for name, a in answers],
             construction_complete=complete, selected_state_counts_complete=not partial))
        save(root/'FINISHED.json', dict(status='complete' if complete else 'construction_incomplete',
             ended=dt.datetime.now().isoformat(), partial_state_counts=len(partial)))
    except Exception as error:
        save(root/'FAILURE.json', dict(error=repr(error))); raise
    finally:
        lock.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true'); p.add_argument('--workers', type=int, default=5)
    p.add_argument('--seconds', type=int, default=600); p.add_argument('--closure-seconds', type=int, default=3)
    p.add_argument('--cap', type=int, default=10000); p.add_argument('--representatives', type=int, default=2)
    p.add_argument('--worker', type=Path); p.add_argument('--lock-fd', type=int)
    args = p.parse_args()
    if args.worker:
        worker(args.worker, args.lock_fd); return
    if not 1 <= args.workers <= 7 or min(args.seconds, args.closure_seconds, args.cap, args.representatives) <= 0:
        p.error('1..7 workers and positive limits required')
    if not args.run:
        print(json.dumps(dict(jobs=JOBS, workers=args.workers, seconds=args.seconds,
              objective='new guaranteed-smooth passports; geometric isotriviality untested'), indent=2)); return
    import sympy  # Fail before acquiring/building if the independent engine is unavailable.
    dest = Path('data')
    if not Path('data').is_dir() or shutil.disk_usage(dest).free < 30*2**30:
        raise RuntimeError('mounted T7 and30GiB free required')
    lock = (PROJECT/'.run-control/.campaign.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = dest/('smooth_special_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')); root.mkdir()
    print('OUTPUT='+str(root), flush=True)
    source = root/'source'; source.mkdir(); hashes = {}
    paths = [PROJECT/'Cargo.toml', PROJECT/'Cargo.lock', *sorted((PROJECT/'src').rglob('*.rs')),
             *[_scripts_path.find(s) for s in ['run_smooth_special.py', 'test_smooth_special.py', 'run_smooth_mission.py',
                'run_s18_linked_roots.py', 'odd_order_spin.py', 'analyze_fixed_braid_constraints.py']],
             *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
             PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv',
             PROJECT/'evidence/smooth_special_20261012/PLAN.md', PROJECT/'evidence/smooth_special_20261012/MATH.md']
    for path in paths:
        rel = path.relative_to(PROJECT); target = source/rel; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target); hashes[str(rel)] = sha(target)
    env = dict(os.environ, CARGO_BUILD_JOBS='1', CARGO_TARGET_DIR=str(PROJECT/'target/smooth-special'), RUSTFLAGS='-C target-cpu=native')
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo', 'test', '--release', '--bin', 'smooth_special', '_tests', '--', '--test-threads=1'], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cargo', 'build', '--release', '--bin', 'smooth_special'], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, str(source/'scripts/smooth/test_smooth_special.py')], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = root/'smooth_special'; shutil.copy2(PROJECT/'target/smooth-special/release/smooth_special', binary)
    known = {}
    if not (KNOWN_ROOT/'FINISHED.json').is_file():
        raise ValueError('prior positive control mission must be complete')
    for folder in known_folders(KNOWN_ROOT):
        answer = audit_component(folder)
        if not answer['independent_six_move_closure']:
            raise ValueError('known catalog must contain closed components')
        digest = key_digest(answer['state_keys'])
        known.setdefault(digest, dict(component_id=digest, folder=str(folder), table_sha256=answer['table_sha256']))
    save(root/'KNOWN_COMPONENTS.json', list(known.values()))
    save(root/'PLAN.json', dict(jobs=JOBS, workers=args.workers, seconds=args.seconds,
         closure_seconds=args.closure_seconds, cap=args.cap, representatives=args.representatives,
         source_hashes=hashes, binary_sha256=sha(binary), known_sha256=sha(root/'KNOWN_COMPONENTS.json'),
         sympy_version=sympy.__version__, quotient='simultaneous S_n conjugation',
         objective='proof-directed special smoothness; not weak counterexamples'))
    with (root/'WORKER.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(source/'scripts/smooth/run_smooth_special.py'),
               '--worker', str(root), '--lock-fd', str(lock.fileno())], cwd=source,
               stdout=log, stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json', dict(pid=child.pid, started=dt.datetime.now().isoformat()))
    lock.close(); print('PID='+str(child.pid), flush=True)
    if shutil.which('caffeinate'):
        subprocess.Popen(['caffeinate', '-i', '-w', str(child.pid)], stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL, start_new_session=True)


if __name__ == '__main__':
    main()
