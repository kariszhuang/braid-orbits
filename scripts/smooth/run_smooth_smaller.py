#!/usr/bin/env python3
"""Exact degree12 arbitrary-S4-voltage family; plan-only unless --run."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import datetime as dt
import fcntl
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_smooth_mission import audit_component, rows, save, sha
from run_smooth_special import key_digest, known_folders, KNOWN_ROOT
from run_s18_linked_roots import state_key
from verify_smooth_geometry import (compose, inverse, cycles, orbits, hurwitz,
                                   parse_tuple, recover_blocks, require)

PROJECT = Path(__file__).resolve().parents[2]
BETA = (1, 2, 0)
GAMMA = BETA


def orbit(p):
    result, x = [], 0
    while x not in result:
        result.append(x)
        x = p[x]
    require(x == 0 and len(result) == len(p), 'not a full block cycle')
    return result


def lift(voltage, tau):
    a = tuple(tau) + tuple(range(4, 12))
    b = tuple(4*BETA[j]+u for j in range(3) for u in range(4))
    c = tuple(4*GAMMA[j]+voltage[j][u] for j in range(3) for u in range(4))
    return a, b, c, inverse(compose(compose(a, b), c))


def construct(free):
    x, y = free
    voltage = (x, y, inverse(compose(y, x)))
    total = tuple(range(4))
    for j in orbit(compose(BETA, GAMMA)):
        total = compose(voltage[j], total)
    return lift(voltage, inverse(total))


def analyze(t):
    # Separate literal verifier, including relabeling-invariant block recovery.
    from verify_smooth_geometry import analyze as geometric_audit
    out = geometric_audit(t, 3, 4)
    require(out['source_genus'] == 2 and out['quotient_genus'] == 1, 'genus')
    require(out['geometry'] == 'isotrivial_elliptic', 'not elliptic mechanism')
    out['geometry_review'] = 'human_review_pending'
    return out


def finite_population():
    local = list(itertools.permutations(range(4)))
    return [construct((local[rank % 24], local[rank // 24])) for rank in range(576)]


def direct_population():
    # Direct cycle-type checks on every unrestricted three-voltage/local-entry frame.
    local = list(itertools.permutations(range(4)))
    out, examined = set(), 0
    for voltage in itertools.product(local, repeat=3):
        for tau in local:
            examined += 1
            t = lift(voltage, tau)
            if sorted(map(len, cycles(t[2]))) == [3]*4 and sorted(map(len, cycles(t[3]))) == [3]*4:
                out.add(t)
    return examined, out


def audit(root):
    started = time.monotonic()
    population = finite_population()
    frames = rows(root/'job/FRAMES.tsv')
    seeds = rows(root/'job/SEEDS.tsv')
    resolutions = rows(root/'job/RESOLUTIONS.tsv')
    coverage = rows(root/'job/CONSTRUCTION.tsv')[0]
    require(len(frames) == len(population) == 576, 'construction coverage gap')
    require(len(seeds) == len(resolutions), 'missing seed resolution')
    known_states, masses = {}, Counter()
    for rank, (t, row) in enumerate(zip(population, frames)):
        require(int(row['rank']) == rank, 'rank ordering')
        status = ('identity' if t[0] == tuple(range(12)) else
                  'intransitive' if len(orbits(t)) != 1 else 'qualifying')
        require(row['status'] == status, 'wrong frame rejection')
        masses[status] += 1
        if status != 'qualifying':
            continue
        analyze(t)
        for m in range(6):
            analyze(hurwitz(t, m))
        key = state_key(t)
        ident = int(row['candidate'])
        require(known_states.setdefault(key, ident) == ident, 'canonical alias mismatch')
        require(parse_tuple(seeds[ident]['seed'], 12) == population[int(seeds[ident]['rank'])],
                'wrong saved construction')
        require(state_key(parse_tuple(seeds[ident]['seed'], 12)) == key, 'wrong seed alias')
    require(len(known_states) == len(seeds) == int(coverage['unique_seed_states']), 'seed count')
    for name in ['identity', 'intransitive', 'qualifying']:
        require(masses[name] == int(coverage[name]), 'mass mismatch '+name)
    from sympy.combinatorics import Permutation, PermutationGroup
    quotient_order = int(PermutationGroup([Permutation(list(BETA)), Permutation(list(GAMMA))]).order())
    require(quotient_order == 3, 'quotient is not C3')
    folders, checked_states, closed_owner = {}, 0, {}
    for folder in sorted([*root.glob('job/component_*'), *root.glob('job/partial_*')]):
        check = audit_component(folder)
        seed = parse_tuple(check['summary']['seed'], 12)
        geometry = analyze(seed)
        qorder = int(PermutationGroup([Permutation(list(q)) for q in geometry['quotient_permutations']]).order())
        require(qorder == 3, 'saved quotient group')
        order = int(PermutationGroup([Permutation(list(p)) for p in seed]).order())
        require(order == int(check['summary']['group_order']), 'group order mismatch')
        for row in rows(folder/'STATES.tsv'):
            analyze(parse_tuple(row['tuple'], 12))
            checked_states += 1
        check['component_id'] = key_digest(check['state_keys']) if check['independent_six_move_closure'] else None
        check['independent_group_order'] = order
        check['geometry'] = geometry
        folders[folder.name] = check
        if check['independent_six_move_closure']:
            for key in check['state_keys']:
                require(key not in closed_owner, 'overlapping supposedly distinct closed components')
                closed_owner[key] = folder.name
    for seed, res in zip(seeds, resolutions):
        require(seed['candidate'] == res['candidate'], 'resolution ordering')
        if res['status'] != 'deferred':
            check = folders[res['folder']]
            require(state_key(parse_tuple(seed['seed'], 12)) in set(check['state_keys']), 'missing resolved seed')
            if res['status'] == 'alias':
                require(check['independent_six_move_closure'], 'alias to partial component')
            else:
                require(res['status'] == check['summary']['status'], 'wrong resolution status')
                require(int(res['states']) == int(check['summary']['states']), 'wrong state count')
    known = {c['component_id'] for c in json.loads((root/'KNOWN_COMPONENTS.json').read_text())}
    catalog, partial = [], []
    for name, check in folders.items():
        item = dict(folder='job/'+name, summary=check['summary'], geometry=check['geometry'],
                    independent_group_order=check['independent_group_order'], table_sha256=check['table_sha256'])
        if check['component_id'] is not None:
            item.update(component_id=check['component_id'], known_component=check['component_id'] in known)
            catalog.append(item)
        else:
            partial.append(item)
    result = dict(coverage=coverage, independent_occurrence_counts=dict(masses),
                  every_qualifying_occurrence_and_six_neighbors_checked=True,
                  independent_closed_components_pairwise_disjoint=True,
                  every_nondeferred_seed_membership_checked=True,
                  quotient_group_order=quotient_order, saved_states_geometry_checked=checked_states,
                  independent_six_move_closed_components=len(catalog),
                  additional_closed_components=sum(not c['known_component'] for c in catalog),
                  partial_state_counts=len(partial), deferred_seed_states=sum(r['status']=='deferred' for r in resolutions),
                  independent_audit_seconds=time.monotonic()-started,
                  quotient='simultaneous S12 conjugation',
                  scope='All576 fixed-quotient solved frames; all seed resolutions attempted; not a degree census')
    save(root/'CATALOG.json', dict(closed_components=catalog, partial_state_counts=partial))
    save(root/'CHECKS.json', result)
    return result


def worker(root, lock_fd):
    lock = os.fdopen(lock_fd, 'a')
    try:
        plan = json.loads((root/'PLAN.json').read_text())
        for name, digest in plan['source_hashes'].items():
            require(sha(root/'source'/name) == digest, 'changed frozen source '+name)
        require(sha(root/'smooth_smaller') == plan['binary_sha256'], 'changed binary')
        require(sha(root/'KNOWN_COMPONENTS.json') == plan['known_sha256'], 'changed known catalog')
        save(root/'STATE.json', dict(status='running', pid=os.getpid(), started=dt.datetime.now().isoformat()))
        cmd = [str(root/'smooth_smaller'), '--out', str(root/'job'), '--cap', str(plan['cap']),
               '--closure-seconds', str(plan['closure_seconds']), '--closures', str(plan['closures'])]
        subprocess.run(cmd, check=True)
        check = audit(root)
        save(root/'FINISHED.json', dict(status='construction_complete_selected_audit_complete',
             ended=dt.datetime.now().isoformat(), partial_state_counts=check['partial_state_counts'],
             deferred_seed_states=check['deferred_seed_states']))
        print(json.dumps(check, indent=2), flush=True)
    except Exception as error:
        save(root/'FAILURE.json', dict(error=repr(error)))
        raise
    finally:
        lock.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true')
    p.add_argument('--closures', type=int, default=100)
    p.add_argument('--closure-seconds', type=int, default=5)
    p.add_argument('--cap', type=int, default=10000)
    p.add_argument('--worker', type=Path)
    p.add_argument('--lock-fd', type=int)
    args = p.parse_args()
    if args.worker:
        worker(args.worker, args.lock_fd)
        return
    if min(args.closures, args.closure_seconds, args.cap) <= 0:
        p.error('positive limits required')
    if not args.run:
        print(json.dumps(dict(degree=12, blocks=3, size=4, solved_choices=576,
              direct_same_family_choices=331776, quotient_group='C3', cpu_workers=1, gpu=False,
              closures=args.closures, closure_seconds=args.closure_seconds, cap=args.cap), indent=2))
        return
    import sympy
    dest = Path('data')
    require(Path('data').is_dir() and shutil.disk_usage(dest).free >= 30*2**30,
            'mounted T7 with30GiB free required')
    lock = (PROJECT/'.run-control/.campaign.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = dest/('smooth_smaller_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir()
    print('OUTPUT='+str(root), flush=True)
    source, hashes = root/'source', {}
    paths = [PROJECT/'Cargo.toml', PROJECT/'Cargo.lock', *sorted((PROJECT/'src').rglob('*.rs')),
             *sorted((PROJECT/'scripts').rglob('*.py')),
             *sorted((PROJECT/'evidence/smooth_smaller_20261012').glob('*.md')),
             PROJECT/'evidence/smooth_geometry_20261012/PROOF.md',
             PROJECT/'evidence/smooth_geometry_20261012/SOURCES.json',
             PROJECT/'evidence/smooth_special_20261012/CATALOG.json',
             PROJECT/'evidence/smooth_noncyclic_20261012/CATALOG.json',
             *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
             PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv']
    for path in paths:
        rel = path.relative_to(PROJECT)
        target = source/rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(rel)] = sha(target)
    env = dict(os.environ, CARGO_BUILD_JOBS='1', CARGO_TARGET_DIR=str(PROJECT/'target/smooth-smaller'),
               RUSTFLAGS='-C target-cpu=native')
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo', 'test', '--release', '--bin', 'smooth_smaller', 'smaller_tests', '--', '--test-threads=1'],
                       cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cargo', 'build', '--release', '--bin', 'smooth_smaller'],
                       cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, 'scripts/smooth/test_smooth_smaller.py'], cwd=source, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = root/'smooth_smaller'
    shutil.copy2(PROJECT/'target/smooth-smaller/release/smooth_smaller', binary)
    prior = json.loads((source/'evidence/smooth_special_20261012/CATALOG.json').read_text())
    require(len(prior['components']) == 52, 'known special catalog boundary changed')
    # Only independently audited, CLOSED state-set identities may exclude novelty.
    known = [dict(component_id=x['component_id'], folder=x['folder'], table_sha256=x['table_sha256'])
             for x in prior['components']]
    for item in known:
        require(sha(Path(prior['raw_root'])/item['folder']/'STATES.tsv') == item['table_sha256'],
                'changed prior closed state table')
    for folder in known_folders(KNOWN_ROOT):
        check = audit_component(folder)
        require(check['independent_six_move_closure'], 'prior controls must be closed')
        known.append(dict(component_id=key_digest(check['state_keys']), folder=str(folder),
                          table_sha256=check['table_sha256']))
    previous = json.loads((source/'evidence/smooth_noncyclic_20261012/CATALOG.json').read_text())
    for item in previous['closed_components']:
        require(sha(Path(previous['raw_root'])/item['folder']/'STATES.tsv') == item['table_sha256'],
                'changed previous noncyclic state table')
        known.append(dict(component_id=item['component_id'], folder=item['folder'], table_sha256=item['table_sha256']))
    save(root/'KNOWN_COMPONENTS.json', known)
    save(root/'PLAN.json', dict(source_hashes=hashes, binary_sha256=sha(binary),
         known_sha256=sha(root/'KNOWN_COMPONENTS.json'), sympy_version=sympy.__version__,
         closures=args.closures, closure_seconds=args.closure_seconds, cap=args.cap,
         cpu_workers=1, gpu=False, quotient='simultaneous S12 conjugation',
         population='fixed beta=gamma=(123),24^2 free S4 voltages',
         objective='smaller positive-genus smooth examples; isotrivial, proper groups, NOT weak failures'))
    with (root/'WORKER.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(source/'scripts/smooth/run_smooth_smaller.py'),
                '--worker', str(root), '--lock-fd', str(lock.fileno())], cwd=source,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(lock.fileno(),))
    save(root/'LAUNCH.json', dict(pid=child.pid, started=dt.datetime.now().isoformat()))
    lock.close()
    print('PID='+str(child.pid), flush=True)
    if shutil.which('caffeinate'):
        subprocess.Popen(['caffeinate', '-i', '-w', str(child.pid)], stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL, start_new_session=True)


if __name__ == '__main__':
    main()
