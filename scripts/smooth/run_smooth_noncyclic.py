#!/usr/bin/env python3
"""Exact A5 block-quotient pilot; plan-only unless --run is supplied."""
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
BETA = (1, 2, 3, 4, 0)
GAMMA = (1, 3, 0, 4, 2)


def orbit(p):
    result, x = [], 0
    while x not in result:
        result.append(x)
        x = p[x]
    require(x == 0 and len(result) == len(p), 'not a full block cycle')
    return result


def lift(voltage, tau):
    a = tuple(tau) + tuple(range(3, 15))
    b = tuple(3*BETA[j]+u for j in range(5) for u in range(3))
    c = tuple(3*GAMMA[j]+voltage[j][u] for j in range(5) for u in range(3))
    return a, b, c, inverse(compose(compose(a, b), c))


def construct(free):
    voltage = [tuple(range(3))]*5
    gamma_order = orbit(GAMMA)
    total = tuple(range(3))
    for j, p in zip(gamma_order, free):
        voltage[j] = p
        total = compose(p, total)
    voltage[gamma_order[-1]] = inverse(total)
    total = tuple(range(3))
    for j in orbit(compose(BETA, GAMMA)):
        total = compose(voltage[j], total)
    return lift(voltage, inverse(total))


def analyze(t):
    """Literal, relabeling-invariant geometry check; no cyclic/Kummer test."""
    n, r, k = 15, 5, 3
    ident = tuple(range(n))
    require(all(sorted(p) == list(ident) for p in t), 'invalid permutation')
    require(compose(compose(compose(t[0], t[1]), t[2]), t[3]) == ident, 'product')
    require(len(orbits(t)) == 1, 'intransitive total group')
    support = [{x for x in range(n) if p[x] != x} for p in t]
    small = [i for i, s in enumerate(support) if 0 < len(s) <= k]
    require(len(small) == 1, 'unique local entry')
    small = small[0]
    blocks = recover_blocks(t, support[small], r, k)
    membership = {x: j for j, part in enumerate(blocks) for x in part}
    require(len(membership) == n, 'block partition')
    quotient, relative = [], []
    for i, p in enumerate(t):
        q = []
        for part in blocks:
            image = {membership[p[x]] for x in part}
            require(len(image) == 1, 'not invariant blocks')
            q.append(next(iter(image)))
        if i == small:
            require(q == list(range(r)) and len({membership[x] for x in support[i]}) == 1,
                    'entry not localized')
        else:
            require(sorted(map(len, cycles(q))) == [r], 'large quotient not r-cycle')
            require(sorted(map(len, cycles(p))) == [r]*k, 'large point type not r^k')
        quotient.append(tuple(q))
        fibers = []
        for qc in cycles(q):
            lengths = [len(pc)//len(qc) for pc in cycles(p) if membership[pc[0]] in qc]
            require(all(len(pc) % len(qc) == 0 for pc in cycles(p) if membership[pc[0]] in qc),
                    'relative index not integral')
            require(sum(lengths) == k, 'relative degree')
            fibers.append(lengths)
        relative.append(fibers)
    require(compose(compose(compose(*quotient[:2]), quotient[2]), quotient[3]) == tuple(range(r)),
            'quotient product')
    branch = [sum(any(e > 1 for e in f) for f in fibers) for fibers in relative]
    require(branch == [int(i == small) for i in range(4)], 'unique relative branch value')
    indices = [n-len(cycles(p)) for p in t]
    genus = 1-n+sum(indices)//2
    h = 1-r+sum(r-len(cycles(q)) for q in quotient)//2
    require(h == 2 and genus == 5, 'genus')
    require(2*genus-2 == k*(2*h-2)+indices[small], 'relative Riemann-Hurwitz')
    ip = n-len(cycles(compose(t[0], t[1])))
    defect = n-ip-len(orbits(t[:2]))-len(orbits(t[2:]))+1
    require(ip in [sum(indices[:2]), sum(indices[2:])] and defect == 0, 'smoothness/defect')
    return dict(source_genus=genus, quotient_genus=h, defect=defect, smooth=True,
                small_position=small+1, relative_branch_value_count=sum(branch),
                quotient_permutations=quotient, block_partition=blocks,
                geometry='nonisotrivial by fixed-genus2 quotient theorem; human review pending')


def finite_population():
    local = list(itertools.permutations(range(3)))
    out = []
    for rank in range(1296):
        value, free = rank, []
        for _ in range(4):
            free.append(local[value % 6])
            value //= 6
        out.append(construct(free))
    return out


def direct_population():
    """Direct actual point cycle-type tests, without solved voltage equations."""
    local = list(itertools.permutations(range(3)))
    out, examined = set(), 0
    for voltage in itertools.product(local, repeat=5):
        for tau in local:
            examined += 1
            t = lift(voltage, tau)
            if sorted(map(len, cycles(t[2]))) == [5]*3 and sorted(map(len, cycles(t[3]))) == [5]*3:
                out.add(t)
    return examined, out


def audit(root):
    started = time.monotonic()
    population = finite_population()
    frames = rows(root/'job/FRAMES.tsv')
    seeds = rows(root/'job/SEEDS.tsv')
    resolutions = rows(root/'job/RESOLUTIONS.tsv')
    coverage = rows(root/'job/CONSTRUCTION.tsv')[0]
    require(len(frames) == len(population) == 1296, 'construction coverage gap')
    require(len(seeds) == len(resolutions), 'missing seed resolution')
    known_states, masses = {}, Counter()
    for rank, (t, row) in enumerate(zip(population, frames)):
        require(int(row['rank']) == rank, 'rank ordering')
        status = ('identity' if t[0] == tuple(range(15)) else
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
        require(parse_tuple(seeds[ident]['seed'], 15) == population[int(seeds[ident]['rank'])],
                'wrong saved construction')
        require(state_key(parse_tuple(seeds[ident]['seed'], 15)) == key, 'wrong seed alias')
    require(len(known_states) == len(seeds) == int(coverage['unique_seed_states']), 'seed count')
    for name in ['identity', 'intransitive', 'qualifying']:
        require(masses[name] == int(coverage[name]), 'mass mismatch '+name)
    from sympy.combinatorics import Permutation, PermutationGroup
    quotient_order = int(PermutationGroup([Permutation(list(BETA)), Permutation(list(GAMMA))]).order())
    require(quotient_order == 60, 'quotient is not A5')
    folders, checked_states = {}, 0
    for folder in sorted([*root.glob('job/component_*'), *root.glob('job/partial_*')]):
        check = audit_component(folder)
        seed = parse_tuple(check['summary']['seed'], 15)
        geometry = analyze(seed)
        qorder = int(PermutationGroup([Permutation(list(q)) for q in geometry['quotient_permutations']]).order())
        require(qorder == 60, 'saved quotient group')
        order = int(PermutationGroup([Permutation(list(p)) for p in seed]).order())
        require(order == int(check['summary']['group_order']), 'group order mismatch')
        for row in rows(folder/'STATES.tsv'):
            analyze(parse_tuple(row['tuple'], 15))
            checked_states += 1
        check['component_id'] = key_digest(check['state_keys']) if check['independent_six_move_closure'] else None
        check['independent_group_order'] = order
        check['geometry'] = geometry
        folders[folder.name] = check
    for seed, res in zip(seeds, resolutions):
        require(seed['candidate'] == res['candidate'], 'resolution ordering')
        if res['status'] == 'alias':
            check = folders[res['folder']]
            require(check['independent_six_move_closure'] and state_key(parse_tuple(seed['seed'], 15)) in set(check['state_keys']),
                    'false component alias')
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
                  quotient_group_order=quotient_order, saved_states_geometry_checked=checked_states,
                  independent_six_move_closed_components=len(catalog),
                  additional_closed_components=sum(not c['known_component'] for c in catalog),
                  partial_state_counts=len(partial), deferred_seed_states=sum(r['status']=='deferred' for r in resolutions),
                  independent_audit_seconds=time.monotonic()-started,
                  quotient='simultaneous S15 conjugation',
                  scope='All1296 fixed-quotient solved frames; selected closures, not all braid components')
    save(root/'CATALOG.json', dict(closed_components=catalog, partial_state_counts=partial))
    save(root/'CHECKS.json', result)
    return result


def worker(root, lock_fd):
    lock = os.fdopen(lock_fd, 'a')
    try:
        plan = json.loads((root/'PLAN.json').read_text())
        for name, digest in plan['source_hashes'].items():
            require(sha(root/'source'/name) == digest, 'changed frozen source '+name)
        require(sha(root/'smooth_noncyclic') == plan['binary_sha256'], 'changed binary')
        require(sha(root/'KNOWN_COMPONENTS.json') == plan['known_sha256'], 'changed known catalog')
        save(root/'STATE.json', dict(status='running', pid=os.getpid(), started=dt.datetime.now().isoformat()))
        cmd = [str(root/'smooth_noncyclic'), '--out', str(root/'job'), '--cap', str(plan['cap']),
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
    p.add_argument('--closures', type=int, default=4)
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
        print(json.dumps(dict(degree=15, blocks=5, size=3, solved_choices=1296,
              direct_same_family_choices=46656, quotient_group='A5', cpu_workers=1, gpu=False,
              closures=args.closures, closure_seconds=args.closure_seconds, cap=args.cap), indent=2))
        return
    import sympy
    dest = Path('data')
    require(Path('data').is_dir() and shutil.disk_usage(dest).free >= 30*2**30,
            'mounted T7 with30GiB free required')
    lock = (PROJECT/'.run-control/.campaign.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = dest/('smooth_noncyclic_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir()
    print('OUTPUT='+str(root), flush=True)
    source, hashes = root/'source', {}
    paths = [PROJECT/'Cargo.toml', PROJECT/'Cargo.lock', *sorted((PROJECT/'src').rglob('*.rs')),
             *sorted((PROJECT/'scripts').rglob('*.py')),
             *sorted((PROJECT/'evidence/smooth_noncyclic_20261012').glob('*.md')),
             PROJECT/'evidence/smooth_geometry_20261012/PROOF.md',
             PROJECT/'evidence/smooth_geometry_20261012/SOURCES.json',
             PROJECT/'evidence/smooth_special_20261012/CATALOG.json',
             *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')),
             PROJECT/'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv']
    for path in paths:
        rel = path.relative_to(PROJECT)
        target = source/rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(rel)] = sha(target)
    env = dict(os.environ, CARGO_BUILD_JOBS='1', CARGO_TARGET_DIR=str(PROJECT/'target/smooth-noncyclic'),
               RUSTFLAGS='-C target-cpu=native')
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo', 'test', '--release', '--bin', 'smooth_noncyclic', 'noncyclic_tests', '--', '--test-threads=1'],
                       cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cargo', 'build', '--release', '--bin', 'smooth_noncyclic'],
                       cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, 'scripts/smooth/test_smooth_noncyclic.py'], cwd=source, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = root/'smooth_noncyclic'
    shutil.copy2(PROJECT/'target/smooth-noncyclic/release/smooth_noncyclic', binary)
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
    save(root/'KNOWN_COMPONENTS.json', known)
    save(root/'PLAN.json', dict(source_hashes=hashes, binary_sha256=sha(binary),
         known_sha256=sha(root/'KNOWN_COMPONENTS.json'), sympy_version=sympy.__version__,
         closures=args.closures, closure_seconds=args.closure_seconds, cap=args.cap,
         cpu_workers=1, gpu=False, quotient='simultaneous S15 conjugation',
         population='fixed beta=(12345),gamma=(12453),six^four free S3 voltages',
         objective='new noncyclic-quotient positive-genus smooth examples; NOT weak failures'))
    with (root/'WORKER.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(source/'scripts/smooth/run_smooth_noncyclic.py'),
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
