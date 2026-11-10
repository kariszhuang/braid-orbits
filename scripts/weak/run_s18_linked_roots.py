#!/usr/bin/env python3
"""Finite linked-root census; shares the discovery lock, preserves raw evidence.

Run explicitly with --run. Enumeration coverage is not closure resolution.
Only an exact frame population, never all S18/A18 tuples, can be completed.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from analyze_fixed_braid_constraints import conjugate, point_orbits, product, power
from odd_order_spin import cycle_lengths, hurwitz, inverse, parse_seed

PROJECT = Path(__file__).resolve().parents[2]
DEST = Path('data')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')


def coords(h, k):
    return [[power(h, i)[power(k, j)[min(b)]] for i in range(3) for j in range(3)]
            for b in point_orbits([h, k])]


def frame_group(h, k):
    blocks = coords(h, k)
    for swap in range(2):
        for a in range(9):
            for b in range(9):
                shifts = [a, b]
                g = list(range(18))
                for block in range(2):
                    i0, j0 = divmod(shifts[block], 3)
                    for i in range(3):
                        for j in range(3):
                            g[blocks[block][3*i+j]] = blocks[block ^ swap][3*((i+i0)%3)+(j+j0)%3]
                yield tuple(g)


def state_key(t):
    """Independent equivariant rooted traversal; full transitive tuples only."""
    n = len(t[0])
    answers = []
    for root in range(n):
        order = [root]
        labels = {root: 0}
        for x in order:
            for p in t[:3]:
                y = p[x]
                if y not in labels:
                    labels[y] = len(order)
                    order.append(y)
        if len(order) != n:
            raise ValueError('canonicalization requires a transitive tuple')
        answers.append(tuple(tuple(labels[p[x]] for x in order) for p in t[:3]))
    return min(answers)


def certified_states(path):
    text = path.read_text().splitlines()
    d = int(next(s[2:] for s in text if s.startswith('d=')))
    if not 4<=d<=25:
        raise ValueError('only historical bit-packed degree 4..25 certificates supported')
    bits = 4 if d<=16 else 5
    for line in text:
        if not line.startswith('STATE\t'):
            continue
        perms = [tuple((int(z,16)>>(bits*i))&((1<<bits)-1) for i in range(d))
                 for z in line.split('\t')[1:]]
        if len(perms)!=3 or any(sorted(p)!=list(range(d)) for p in perms):
            raise ValueError('invalid certificate state encoding')
        yield tuple(perms+[inverse(product(*perms))])


def audit_output(root, audit):
    plan = json.loads((root/'PLAN.json').read_text())
    for path, expected in plan['source_hashes'].items():
        if sha(root/'source'/path)!=expected:
            raise ValueError(f'frozen source hash mismatch: {path}')
    if sha(root/'FRAME.tsv')!=plan['frame_sha256']:
        raise ValueError('frame input hash mismatch')
    execution = json.loads((root/'EXECUTION.json').read_text())
    if sha(root/'s18_linked_roots')!=execution['binary_sha256']:
        raise ValueError('frozen binary hash mismatch')
    rows = list(csv.DictReader((root/'FRAME_CLASSES.tsv').open(), delimiter='\t'))
    mass = sum(int(r['labeled_mass']) for r in rows)
    if mass != 856**2 or len({r['class'] for r in rows}) != len(rows):
        raise ValueError('frame class mass/count mismatch')
    seed_case = next(c for c in audit['cases'] if c['certificate']=='CERTIFIED_S18_XXXX_genus1_weak.txt')
    frame = {n: tuple(x-1 for x in p) for n,p in seed_case['seed']['ends']['maps'].items()}
    h, k = frame['h'], frame['k']
    group = list(frame_group(h,k))
    if len(set(group)) != 162 or any(conjugate(h,g)!=h or conjugate(k,g)!=k for g in group):
        raise ValueError('independent frame centralizer failed')
    def key(u,v):
        return min((conjugate(u,g),conjugate(v,g)) for g in group)
    # Independently check every reconstructed pair, literal algebra and defect.
    classes = {}
    for row in rows:
        t = parse_seed(row['seed'],18)
        u, v = product(t[0],inverse(h)), product(t[2],inverse(k))
        if power(u,2)!=k or power(v,2)!=h or product(*t)!=tuple(range(18)):
            raise ValueError('literal root reconstruction mismatch')
        genus = 1-18+sum(18-len(cycle_lengths(p)) for p in t)//2
        defect = len(cycle_lengths(product(*t[:2])))-len(point_orbits(t[:2]))-len(point_orbits(t[2:]))+1
        if genus!=int(row['genus']) or defect!=int(row['defect']):
            raise ValueError('genus/defect mismatch')
        # Rust's sorted root indices choose the minimum image-pair in each orbit.
        if (u,v) != key(u,v) or (u,v) in classes:
            raise ValueError('nonminimal or duplicate frame representative')
        stabilizer = sum(conjugate(u,g)==u and conjugate(v,g)==v for g in group)
        if 162//stabilizer != int(row['labeled_mass']):
            raise ValueError('orbit-stabilizer mass mismatch')
        classes[(u,v)] = row
    controls = []
    target_blocks = coords(h,k)
    for case in audit['cases']:
        if case['degree'] != 18:
            continue
        checks = [('seed',case['seed'])]+[(str(i),s) for i,s in enumerate(case['state_checks'])]
        for label, state in checks:
            ends = state.get('ends')
            if not ends or not ends['elementary_abelian_three_frame']:
                continue
            maps = {n:tuple(x-1 for x in p) for n,p in ends['maps'].items()}
            if maps['p'] != product(power(maps['h'],2),maps['k']):
                continue
            source_blocks = coords(maps['h'],maps['k'])
            g = list(range(18))
            for src,dst in zip(source_blocks,target_blocks):
                for x,y in zip(src,dst):
                    g[x]=y
            control = classes[key(conjugate(maps['u'],tuple(g)),conjugate(maps['v'],tuple(g)))]
            if control['status']!='closure_S':
                raise ValueError('known weak control pruned incorrectly')
            controls.append({'certificate':case['certificate'],'state':label,'frame_class':int(control['class'])})
    for name in ['CERTIFIED_S18_XXXX_genus1_weak.txt','CERTIFIED_S18_XXYY_genus1_weak.txt']:
        if not any(c['certificate']==name for c in controls):
            raise ValueError(f'missing positive control {name}')
    status_counts = Counter(r['status'] for r in rows)
    masses = Counter()
    for r in rows:
        masses[r['status']]+=int(r['labeled_mass'])
    resolved = True
    closure_checks = {}
    known = {}
    for path in (root/'source/evidence/certificates').glob('CERTIFIED_S18*weak.txt'):
        for t in certified_states(path):
            known[state_key(t)] = path.name
    for ambient in ['A','S']:
        expected = status_counts[f'closure_{ambient}']
        if not expected:
            continue
        directory = root/f'closure_{ambient}'
        closed = dict(known)
        for certificate in directory.glob('CERTIFIED_*.txt'):
            for t in certified_states(certificate):
                closed[state_key(t)]='new:'+certificate.name
        summary_path = directory/'SUMMARY.txt'
        if not summary_path.exists():
            resolved = False
            continue
        summary = dict(line.split(' = ',1) for line in summary_path.read_text().splitlines() if ' = ' in line)
        if int(summary['replay inputs'])!=expected:
            raise ValueError('unexpected S_d replay deduplication: investigate frame uniqueness')
        witnesses = {}
        lengths = Counter()
        for record in csv.DictReader((directory/'POSITIVE_WITNESSES.tsv').open(),delimiter='\t'):
            if record['degree']!='18':
                raise ValueError('wrong witness degree')
            seed = parse_seed(record['seed'],18)
            moves = [int(s) for s in record['moves'].split(',')] if record['moves'] else []
            t = seed
            for mv in moves:
                if mv not in range(6):
                    raise ValueError('invalid witness move')
                t = hurwitz(t,mv)
            if len(point_orbits(t[:2]))!=1:
                raise ValueError('false positive braid witness')
            witnesses[state_key(seed)] = moves
            lengths[len(moves)]+=1
        outcomes = Counter()
        for line in (root/f'REPLAY_{ambient}.txt').read_text().splitlines():
            key = state_key(parse_seed(line,18))
            if key in witnesses:
                outcomes['literal_positive_witness']+=1
            elif key in closed:
                outcomes[closed[key]]+=1
            else:
                outcomes['unresolved']+=1
        unresolved_promoted = 0
        for line in (directory/'survivors.txt').read_text().splitlines():
            key = state_key(parse_seed(line.split('seed=',1)[1],18))
            unresolved_promoted += key not in witnesses and key not in closed
        complete = int(summary['attempted seeds'])==expected and int(summary['hard unsolved candidates'])==0
        if (directory/'hard_candidates.txt').stat().st_size or unresolved_promoted or outcomes['unresolved']:
            complete = False
        closure_checks[ambient]={'expected':expected,'complete':complete,'summary':summary,
                                'sha256':sha(summary_path),'independent_outcomes':dict(outcomes),
                                'independent_witness_length_counts':dict(lengths),
                                'unresolved_promoted_records':unresolved_promoted,
                                'note':'survivors.txt durably logs promoted starts, including subsequently resolved ones'}
        resolved &= complete
    return {'population':'fixed C3^2 two-regular-block frame with p=h^2*k, u^2=k, v^2=h',
            'labeled_pairs':mass,'frame_classes':len(rows),'classification_complete':True,
            'fully_resolved':resolved,'classes_by_status':dict(status_counts),
            'labeled_pairs_by_status':dict(masses),'known_component_controls':controls,
            'closure_checks':closure_checks,'table_sha256':sha(root/'FRAME_CLASSES.tsv'),
            'warning':'Not a whole-degree census; proper groups/identity triples remain separate strata.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',action='store_true')
    parser.add_argument('--threads',type=int,default=7)
    parser.add_argument('--seconds',type=int,default=3600)
    parser.add_argument('--candidate-seconds',type=int,default=30)
    parser.add_argument('--audit-root',type=Path,help='independent completion audit of an existing frozen run; no new search')
    parser.add_argument('--checks-out',type=Path,help='fresh audit JSON path (never overwritten)')
    args = parser.parse_args()
    if args.audit_root:
        root=args.audit_root.resolve()
        audit=json.loads((root/'source/evidence/fixed_braid_math_20261011/AUDIT.json').read_text())
        checks=audit_output(root,audit)
        checks['completion_verifier_sha256']=sha(Path(__file__).resolve())
        checks['helper_hashes']={name:sha(_scripts_path.find(name)) for name in ['analyze_fixed_braid_constraints.py','odd_order_spin.py']}
        destination=args.checks_out or root/('CHECKS_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S')+'.json')
        save(destination,checks)
        print(json.dumps(checks,indent=2))
        return
    if not args.run:
        print('Plan: 856^2 exact linked-root pairs, quotient frame centralizer, generation, exact closure; use --run.')
        return
    if not 1<=args.threads<=7 or args.seconds<1 or args.candidate_seconds<1:
        parser.error('require 1..7 workers and positive time limits')
    if not Path('data').is_dir() or not DEST.is_dir():
        raise RuntimeError('T7 must be mounted')
    if shutil.disk_usage(DEST).free < 30*2**30:
        raise RuntimeError('require 30 GiB free on T7')
    lock = (PROJECT/'.run-control/.campaign.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    root = DEST/('linked_roots_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir()
    print(f'OUTPUT={root}',flush=True)
    source = root/'source'; source.mkdir()
    paths = [PROJECT/'Cargo.toml',PROJECT/'Cargo.lock',Path(__file__).resolve(),
             *sorted((PROJECT/'src').rglob('*.rs')),
             PROJECT/'evidence/fixed_braid_math_20261011/AUDIT.json',
             PROJECT/'scripts/weak/analyze_fixed_braid_constraints.py',PROJECT/'scripts/weak/odd_order_spin.py',
             *sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*_weak.txt'))]
    hashes = {}
    for path in paths:
        relative = path.relative_to(PROJECT)
        target = source/relative; target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target); hashes[str(relative)]=sha(target)
    audit_path = source/'evidence/fixed_braid_math_20261011/AUDIT.json'
    audit = json.loads(audit_path.read_text())
    maps = next(c for c in audit['cases'] if c['certificate']=='CERTIFIED_S18_XXXX_genus1_weak.txt')['seed']['ends']['maps']
    with (root/'FRAME.tsv').open('x') as f:
        for name in ['h','k','p']:
            f.write(name+'\t'+','.join(map(str,maps[name]))+'\n')
    save(root/'PLAN.json',{'population':'856^2 root pairs in fixed S18 C3^2 frame',
         'quotient':'centralizer of h,k,p (order 162); full tuples certified modulo S18',
         'source_hashes':hashes,'frame_sha256':sha(root/'FRAME.tsv'),
         'cpu_workers':args.threads,'gpu_lanes':0,'classification':'complete finite enumeration',
         'closure_limits':{'seconds':args.seconds,'candidate_seconds':args.candidate_seconds,'states':200000},
         'completion':'all masses reconciled, controls recovered, replay inputs attempted, promoted records reconciled, zero hard leads, certificates verified',
         'previous_campaign':'data/campaign_20261010_212021_242369'})
    env = dict(os.environ,CARGO_TARGET_DIR=str(PROJECT/'target/linked-roots'),CARGO_BUILD_JOBS='1',RUSTFLAGS='-C target-cpu=native')
    with (root/'BUILD.log').open('x') as log:
        subprocess.run(['cargo','test','--release','--bin','s18_linked_roots','linked_roots_tests','--','--test-threads=1'],cwd=PROJECT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run(['cargo','build','--release','--bin','s18_linked_roots'],cwd=PROJECT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    binary = root/'s18_linked_roots'
    shutil.copy2(PROJECT/'target/linked-roots/release/s18_linked_roots',binary)
    command = [str(binary),'--frame',str(root/'FRAME.tsv'),'--out',str(root),
               '--threads',str(args.threads),'--seconds',str(args.seconds),
               '--candidate-seconds',str(args.candidate_seconds),
               '--known-certificates',str(source/'evidence/certificates')]
    save(root/'EXECUTION.json',{'binary_sha256':sha(binary),'command':command})
    with (root/'RUN.log').open('x') as log:
        subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,check=True)
    checks = audit_output(root,audit)
    save(root/'CHECKS.json',checks)
    print(json.dumps(checks,indent=2),flush=True)
    lock.close()


if __name__=='__main__':
    main()
