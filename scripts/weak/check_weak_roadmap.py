#!/usr/bin/env python3
"""Finite, locked validation of the seven-direction weak-search roadmap.

This is a sequential evidence harness, not a second discovery scheduler.
Requires a stopped campaign, a tested release binary, and mounted T7.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from run_campaign import CONTROL, RESULTS, atomic_json, require_results_volume


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--frontier',action='store_true',help='validate policy-five finite S3, pruning, positive-genus gluing and ninth certificate only')
    parser.add_argument('--patterns',action='store_true',help='bounded policy-six marked-pair/repair pilots and paired construction benchmark')
    parser.add_argument('--fibers',action='store_true',help='bounded policy-seven cut/join parity-bridge pilots and identical-input work-order benchmark')
    parser.add_argument('--refinement',action='store_true',help='finite S18-to-A18 one-cut shards, post-cut gluing pilots and nine controls')
    args = parser.parse_args()
    if sum((args.patterns,args.frontier,args.fibers,args.refinement))>1:parser.error('choose one validation mode')
    root = args.out.resolve()
    require_results_volume(root)
    if shutil.disk_usage(RESULTS).free < 30 * (1 << 30):
        raise RuntimeError('T7 reserve below 30 GiB')
    CONTROL.mkdir(exist_ok=True)
    with (CONTROL / '.campaign.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        root.mkdir(exist_ok=False)
        sources = ['src/main.rs', 'src/smart.rs', 'src/smart/weak_discovery.rs', 'src/smart/pair_fiber.rs', 'src/smart/braid_frames.rs', 'src/smart/positive_witness.rs',
                   'src/smart/defect_mission.rs', 'src/smart/tests.rs', *[str(p.relative_to(PROJECT)) for p in sorted((PROJECT/'src/smart/weak_discovery').glob('*.rs'))], 'Cargo.toml', 'Cargo.lock',
                   'scripts/weak/run_campaign.py', 'scripts/weak/test_campaign_runner.py',
                   'scripts/weak/analyze_weak_mechanisms.py', 'scripts/weak/check_weak_roadmap.py']
        (root / 'source').mkdir()
        hashes = {}
        for name in sources:
            path = PROJECT / name
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            hashes[name] = digest
            shutil.copy2(path, root / 'source' / (digest + '_' + path.name))
        binary = root / 'hurwitz_breaker'
        shutil.copy2(PROJECT / 'target/smart/release/hurwitz_breaker', binary)
        certs = sorted((PROJECT / 'evidence/certificates').glob('CERTIFIED_*weak.txt'))
        (root / 'known').mkdir()
        for cert in certs:
            shutil.copy2(cert, root / 'known' / cert.name)
            hashes[str(cert.relative_to(PROJECT))] = hashlib.sha256(cert.read_bytes()).hexdigest()
        plan = dict(started=dt.datetime.now().isoformat(), source_sha256=hashes,
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    resources='one CPU process at a time; no competing campaign', jobs=[],
                    scope='bounded technical pilots and selected exact envelope/certificate analyses; not a census')
        atomic_json(root / 'PLAN.json', plan)

        def run(name, command):
            print('CHECK', name, flush=True)
            record = dict(name=name, command=list(map(str, command)), started=dt.datetime.now().isoformat())
            plan['jobs'].append(record)
            atomic_json(root / 'PLAN.json', plan)
            started = time.monotonic()
            with (root / (name + '.log')).open('x') as log:
                process = subprocess.run(record['command'], cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT)
            record.update(returncode=process.returncode, wall_seconds=time.monotonic() - started)
            atomic_json(root / 'PLAN.json', plan)
            if process.returncode:
                raise RuntimeError(f'{name} failed; original log preserved')
            return root / name

        known = root / 'known'
        common = [binary, 'smart-hunt', '--threads', 1, '--weak-only', '--all-pairs-intransitive',
                  '--known-certificates', known, '--seconds', 3, '--candidate-seconds', 1,
                  '--cap', 5000, '--steps', 64, '--seed', 20261002]
        if args.refinement:
            from run_campaign import weak_symmetry_rows
            plan['scope']='complete specified one-cut S18-to-A18 constructions; bounded broadened gluing pilots; independent nine-control verification; not an A18 census'
            common[common.index('--seed')+1]=20261004
            common[common.index('--seconds')+1]=30
            common[common.index('--candidate-seconds')+1]=5
            common[common.index('--cap')+1]=200000
            common+=['--screen-steps',64,'--skip-known-starts']
            rows=weak_symmetry_rows(0,known,include_fibers=True,include_refinement=True)
            for row in rows:
                if row['family']!='pair-fiber-exact':continue
                name=row['name']
                run(name,common+['--family',row['family'],'--group','A','--d',18,
                    '--certificate',row['source'],'--fiber-state-start',row['fiber_state_start'],
                    '--fiber-state-count',row['fiber_state_count'],'--out',root/name])
            for group,d,genus in [('A',9,1),('S',9,0),('A',10,1),('S',10,1),('A',11,1),('S',11,3),('A',18,1),('S',14,3)]:
                for expanded in (False,True):
                    name=f'base_{"postcut" if expanded else "original"}_{group}{d}_g{genus}'
                    command=common.copy();command[command.index('--seconds')+1]=2
                    run(name,command+['--family','cut-join-glue','--group',group,'--d',d,'--glue-genus',genus,
                        '--pair-genus-max',3,'--expand-pair-actions',*(['--post-cut-base-frames'] if expanded else []),'--out',root/name])
            for cert in certs:run('verify_'+cert.stem,[binary,'verify','--certificate',known/cert.name])
            cargo=['env','CARGO_BUILD_JOBS=1','CARGO_TARGET_DIR=target/smart','RUSTFLAGS=-C target-cpu=native','cargo','test','--release']
            run('cut_listing_benchmark',cargo+['benchmark_finite_cut_listing','--','--ignored','--nocapture','--test-threads=1'])
            for candidate in sorted(root.glob('*/hard_candidates.txt')):
                if not candidate.read_text().strip():continue
                summary=dict(l.split(' = ',1) for l in (candidate.parent/'SUMMARY.txt').read_text().splitlines() if ' = ' in l)
                group,degree=summary['group'].split('_');name='replay_'+candidate.parent.name
                run(name,[binary,'smart-hunt','--group',group,'--d',degree,'--threads',1,'--weak-only',
                    '--known-certificates',known,'--replay',candidate,'--replay-only','--seconds',30,
                    '--candidate-seconds',15,'--steps',4096,'--cap',1000000,'--seed',20261004,'--out',root/name])
            for path in sorted(root.glob('*/CERTIFIED_*.txt')):
                if path.parent==known:continue
                run('verify_new_'+hashlib.sha256(path.read_bytes()).hexdigest()[:16],[binary,'verify','--certificate',path])
            checks=[]
            for job in plan['jobs']:
                path=root/job['name']/'SUMMARY.txt'
                values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l) if path.exists() else {}
                checks.append(dict(**job,summary=values))
            atomic_json(root/'CHECKS.json',dict(scope=plan['scope'],binary_sha256=plan['binary_sha256'],checks=checks))
            plan['completed']=dt.datetime.now().isoformat();atomic_json(root/'PLAN.json',plan)
            print('REFINEMENT CHECKS COMPLETE',root,flush=True)
            return
        if args.fibers:
            plan['scope']='bounded fixed-product fiber/parity-bridge/gluing technical pilots; nine exact controls; identical-input paired work-order benchmark; not a census'
            common[common.index('--seed')+1]=20261004
            common+=['--screen-steps',64,'--skip-known-starts']
            source_cases=[]
            for stem in ('S18_genus1','S18_XXYY_genus1','S18_XXXX_genus1'):
                source_cases.extend(('A',18,stem,r) for r in (0,2))
            source_cases.extend(('S',18,'S18_genus1',r) for r in (0,2))
            for group in ('A','S'):
                source_cases.extend((group,25,'A25_XXYY_genus4',r) for r in (0,2))
                source_cases.append((group,25,'A25',0))
            source_cases.append(('S',25,'A25_order12_genus12',0))
            for group,d,stem,radius in source_cases:
                name=f'fiber_{group}{d}_{stem}_r{radius}'
                run(name,common+['--family','pair-fiber','--group',group,'--d',d,'--fiber-radius',radius,
                    '--certificate',known/f'CERTIFIED_{stem}_weak.txt','--out',root/name])
            for group,d,genus in [('S',9,0),('S',9,1),('A',9,1),('S',10,1),('A',10,0),('A',11,1),
                ('S',11,3),('S',14,1),('S',14,3),('A',18,1),('A',25,4)]:
                name=f'fiber_glue_{group}{d}_g{genus}'
                run(name,common+['--family','cut-join-glue','--group',group,'--d',d,'--glue-genus',genus,
                    '--pair-genus-max',3,'--expand-pair-actions','--fiber-radius',0,'--out',root/name])
            for cert in certs:run('verify_'+cert.stem,[binary,'verify','--certificate',known/cert.name])
            # Job-deadline tails are leads, not failures. Replay saved caps
            # sequentially with a larger budget and keep any remaining files.
            for candidate in sorted(root.glob('fiber_*/hard_candidates.txt')):
                if not candidate.read_text().strip():continue
                summary=dict(l.split(' = ',1) for l in (candidate.parent/'SUMMARY.txt').read_text().splitlines() if ' = ' in l)
                group,degree=summary['group'].split('_')
                name='replay_'+candidate.parent.name
                run(name,[binary,'smart-hunt','--group',group,'--d',degree,'--threads',1,'--weak-only',
                    '--known-certificates',known,'--replay',candidate,'--replay-only','--seconds',30,
                    '--candidate-seconds',15,'--steps',4096,'--cap',1000000,'--seed',20261004,'--out',root/name])
            cargo=['env','CARGO_BUILD_JOBS=1','CARGO_TARGET_DIR=target/smart','RUSTFLAGS=-C target-cpu=native','cargo','test','--release']
            run('order_benchmark',cargo+['benchmark_witness_before_generation','--','--ignored','--nocapture','--test-threads=1'])
            for path in sorted(root.glob('*/CERTIFIED_*.txt')):
                if path.parent==known:continue
                run('verify_new_'+hashlib.sha256(path.read_bytes()).hexdigest()[:16],[binary,'verify','--certificate',path])
            checks=[]
            for job in plan['jobs']:
                path=root/job['name']/'SUMMARY.txt'
                values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l) if path.exists() else {}
                checks.append(dict(**job,summary=values))
            atomic_json(root/'CHECKS.json',dict(scope=plan['scope'],binary_sha256=plan['binary_sha256'],checks=checks))
            plan['completed']=dt.datetime.now().isoformat();atomic_json(root/'PLAN.json',plan)
            print('FIBER CHECKS COMPLETE',root,flush=True)
            return
        if args.patterns:
            plan['scope']='bounded marked pair/dihedral/resistant-glue technical pilots, exact controls and identical-start paired benchmark; no census or yield forecast'
            common[common.index('--seconds')+1]=2
            common[common.index('--seed')+1]=20261003
            for group,d,model,power in [('S',9,'A4',1),('S',10,'S4',2),('A',11,'A5',1),('S',11,'D5',1),
                ('S',12,'D7',2),('S',20,'D9',4),('S',24,'D11',5),('S',28,'D13',6)]:
                name=f'marked_{group}{d}_{model}_r{power}'
                extra=['--novel-power'] if power>3 else []
                run(name,common+['--family','marked-symmetry','--group',group,'--d',d,'--pair-model',model,
                    '--pure-power',power,'--marked-pairs','all' if model.startswith('D') else 'novel',*extra,'--out',root/name])
            for group,d,genus in [('S',10,1),('S',10,3),('A',11,5),('S',14,3),('S',14,5),('S',18,1),('A',25,4),('A',25,12)]:
                name=f'resistant_{group}{d}_g{genus}'
                run(name,common+['--family','resistant-glue','--group',group,'--d',d,'--glue-genus',genus,
                    '--pair-genus-max',3,'--expand-pair-actions','--glue-repair-steps',16,'--out',root/name])
            for cert in certs:run('verify_'+cert.stem,[binary,'verify','--certificate',known/cert.name])
            cargo=['env','CARGO_BUILD_JOBS=2','CARGO_TARGET_DIR=target/smart','RUSTFLAGS=-C target-cpu=native','cargo','test','--release']
            run('marked_model_counts',cargo+['marked_groups_and_new_dihedral_actions_are_exact','--','--nocapture','--test-threads=1'])
            run('boundary_benchmark',cargo+['benchmark_boundary_resistance','--','--ignored','--nocapture','--test-threads=1'])
            for path in sorted(root.glob('*/CERTIFIED_*.txt')):
                run('verify_new_'+hashlib.sha256(path.read_bytes()).hexdigest()[:16],[binary,'verify','--certificate',path])
            checks=[]
            for job in plan['jobs']:
                path=root/job['name']/'SUMMARY.txt'
                values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l) if path.exists() else {}
                checks.append(dict(**job,summary=values))
            atomic_json(root/'CHECKS.json',dict(scope=plan['scope'],binary_sha256=plan['binary_sha256'],checks=checks))
            plan['completed']=dt.datetime.now().isoformat();atomic_json(root/'PLAN.json',plan)
            print('PATTERN CHECKS COMPLETE',root,flush=True)
            return
        if args.frontier:
            from run_campaign import finite_passport_exclusions
            exclusions=root/'FINITE_PASSPORT_EXCLUSIONS.tsv'
            exclusions.write_text('# group\tdegree\tpassport\tverified provenance\n'+''.join(
                f'{g}\t{d}\t{"/".join(p)}\t{source}\n' for (g,d,p),source in sorted(finite_passport_exclusions().items())))
            plan['scope']='first-frame finite S3 technical pilots in A/S9..11 (not full construction censuses); bounded positive-genus glue/source checks'
            plan['excluded_passports_sha256']=hashlib.sha256(exclusions.read_bytes()).hexdigest()
            for d in (9,10,11):
                for group in ('A','S'):
                    name=f's3_exact_{group}{d}'
                    run(name,common+['--family','s3-exact','--group',group,'--d',d,'--frame-cap',1,'--excluded-passports',exclusions,'--out',root/name])
            for group,d,genus in [('S',10,2),('A',11,4),('S',14,2),('A',25,4)]:
                name=f'pair_glue_{group}{d}_g{genus}'
                run(name,common+['--family','pair-glue','--group',group,'--d',d,'--glue-genus',genus,'--pair-genus-max',2,'--out',root/name])
            cert=known/'CERTIFIED_A25_XXYY_genus4_weak.txt'
            run('verify_ninth',[binary,'verify','--certificate',cert])
            run('profile_ninth',[binary,'certificate-profile','--certificate',cert,'--out',root/'NINTH_PROFILE.tsv'])
            for power in (1,2,3):
                name=f'ninth_source_r{power}'
                run(name,common+['--family','pure-centralizer','--group','A','--d',25,'--certificate',cert,'--pure-power',power,'--out',root/name])
            for path in sorted(root.glob('*/CERTIFIED_*.txt')):
                run('verify_'+hashlib.sha256(path.read_bytes()).hexdigest()[:16],[binary,'verify','--certificate',path])
            checks=[]
            for job in plan['jobs']:
                path=root/job['name']/'SUMMARY.txt'
                values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l) if path.exists() else {}
                checks.append(dict(**job,summary=values))
            atomic_json(root/'CHECKS.json',dict(scope=plan['scope'],binary_sha256=plan['binary_sha256'],checks=checks))
            plan['completed']=dt.datetime.now().isoformat();atomic_json(root/'PLAN.json',plan)
            print('FRONTIER CHECKS COMPLETE',root,flush=True)
            return
        for group, d, model in [('A', 9, 'A4'), ('S', 10, 'S4'), ('A', 11, 'A5'), ('A', 14, 'PSL27')]:
            name = f'audit_{group}{d}_{model}'
            run(name, [binary, 'smart-hunt', '--audit-construction', '--family', 'pure-symmetry',
                       '--group', group, '--d', d, '--pair-model', model, '--frame-cap', 64,
                       '--twist-cap', 256, '--out', root / name])
        for group, d, model, power in [('S', 9, 'A4', 2), ('A', 10, 'S4', 3), ('A', 11, 'A5', 2),
                                      ('A', 9, 'PSL27', 1), ('A', 14, 'PSL27', 2), ('S', 15, 'PSL27', 3)]:
            name = f'pure_{group}{d}_{model}_r{power}'
            run(name, common + ['--family', 'pure-symmetry', '--group', group, '--d', d,
                                '--pair-model', model, '--pure-power', power, '--out', root / name])
        for group, d, genus in [('A', 9, 1), ('S', 10, 0), ('A', 11, 1), ('S', 14, 1), ('A', 17, 1)]:
            name = f'glue_{group}{d}_g{genus}'
            run(name, common + ['--family', 'dessin-glue', '--group', group, '--d', d,
                                '--glue-genus', genus, '--out', root / name])
        for ambient, filename in [('S', 'CERTIFIED_S18_XXXX_genus1_weak.txt'), ('A', 'CERTIFIED_A25_mixed_genus0_weak.txt')]:
            d = 18 if ambient == 'S' else 25
            for power in (2, 3):
                name = f'known_{ambient}{d}_r{power}'
                run(name, common + ['--family', 'pure-centralizer', '--group', ambient, '--d', d,
                                    '--certificate', known / filename, '--pure-power', power, '--out', root / name])
        # Paired technical throughput controls; alternating order reduces
        # session drift. These are sample-occurrence rates, not discovery odds.
        for repetition in range(3):
            for enabled in ((True, False) if repetition % 2 == 0 else (False, True)):
                name = f'perf_S12_rep{repetition}_envelope{int(enabled)}'
                run(name, common + ['--family', 'pure-symmetry', '--group', 'S', '--d', 12,
                                    '--pair-model', 'S4', '--out', root / name] + ([] if enabled else ['--no-envelope']))
        audit = run('mechanisms', [binary, 'mechanism-audit', '--certificates', ','.join(map(str, sorted(known.glob('CERTIFIED_*weak.txt')))), '--out', root / 'mechanisms'])
        run('exact_braid_images', ['uv', 'run', '--offline', '--isolated', '--no-project', '--python', '3.13',
                                  '--with', 'sympy==1.14.0', 'python', PROJECT / 'scripts/weak/analyze_weak_mechanisms.py',
                                  '--audit', audit, '--out', root / 'BRAID_IMAGES.json'])
        for d in (9, 10, 11):
            for group in ('A', 'S'):
                name = f'zero_defect_{group}{d}'
                run(name, [binary, 'defect-mission', '--mode', 'walk', '--group', group, '--d', d,
                           '--threads', 1, '--draws', 8, '--walk-steps', 1024, '--rng', 20261002 + 2 * d + (group == 'S'),
                           '--seconds', 2, '--cap', 5000, '--family', 'low-index', '--out', root / name])
        for cert in sorted(root.glob('*/CERTIFIED_*.txt')):
            run('verify_' + cert.parent.name + '_' + cert.stem, [binary, 'verify', '--certificate', cert])
        summaries = {}
        for path in sorted(root.glob('*/SUMMARY.txt')):
            summaries[str(path.relative_to(root))] = dict(line.split(' = ', 1) for line in path.read_text().splitlines() if ' = ' in line)
        zero = []
        for path in sorted(root.glob('zero_defect_*/draw_*/SUMMARY.tsv')):
            with path.open() as stream:
                zero.extend(csv.DictReader(stream, delimiter='\t'))
        checks = dict(binary_sha256=plan['binary_sha256'], source_sha256=hashes,
                      jobs=plan['jobs'], summaries=summaries, zero_defect_draws=zero,
                      output_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in sorted(root.rglob('*')) if p.is_file() and 'source' not in p.parts and p != binary},
                      limitations='Technical pilots are biased. Missing witnesses and capped closure are not counterexamples. Exact polynomials do not establish a structural classification.')
        atomic_json(root / 'CHECKS.json', checks)
        (root / 'DONE').write_text('All selected checks exited successfully. Not a degree census.\n')
        print('DONE', root, flush=True)


if __name__ == '__main__':
    main()
