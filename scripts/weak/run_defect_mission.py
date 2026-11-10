#!/usr/bin/env python3
"""Plan first; run only with --run. Frozen, locked, resumable defect questions.

No new census or parity walk is scheduled for already answered S4-S7 questions.
The S8 zero screen stops at a positive witness instead of rebuilding spectra.
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
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from run_campaign import CONTROL, RESULTS, atomic_json, require_results_volume, run_event
from compute_spectra import component_files


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1 << 20), b''):
            h.update(data)
    return h.hexdigest()


def s8_selection(limit=24):
    """Deterministic stratification, not random population sampling."""
    from run_campaign import seed_profile
    files, missing = component_files('S', 8, 8)
    if missing:
        raise RuntimeError('Completed S8 catalogs are unavailable')
    blocks = {}
    for file in files:
        shard = file.parent.name
        with file.open(newline='') as stream:
            for row in csv.DictReader(stream):
                _, genus = seed_profile(row['seed'], 8)
                n = int(row['quotient_states'])
                identity = any(part.strip() == '()' for part in row['seed'][1:-1].split(','))
                if genus <= 0 or identity or n > 512:
                    continue
                bucket = 0 if n <= 24 else 1 if n <= 128 else 2
                case = f'S8:{shard}:component-{row["component"]}'
                blocks.setdefault((genus, bucket), []).append((n, case))
    selected = []
    # Interleave strata; smallest and largest eligible graph from each stratum.
    for offset in (0, -1):
        for block in sorted(blocks):
            case = sorted(blocks[block])[offset][1]
            if case not in selected:
                selected.append(case)
            if len(selected) == limit:
                return files, selected
    return files, selected


def mission_plan(threads=5, draws=128, seconds=30, cap=250000):
    files, missing = component_files('A', 4, 8)
    if missing:
        raise RuntimeError(f'Missing alternating catalogs: {missing}')
    small = [p for p in files if any(parent.name in ('A4', 'A5', 'A6') for parent in p.parents)]
    larger = [p for p in files if p not in small]
    s8, ids = s8_selection()
    certs = sorted((PROJECT / 'evidence/certificates').glob('*weak.txt'))
    jobs = [dict(name='known-inner', mode='profile', certificates=list(map(str, certs)), state_rows=True, inner_spectra=True),
            dict(name='alternating-defect-profiles', mode='profile', components=list(map(str, small)), positive_nonidentity=True),
            dict(name='alternating-zero-screen', mode='screen', components=list(map(str, larger))),
            dict(name='S8-zero-screen', mode='screen', components=list(map(str, s8))),
            dict(name='S8-selected-profiles', mode='profile', components=list(map(str, s8)), selected=ids, positive_nonidentity=True)]
    for d in (9, 10, 11, 12, 13, 14, 18, 25):
        for group in ('A', 'S'):
            jobs.append(dict(name=f'zero-walk-{group}{d}', mode='walk', group=group, degree=d,
                             draws=draws, family='low-index', rng=20261002 + 2*d + (group == 'S')))
    jobs += [dict(name='S8-class-nulls', mode='nulls', components=list(map(str, s8)), selected=ids),
             dict(name='completed-S8-summary', mode='summary')]
    sources = set(files + s8 + certs)
    sources.update(PROJECT / p for p in ('Cargo.toml', 'Cargo.lock', 'src/main.rs', 'src/smart.rs',
        'src/smart/defect_mission.rs', 'src/smart/weak_discovery.rs', 'src/smart/tests.rs', 'src/hybrid.rs', 'scripts/weak/run_campaign.py',
        'scripts/spectra/compute_spectra.py', 'scripts/weak/run_defect_mission.py', 'scripts/weak/analyze_defect_mission.py'))
    sources.update(PROJECT / f'runs/spectra/S8_exact_shards/shard_{i}/SPECTRUM.tsv' for i in range(16))
    sources.update((PROJECT/'src').rglob('*.rs'))
    return dict(version=1, threads=threads, cap=cap, seconds=seconds, jobs=jobs,
                input_hashes={str(p): sha(p) for p in sorted(sources)},
                limitations='Selected catalogs and biased zero-start discovery. Caps are unresolved. No parity-only run; even-defect is already disproved. No new counterexample-search campaign is started.')


def rust_command(binary, job, output, plan):
    cmd = [str(binary), 'defect-mission', '--mode', job['mode'], '--out', str(output),
           '--threads', str(plan['threads']), '--cap', str(plan['cap']), '--seconds', str(plan['seconds'])]
    for key in ('components', 'certificates'):
        if key in job:
            cmd += [f'--{key}', ','.join(job[key])]
    if 'selected' in job:
        selection = output.parent / 'SELECTED.txt'
        selection.write_text('\n'.join(job['selected']) + '\n')
        cmd += ['--select-case-ids', str(selection)]
    for key in ('state_rows', 'inner_spectra', 'positive_nonidentity'):
        if job.get(key):
            cmd.append('--' + key.replace('_', '-'))
    if job['mode'] == 'walk':
        cmd += ['--d', str(job['degree']), '--group', job['group'], '--draws', str(job['draws']),
                '--family', job['family'], '--rng', str(job['rng']), '--walk-steps', '4096']
    return cmd


def verify_hashes(hashes):
    for file, expected in hashes.items():
        if sha(file) != expected:
            raise RuntimeError(f'Frozen input/output changed: {file}')


def output_hashes(output):
    return {str(p): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}


def preserve_incomplete(output):
    """Move only unfinished case/draw folders aside; never erase original files."""
    if not output.exists():
        return
    for folder in sorted(output.iterdir()):
        if folder.is_dir() and folder.name.startswith(('case_', 'draw_')) and not (folder / 'DONE').exists():
            archive = output.parent / 'interrupted' / dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            archive.mkdir(parents=True, exist_ok=False)
            folder.rename(archive / folder.name)


def completed_hashes(output):
    hashes = {}
    for folder in sorted(output.glob('*_*')):
        if folder.is_dir() and (folder / 'DONE').is_file():
            hashes.update(output_hashes(folder))
    return hashes


def run(root, plan, resume=False):
    require_results_volume(root)
    CONTROL.mkdir(exist_ok=True)
    with (CONTROL / '.campaign.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify_hashes(plan['input_hashes'])
        if shutil.disk_usage(RESULTS).free < 30 * (1 << 30):
            raise RuntimeError('Less than 30 GiB free on T7')
        if not resume:
            root.mkdir(parents=True, exist_ok=False)
            atomic_json(root / 'PLAN.json', plan)
        if 'binary_hash' not in plan:
            env = dict(os.environ, CARGO_TARGET_DIR='target/smart', RUSTFLAGS='-C target-cpu=native')
            subprocess.run(['cargo', 'test', '--release'], cwd=PROJECT, env=env, check=True)
            subprocess.run(['cargo', 'build', '--release'], cwd=PROJECT, env=env, check=True)
            binary = root / 'bin/hurwitz_breaker'
            binary.parent.mkdir(exist_ok=True)
            shutil.copy2(PROJECT / 'target/smart/release/hurwitz_breaker', binary)
            plan['binary_hash'] = sha(binary)
            atomic_json(root / 'PLAN.json', plan)
        binary = root / 'bin/hurwitz_breaker'
        verify_hashes({str(binary): plan['binary_hash']})
        state_file = root / 'STATE.json'
        state = json.loads(state_file.read_text()) if state_file.exists() else dict(completed={}, case_checkpoints={})
        child = None
        def stop(signum, frame):
            raise KeyboardInterrupt(f'signal {signum}')
        signal.signal(signal.SIGTERM, stop)
        try:
            for job in plan['jobs']:
                name = job['name']
                if name in state['completed']:
                    verify_hashes(state['completed'][name])
                    continue
                require_results_volume(root)
                if shutil.disk_usage(root).free < 30 * (1 << 30):
                    raise RuntimeError('T7 disk reserve reached; resume after freeing space')
                folder = root / 'jobs' / name
                folder.mkdir(parents=True, exist_ok=True)
                output = folder / 'output'
                verify_hashes(state['case_checkpoints'].get(name, {}))
                preserve_incomplete(output)
                if job['mode'] == 'summary':
                    if output.exists():
                        output.rename(folder / ('interrupted-summary-' + dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')))
                    command = [sys.executable, str(PROJECT / 'scripts/weak/analyze_defect_mission.py'), '--root', str(root), '--out', str(output)]
                elif job['mode'] == 'nulls':
                    profiles = root / 'jobs/S8-selected-profiles/output'
                    ids, codes = [], set()
                    for file in sorted(profiles.glob('case_*/SUMMARY.tsv')):
                        code_file = file.parent / 'GRAPH_CODE.txt'
                        if not (file.parent / 'DONE').exists() or not code_file.exists():
                            continue
                        code = code_file.read_text()
                        if not code or code in codes:
                            continue
                        codes.add(code)
                        with file.open() as stream:
                            ids.append(next(csv.DictReader(stream, delimiter='\t'))['case_id'])
                    if not ids:
                        raise RuntimeError('No closed selected graph classes for null comparison')
                    # Reuse exact class identification, all four nulls and the
                    # existing class-balanced analyzer; no duplicated null engine.
                    command = [str(binary), 'research', '--components', ','.join(job['components']), '--only', ','.join(ids),
                        '--out', str(output), '--max-lift', '1', '--standard-max', '1', '--null-max', '512',
                        '--null-reps', '128', '--full-max', '512', '--lanczos-steps', '384', '--trace-probes', '8']
                    atomic_json(folder / 'CLASS_SELECTION.json', dict(ids=ids, classes=len(codes)))
                    if output.exists():
                        preserved = folder / ('interrupted-nulls-' + dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
                        output.rename(preserved)
                else:
                    command = rust_command(binary, job, output, plan)
                atomic_json(folder / 'COMMAND.json', command)
                state.update(status='running', job=name, pid=os.getpid())
                atomic_json(state_file, state)
                run_event(root, 'job_started', job=name, command=command)
                with (folder / 'RUN.log').open('a') as log:
                    child = subprocess.Popen(command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT)
                    child.wait()
                    if child.returncode:
                        raise RuntimeError(f'{name} exited {child.returncode}; see {folder / "RUN.log"}')
                    child = None
                if job['mode'] == 'nulls':
                    subprocess.run([str(binary), 'analyze-nulls', '--inputs', str(output), '--out', str(output / 'analysis')], cwd=PROJECT, check=True)
                state['completed'][name] = output_hashes(output)
                atomic_json(state_file, state)
                run_event(root, 'job_verified', job=name, files=len(state['completed'][name]))
            state.update(status='complete', job=None)
            atomic_json(state_file, state)
        except BaseException as error:
            if child is not None and child.poll() is None:
                child.terminate()
                child.wait()
            current = state.get('job')
            if current:
                state['case_checkpoints'][current] = completed_hashes(root / 'jobs' / current / 'output')
            state.update(status='stopped' if isinstance(error, KeyboardInterrupt) else 'failed', error=repr(error))
            atomic_json(state_file, state)
            run_event(root, 'interrupted', error=repr(error), job=current)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--plan', action='store_true')
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--threads', type=int, default=5)
    parser.add_argument('--draws', type=int, default=128)
    parser.add_argument('--seconds', type=int, default=30)
    parser.add_argument('--cap', type=int, default=250000)
    args = parser.parse_args()
    if not 1 <= args.threads <= 128 or min(args.draws, args.seconds, args.cap) < 1:
        parser.error('positive draws/seconds/cap and 1..128 threads required')
    root = (args.resume or args.output or RESULTS / ('defect_mission_' + dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))).resolve()
    plan = json.loads((root / 'PLAN.json').read_text()) if args.resume else mission_plan(args.threads, args.draws, args.seconds, args.cap)
    if not args.run:
        print(json.dumps(dict(output=str(root), planned_jobs=len(plan['jobs']), threads=plan['threads'], cap=plan['cap'],
            seconds=plan['seconds'], jobs=plan['jobs'], inputs=len(plan['input_hashes']), limitations=plan['limitations']), indent=2))
        return
    run(root, plan, resume=bool(args.resume))
    print(root)


if __name__ == '__main__':
    main()
