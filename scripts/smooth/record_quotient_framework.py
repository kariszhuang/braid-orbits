"""Curate source-bound quotient/pilot/cleanup evidence without starting a search."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys

from compact_smooth_outputs import identity, sha
from investigate_smooth_quotients import index, index_regime
from verify_smooth_geometry import parse_tuple


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'analysis', 'cleanup', 'out']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--replace-unpublished', action='store_true')
    p.add_argument('--tests', action='store_true', help='run focused safety, quotient and metadata regressions')
    args = p.parse_args(); root = args.root.resolve()
    if args.out.exists():
        manifest = json.loads((Path(__file__).resolve().parents[2] / 'evidence/manifest.json').read_text())
        project = Path(__file__).resolve().parents[2]
        if not args.replace_unpublished or any((project / e['path']).resolve().is_relative_to(args.out.resolve()) for e in manifest['sources']):
            raise ValueError('do not overwrite published evidence')
    args.out.mkdir(exist_ok=True)
    state = json.loads((root / 'STATE.json').read_text())
    if state['status'] != 'failed' or state['workers'] != 0:
        raise ValueError('expected preserved failed/inactive run')
    launch = json.loads((root / 'LAUNCH.json').read_text())
    try:
        os.kill(launch['pid'], 0)
    except ProcessLookupError:
        pass
    else:
        raise ValueError('controller alive or PID reused')
    queue = json.loads((root / 'PROFILE_INPUT.json').read_text())['components']
    new_entries = json.loads((root / 'NEW_COMPONENTS.json').read_text())
    old_keys = {identity(e) for e in queue if 'component_identity_key' in e}
    new_keys = {identity(e) for e in new_entries}
    if len(new_keys) != len(new_entries) or old_keys & new_keys:
        raise ValueError('new component identity overlap')
    jobs = {}; sources = []; lanes = Counter(); degrees = Counter(); genera = Counter(); regimes = Counter(); mass = 0
    for e in new_entries:
        job = root / e['job']; path = job / 'CHECKS.json'
        if e['job'] not in jobs:
            if not (job / 'VERIFIED').is_file():
                raise ValueError('new component in an unverified job')
            jobs[e['job']] = json.loads(path.read_text())
        report = jobs[e['job']]
        matches = [r for r in report['closed_identities'] if identity(r) == identity(e) and r['folder'] == e['folder']]
        if len(matches) != 1 or matches[0]['summary'] != e['summary']:
            raise ValueError('component index/check mismatch')
        if not (job / e['folder'] / 'STATES.tsv').is_file():
            raise ValueError('unique certificate was not preserved loose')
        lane = report['construction']['lane']; lanes[lane] += 1
        degrees[e['degree']] += 1; genera[int(e['summary']['genus'])] += 1
        mass += int(e['summary']['states'])
        t = parse_tuple(e['summary']['seed'], e['degree'])
        regimes[index_regime([index(g) for g in t], e['degree'])] += 1
    jobs_by_lane = Counter(); closed_by_lane = Counter(); unresolved_by_lane = Counter()
    for job in sorted(root.iterdir()):
        if not job.is_dir() or not (job / 'VERIFIED').is_file() or not (job / 'CHECKS.json').is_file():
            continue
        path = job / 'CHECKS.json'; report = json.loads(path.read_text())
        c = report['construction']; lane = c['lane']
        jobs_by_lane[lane] += 1; closed_by_lane[lane] += int(c['closed'])
        unresolved_by_lane[lane] += int(c['unresolved'])
        sources.append(dict(path=str(path), sha256=sha(path)))
    runtime = dict(root=str(root), state=state, controller_alive=False,
        index_distinct=len(new_entries), by_lane=dict(lanes), by_degree=dict(degrees), by_genus=dict(genera),
        states=mass, previous_profile_population=len(queue), recorded_component_total=len(queue)+len(new_entries),
        index_budget_regimes=dict(regimes),
        verified_jobs=sum(jobs_by_lane.values()), jobs_by_lane=dict(jobs_by_lane),
        closed_occurrences_by_lane=dict(closed_by_lane), unresolved_occurrences_by_lane=dict(unresolved_by_lane),
        job_check_bindings=sources,
        limits='Recorded complete identities, not distinct mechanisms or geometric classifications. '
        'The 1574 live discoveries are outside the finished frozen 5595-profile queue. '
        'Interrupted unverified jobs and all capped starts remain unresolved. '
        'Explicit queue identity overlap checked; initial curated controls were retired by the runner bound known-state bank.')
    save(args.out / 'RUNTIME.json', runtime)
    qpath = args.analysis / 'QUOTIENT_GENERA.json'; q = json.loads(qpath.read_text())
    if not q['passed'] or q['components'] != len(queue) or q['input_queue_sha256'] != sha(root / 'PROFILE_INPUT.json'):
        raise ValueError('quotient audit binding')
    q['full_report'] = dict(path=str(qpath), sha256=sha(qpath))
    q.pop('source_reports')
    save(args.out / 'QUOTIENT_GENERA.json', q)
    pilots = {}
    for name in ['TWO_ACTIVE_BLOCKS', 'GENUS_ZERO_QUOTIENT']:
        path = args.analysis / (name + '.json'); result = json.loads(path.read_text())
        audit = result['independent_negative_audit']
        if audit['records'] != result['counts']['nonsmooth'] or audit['unresolved'] or audit['closed']:
            raise ValueError('independent pilot negative audit')
        result['word_length_counts'] = dict(Counter(len(e['word']) for e in result['results']))
        result['source'] = dict(path=str(path), sha256=sha(path))
        pilots[name] = result
    save(args.out / 'PILOTS.json', pilots)
    cpath = args.cleanup / 'RESULTS.json'; c = json.loads(cpath.read_text()); c.pop('members')
    if not c['all_members_verified_before_removal'] or sha(Path(c['archive'])) != c['archive_sha256']:
        raise ValueError('cleanup recovery binding')
    c['net_free_increase_bytes'] = c['free_bytes_after'] - c['free_bytes_before']
    c['source'] = dict(path=str(cpath), sha256=sha(cpath))
    c['recovery'] = dict(path=str(args.cleanup / 'RECOVERY.json'), sha256=sha(args.cleanup / 'RECOVERY.json'))
    save(args.out / 'CLEANUP.json', c)
    tests = []
    if args.tests:
        project = Path(__file__).resolve().parents[2]
        commands = [
            [sys.executable, 'scripts/smooth/test_compact_smooth_outputs.py'],
            [sys.executable, 'scripts/smooth/test_investigate_smooth_quotients.py'],
            ['uv', 'run', '--offline', '--isolated', '--no-project', '--python', '3.13', '--with', 'sympy==1.14.0', 'python', 'scripts/smooth/test_mechanism_breaker.py'],
            [sys.executable, '-c', 'import sys,tempfile,runpy; sys.path.insert(0,"scripts/smooth"); tempfile.tempdir="data"; runpy.run_path("scripts/smooth/test_compact_smooth_outputs.py",run_name="__main__")'],
        ]
        for command in commands:
            result = subprocess.run(command, cwd=project, capture_output=True, text=True)
            tests.append(dict(command=command, exit_code=result.returncode, output=result.stdout+result.stderr))
            if result.returncode:
                save(args.out / 'TESTS.json', dict(passed=False, runs=tests))
                raise ValueError('focused regression failed')
        save(args.out / 'TESTS.json', dict(passed=True, runs=tests,
            unique_checks=19, t7_round_trip_safety_checks_repeated=5,
            source_hashes={name: sha(_scripts_path.find(name)) for name in [
                'compact_smooth_outputs.py', 'test_compact_smooth_outputs.py',
                'investigate_smooth_quotients.py', 'test_investigate_smooth_quotients.py',
                'profile_mechanism_breaker.py', 'test_mechanism_breaker.py', 'record_quotient_framework.py']}))
    binding_paths = [root / n for n in ['STATE.json', 'FAILURE.json', 'LAUNCH.json', 'NEW_COMPONENTS.json', 'PROFILE_INPUT.json', 'PROFILE_SUMMARY.json', 'KNOWN_CACHE.json']]
    binding_paths += [qpath, args.analysis / 'TWO_ACTIVE_BLOCKS.json', args.analysis / 'GENUS_ZERO_QUOTIENT.json', cpath, args.cleanup / 'RECOVERY.json']
    if tests:
        binding_paths.append(args.out / 'TESTS.json')
    save(args.out / 'CHECKS.json', dict(passed=True, controller_stopped=True, no_long_run_started=True,
        quotient_components=q['components'], pilot_native_negatives=sum(e['independent_negative_audit']['records'] for e in pilots.values()),
        bindings=[dict(path=str(path), sha256=sha(path)) for path in binding_paths],
        proof_review='human review pending; no universal classification or geometric/literature claim'))
    print('RECORDED', len(new_entries), 'new identities;', q['components'], 'quotient audits;', c['net_free_increase_bytes'], 'bytes net free')


if __name__ == '__main__':
    main()
