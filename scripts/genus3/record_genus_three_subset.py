"""Bind the C215 implementation, recovery, controls and paired measurement."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from run_genus_three import PROJECT, sha, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['build', 'benchmark', 'recovered', 'temporary', 'successor', 'out']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    validation = json.loads((args.build/'VALIDATION.json').read_text())
    benchmark = json.loads((args.benchmark/'CHECKS.json').read_text())
    recovery = json.loads((args.recovered/'RECOVERY/CHECKS.json').read_text())
    if not all(value['passed'] for value in [validation, benchmark, recovery]):
        raise ValueError('an admission receipt did not pass')
    if benchmark['screen'] != 'forest':
        raise ValueError('not a whole-subset comparison')
    for name, key in [('wide', 'binary_sha256'), ('compact', 'compact_binary_sha256'), ('auditor', 'audit_binary_sha256')]:
        if sha(args.build/'bin'/name) != validation[key]:
            raise ValueError('validated build changed')
    for name in ['src/smooth_cyclic_screen.rs', 'src/smooth_genus_three.rs']:
        if sha(PROJECT/name) != validation['compiled_source_files'][name]:
            raise ValueError('scientific source changed after validation')
    for name, digest in json.loads((args.successor/'BINDINGS.json').read_text()).items():
        if sha(args.successor/name) != digest:
            raise ValueError('frozen successor binding changed: '+name)
    temporary = json.loads((args.temporary/'FINISHED.json').read_text())
    if temporary['status'] != 'stopped' or temporary['counts']['unresolved']:
        raise ValueError('temporary predecessor did not finalize safely')
    tests = []
    for filename in ['test_recover_genus_three.py', 'test_genus_three.py',
                     'test_genus_three_boundary.py', 'test_watch_genus_three.py']:
        run = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', 'scripts/genus3',
            '-p', filename], cwd=PROJECT, capture_output=True, text=True, check=True,
            env=dict(os.environ, G3_CYCLIC_VALIDATION=str(args.build.resolve())))
        if 'skipped' in run.stderr:
            raise ValueError('required lifecycle/native test was skipped')
        tests.append(dict(file=filename, output=run.stdout+run.stderr))
    current = json.loads((PROJECT/'.run-control/genus-three-current.json').read_text())
    if Path(current['root']) != args.successor.resolve():
        raise ValueError('successor is not current')
    status = json.loads(subprocess.check_output([sys.executable, PROJECT/'scripts/genus3/watch_genus_three.py',
        '--json'], text=True))
    if status['status'] != 'RUNNING' or not 1 <= status['processes']['search_workers'] <= 6:
        raise ValueError('owned current six-slot continuation is not running')
    files = ['scripts/genus3/recover_genus_three.py', 'scripts/genus3/test_recover_genus_three.py',
        'scripts/genus3/benchmark_genus_three_cyclic.py', 'scripts/genus3/record_genus_three_subset.py',
        'src/smooth_cyclic_screen.rs', 'src/smooth_genus_three.rs']
    result = dict(claim='C215', passed=True, quotient='simultaneous S_d conjugacy',
        recovery=dict(path=str(args.recovered/'RECOVERY/CHECKS.json'), sha256=sha(args.recovered/'RECOVERY/CHECKS.json'),
            completed_jobs=recovery['completed_jobs'],
            durable_interrupted_leads=sum(row['durable_leads'] for row in recovery['interrupted']),
            pending_recovered=sum(row['recovered_pending'] for row in recovery['interrupted']),
            pending=recovery['pending']),
        lifecycle_tests=tests,
        validation=dict(path=str(args.build/'VALIDATION.json'), sha256=sha(args.build/'VALIDATION.json'),
            cyclic_tests=(args.build/'logs/cyclic_regressions.out').read_text()),
        benchmark=dict(path=str(args.benchmark/'CHECKS.json'), sha256=sha(args.benchmark/'CHECKS.json'),
            summaries=benchmark['summaries'],
            forest_tested=sum(int(pair['results']['cyclic']['forest']['tested']) for pair in benchmark['paired']),
            forest_rejected=sum(int(pair['results']['cyclic']['forest']['rejected']) for pair in benchmark['paired']),
            qualification='matched fixed proposals/RNG; bounded throughput, not discovery probability'),
        source_sha256={name:sha(PROJECT/name) for name in files},
        scope='exact literal negative-only opt-in whole-subset cyclic sieve; no geometry rule or global nonexistence',
        admission='opt-in only; no production activation without a measured whole-job benefit',
        no_new_discovery_claim=True)
    write(args.out/'CHECKS.json', result)
    write(args.out/'HANDOFF.json', dict(passed=True, previous=str(args.recovered),
        temporary=str(args.temporary), successor=str(args.successor), live_snapshot=status,
        preserved_epochs=json.loads((args.successor/'INHERITED_CHECKPOINT.json').read_text())['next_epochs'],
        bindings_sha256=sha(args.successor/'BINDINGS.json'),
        qualification='launch/status snapshot, not a future uptime guarantee; one shared pool; no GPU'))
    print(json.dumps(dict(passed=True, out=str(args.out))))


if __name__ == '__main__':
    main()
