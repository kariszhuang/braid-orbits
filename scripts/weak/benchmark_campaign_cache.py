#!/usr/bin/env python3
"""Paired old/full-scan versus RAM-cached progress export; synthetic exact files.

No discovery computations. Baseline is the previous committed real refresh
function, not a hand-written estimate. Equality checks cover both progress tables.
The interval describes this technical workload, not a discovery-rate multiplier.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import ast
import csv
import itertools
import json
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile
import time
from unittest.mock import patch

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
import run_campaign as runner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, default=1000)
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--baseline-ref', default='6dbdaddb3b16a5987a6a49f39a7785d832184de7')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if not 10 <= args.rows <= 40000 or not 3 <= args.repeats <= 25:
        parser.error('rows must be 10..40000; repeats 3..25')
    old = subprocess.check_output(['git', 'show', f'{args.baseline_ref}:run_campaign.py'], cwd=PROJECT, text=True)
    node = next(n for n in ast.walk(ast.parse(old)) if isinstance(n, ast.FunctionDef) and n.name == 'discovery_refresh')
    with tempfile.TemporaryDirectory(prefix='hurwitz-progress-benchmark-') as temporary:
        root = Path(temporary)
        rows = []
        for i in range(args.rows):
            row = dict(name=f'passport_{i}', kind='passport', lane='exact', group='S', d=9,
                       genus=1, passport='3/3/3/3')
            rows.append(row)
            folder = root / 'jobs' / row['name'] / 'attempt_000'
            folder.mkdir(parents=True)
            (folder / 'VERIFIED').touch()
            (folder / 'JOB_INPUT.json').write_text(json.dumps(row))
            (folder / 'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='fixture', passport=row['passport'], enumeration='centralizer-v1')))
            (folder / 'SUMMARY.txt').write_text('scan start = 0\nscan next = 10\nscan total candidate pairs = 10\nscan enumeration = centralizer-v1\nhard unsolved candidates = 0\n')
        schedule = dict(lanes={'exact': rows})
        index = runner.VerifiedDiscoveryIndex(root)
        cache = runner.VerifiedCoverageCache(index)
        env = dict(root=root, schedule=schedule, digest='fixture', csv=csv, itertools=itertools,
                   Path=Path, exact_mission=True, live_index=index, atomic_json=runner.atomic_json,
                   discovery_coverage=runner.discovery_coverage, discovery_prefixes=runner.discovery_prefixes,
                   pending_paths=lambda: [], candidates_from_files=runner.candidates_from_files)
        exec(compile(ast.Module(body=[node], type_ignores=[]), '<committed-refresh>', 'exec'), env)
        baseline = env['discovery_refresh']
        optimized = lambda: runner.write_discovery_progress(root, schedule, 'fixture', cache, index, True)
        baseline()
        names = ('DISCOVERY_PROGRESS.tsv', 'EXACT_PASSPORT_PROGRESS.tsv')
        expected = {name: (root / name).read_bytes() for name in names}
        optimized()  # validate/warm each historical row once
        assert expected == {name: (root / name).read_bytes() for name in names}
        trials = []
        rng = random.Random(20261010)
        for repeat in range(args.repeats):
            arms = [('full_scan', baseline), ('ram_cache', optimized)]
            rng.shuffle(arms)
            values = {}
            for label, function in arms:
                opens = [0]
                original = Path.open
                def counted(path, *a, **kw):
                    opens[0] += 1
                    return original(path, *a, **kw)
                with patch.object(Path, 'open', counted):
                    begin = time.perf_counter()
                    function()
                    elapsed = time.perf_counter() - begin
                assert expected == {name: (root / name).read_bytes() for name in names}
                values[label] = dict(seconds=elapsed, path_open_calls=opens[0])
            values['ratio'] = values['full_scan']['seconds'] / values['ram_cache']['seconds']
            trials.append(values)
        ratios = [t['ratio'] for t in trials]
        bootstrap = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(4000))
        result = dict(rows=args.rows, repeats=args.repeats, equal_progress_tables=True,
                      full_scan_seconds_median=statistics.median(t['full_scan']['seconds'] for t in trials),
                      cached_seconds_median=statistics.median(t['ram_cache']['seconds'] for t in trials),
                      paired_ratio_median=statistics.median(ratios),
                      paired_ratio_bootstrap_95_percent=[bootstrap[100], bootstrap[3899]],
                      cache_hits=cache.hits, cache_misses=cache.misses, trials=trials,
                      caveat='Synthetic warmed internal-disk progress export; not live T7 timing or discovery-rate acceleration.')
        args.out.parent.mkdir(parents=True, exist_ok=True)
        runner.atomic_json(args.out, result)
        print(json.dumps({k: v for k, v in result.items() if k != 'trials'}, indent=2))


if __name__ == '__main__':
    main()
