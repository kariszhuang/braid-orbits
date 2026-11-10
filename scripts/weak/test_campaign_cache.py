"""Exact-cache and scheduling regressions; no research computations."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import ast
import csv
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import run_campaign as runner

SEED = '[(1 2)(3 4), (1 3)(2 4), (1 4)(2 3), ()]'
OTHER = '[(1 2 3), (1 3 2), (), ()]'


class CacheTests(unittest.TestCase):
    def test_periodic_report_never_rebuilds_all_coverage(self):
        tree = ast.parse(Path(runner.__file__).read_text())
        block = next(n for n in ast.walk(tree) if isinstance(n, ast.If)
                     and 'last_report' in ast.unparse(n.test) and '300' in ast.unparse(n.test))
        calls = []
        env = dict(time=SimpleNamespace(monotonic=lambda: 301), last_report=0,
                   live_index=SimpleNamespace(live_report=lambda: calls.append('live')),
                   exact_mission=True, discovery_refresh=lambda: calls.append('full'),
                   report=lambda root: calls.append('audit'), root=None)
        exec(compile(ast.Module(body=[block], type_ignores=[]), '<scheduler>', 'exec'), env)
        self.assertEqual(calls, ['live'])

    def attempt(self, root, name, seeds=(), statuses=(), verified=True):
        path = root / 'jobs' / name / 'attempt_000'
        path.mkdir(parents=True)
        (path / 'JOB_INPUT.json').write_text(json.dumps(dict(lane='test', kind='replay' if name=='replay' else 'sample')))
        (path / 'SUMMARY.txt').write_text('accepted generating seeds = 2\nsolved seeds = 1\n')
        (path / 'hard_candidates.txt').write_text(''.join(f'group=A_4 seed={s}\n' for s in seeds))
        (path / 'RESULTS.tsv').write_text('status\tseed\n' + ''.join(f'{status}\t{seed}\n' for status, seed in statuses))
        if verified:
            (path / 'VERIFIED').touch()
        return path

    def test_pending_matches_literal_reconciliation_and_replay(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = self.attempt(root, 'first', [SEED, OTHER])
            index = runner.VerifiedDiscoveryIndex(root)
            self.assertEqual(index.pending(), runner.candidates_from_files(index.paths, weak_only=True))
            second = self.attempt(root, 'second', [SEED], [('weak_found', SEED)])
            index.register(second)
            self.assertEqual(index.pending(), {'A_4': {OTHER}})
            self.assertEqual(index.pending(), runner.candidates_from_files(index.paths, weak_only=True))
            index.register(first)
            self.assertEqual(index.totals, [4, 2, 0])
            # Neither status nor replay touches previously closed raw files.
            with patch.object(Path, 'read_text', side_effect=AssertionError('raw reread')), \
                 patch.object(Path, 'open', side_effect=AssertionError('raw reopen')):
                self.assertEqual(index.pending(), {'A_4': {OTHER}})

    def test_active_attempt_is_never_parsed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self.attempt(root, 'active', verified=False)
            (path / 'hard_candidates.txt').write_text('group=')
            index = runner.VerifiedDiscoveryIndex(root)
            index.register(path)
            self.assertEqual(index.pending(), {})
            (path / 'VERIFIED').touch()
            with self.assertRaises(ValueError):
                index.register(path)

    def test_construction_survivors_and_imported_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self.attempt(root, 'construction')
            (path / 'construction_survivors.txt').write_text(f'group=A_4 seed={SEED}\n')
            (root / 'inputs').mkdir()
            (root / 'inputs' / 'RECOVERED_PENDING.txt').write_text(f'group=A_4 seed={OTHER}\n')
            index = runner.VerifiedDiscoveryIndex(root)
            self.assertEqual(index.pending(), {'A_4': {SEED, OTHER}})
            self.assertEqual(index.local_pending([path.parent], index.pending()), 1)

    def test_inherited_pending_is_associated_with_old_prefix_and_cached(self):
        with tempfile.TemporaryDirectory() as temporary:
            old = Path(temporary) / 'old'
            source = self.attempt(old, 'prefix', [SEED])
            root = Path(temporary) / 'new'
            (root / 'inputs').mkdir(parents=True)
            (root / 'inputs' / 'INHERITED_PENDING.txt').write_text(f'group=A_4 seed={SEED}\n')
            index = runner.VerifiedDiscoveryIndex(root)
            self.assertEqual(index.local_pending([source.parent], index.pending()), 1)
            with patch.object(Path, 'glob', side_effect=AssertionError('prefix rescan')):
                self.assertEqual(index.local_pending([source.parent], index.pending()), 1)

    def test_closed_witness_without_candidate_file_resolves_old_lead(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.attempt(root, 'first', [SEED])
            index = runner.VerifiedDiscoveryIndex(root)
            path = self.attempt(root, 'replay', statuses=[('weak_found', SEED)])
            index.register(path)
            self.assertEqual(index.pending(), {})

    def test_pending_callers_cannot_mutate_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.attempt(root, 'first', [SEED])
            index = runner.VerifiedDiscoveryIndex(root)
            index.pending()['A_4'].clear()
            self.assertEqual(index.pending(), {'A_4': {SEED}})

    def test_positive_cache_eviction_never_resurrects_old_pending(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.attempt(root, 'first', [SEED])
            index = runner.VerifiedDiscoveryIndex(root)
            index._resolved_limit = 1
            solved = self.attempt(root, 'replay', statuses=[('weak_found', SEED), ('weak_found', OTHER)])
            index.register(solved)
            self.assertEqual(index.pending(), {})
            # A forgotten positive may require replay when a NEW source repeats
            # it; no unknown candidate is suppressed to save memory.
            later = self.attempt(root, 'later', [SEED])
            index.register(later)
            self.assertIn(index.pending(), ({}, {'A_4': {SEED}}))

    def test_real_exact_cache_keeps_binary_and_gap_guards(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = dict(name='p', kind='passport', passport='3/3/3/3')
            folder = root / 'p'
            cache = runner.VerifiedCoverageCache()
            self.assertEqual(cache.get(row, folder, 'sha')['next'], 0)
            attempt = folder / 'attempt_000'
            attempt.mkdir(parents=True)
            (attempt / 'VERIFIED').touch()
            (attempt / 'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='sha', passport=row['passport'], enumeration='centralizer-v1')))
            summary = attempt / 'SUMMARY.txt'
            summary.write_text('scan start = 0\nscan next = 4\nscan total candidate pairs = 8\nscan enumeration = centralizer-v1\n')
            cache.invalidate(folder)
            self.assertEqual(cache.get(row, folder, 'sha')['next'], 4)
            with self.assertRaises(RuntimeError):
                cache.get(row, folder, 'other-binary')
            summary.write_text(summary.read_text().replace('start = 0', 'start = 1'))
            cache.invalidate(folder)
            with self.assertRaises(RuntimeError):
                cache.get(row, folder, 'sha')

    def test_indexed_exact_coverage_matches_raw_without_raw_rereads(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = dict(name='p', kind='passport', lane='exact', passport='3/3/3/3')
            folder = root / 'jobs' / 'p'
            for i in range(5):
                attempt = folder / f'attempt_{i:03d}'
                attempt.mkdir(parents=True)
                (attempt / 'VERIFIED').touch()
                (attempt / 'JOB_INPUT.json').write_text(json.dumps(row))
                (attempt / 'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='sha', passport=row['passport'], enumeration='centralizer-v1')))
                (attempt / 'SUMMARY.txt').write_text(f'scan start = {i*2}\nscan next = {(i+1)*2}\nscan total candidate pairs = 10\nscan enumeration = centralizer-v1\nhard unsolved candidates = {int(i==2)}\n')
            expected = runner.exact_coverage(folder, 'sha', row['passport'])
            index = runner.VerifiedDiscoveryIndex(root)
            cache = runner.VerifiedCoverageCache(index)
            with patch.object(Path, 'read_text', side_effect=AssertionError('raw reread')), \
                 patch.object(Path, 'open', side_effect=AssertionError('raw reopen')):
                self.assertEqual(cache.get(row, folder, 'sha'), expected)
                self.assertFalse(cache.get(row, folder, 'sha')['search_complete'])
                with self.assertRaises(RuntimeError):
                    cache.get(row, folder, 'wrong-sha')

    def test_coverage_cache_revalidates_only_changed_job_and_contract(self):
        cache = runner.VerifiedCoverageCache()
        row = dict(name='p', kind='passport', passport='3/3/3/3')
        with patch.object(runner, 'discovery_coverage', return_value={'next': 10}) as loader:
            self.assertEqual(cache.get(row, Path('job'), 'sha'), {'next': 10})
            cache.get(row, Path('job'), 'sha')['next'] = 999
            self.assertEqual(cache.get(row, Path('job'), 'sha')['next'], 10)
            self.assertEqual(loader.call_count, 1)
            cache.invalidate(Path('job'))
            cache.get(row, Path('job'), 'sha')
            self.assertEqual(loader.call_count, 2)
            cache.get(dict(row, passport='4/4/4/4'), Path('job'), 'sha')
            cache.get(row, Path('job'), 'changed-sha')
            self.assertEqual(loader.call_count, 4)

    def test_progress_build_uses_one_coverage_lookup_per_row(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            index = runner.VerifiedDiscoveryIndex(root)
            row = dict(name='p', kind='passport', lane='exact', group='S', d=9,
                       genus=1, passport='3/3/3/3')
            coverage = dict(next=4, total=4, enumeration_complete=True,
                            traversal_exhausted=True, search_complete=True)
            with patch.object(runner, 'discovery_coverage', return_value=coverage) as loader:
                cache = runner.VerifiedCoverageCache()
                runner.write_discovery_progress(root, {'lanes': {'exact': [row]}}, 'sha', cache, index, True)
                self.assertEqual(loader.call_count, 1)
                with (root / 'EXACT_PASSPORT_PROGRESS.tsv').open() as stream:
                    actual = list(csv.DictReader(stream, delimiter='\t'))
                self.assertEqual(actual[0]['search_complete'], 'True')
                runner.write_discovery_progress(root, {'lanes': {'exact': [row]}}, 'sha', cache, index, True)
                self.assertEqual(loader.call_count, 1)

    def test_indexed_audit_report_matches_raw_report_without_summary_reads(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.attempt(root, 'first', [SEED], [('weak_found', SEED)])
            index = runner.VerifiedDiscoveryIndex(root)
            runner.report(root)
            expected = {n: (root / n).read_text() for n in ('COVERAGE.tsv', 'PASSPORT_SCANS.tsv', 'CATALOG.tsv')}
            original = Path.read_text
            def checked(path, *args, **kwargs):
                if path.name == 'SUMMARY.txt':
                    raise AssertionError('audit reread closed summary')
                return original(path, *args, **kwargs)
            with patch.object(Path, 'read_text', checked):
                runner.report(root, index=index)
            self.assertEqual(expected, {n: (root / n).read_text() for n in expected})


if __name__ == '__main__':
    unittest.main()
