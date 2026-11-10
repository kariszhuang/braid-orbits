try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from extend_public_smooth_catalog import extend, merge, sha, verify_bindings


class MergeTests(unittest.TestCase):
    def entry(self, identity='a'):
        return dict(id=identity, degree=9, genus=2, states=12, group_order=324)

    def test_existing_representative_and_metadata_are_preserved(self):
        old = dict(self.entry(), tuple=['old'], revision=4, last_updated='old')
        result, added = merge(dict(component_count=1, components=[old]),
                              [dict(self.entry(), tuple=['new'])])
        self.assertEqual(result['components'], [old])
        self.assertEqual(added, [])

    def test_addition_and_duplicate_rejection(self):
        old = dict(component_count=1, components=[self.entry()])
        result, added = merge(old, [self.entry('b')])
        self.assertEqual(result['component_count'], 2)
        self.assertEqual(added, ['b'])
        with self.assertRaises(ValueError):
            merge(old, [self.entry('b'), self.entry('b')])

    def test_conflicting_identity_fails_closed(self):
        old = dict(component_count=1, components=[self.entry()])
        with self.assertRaises(ValueError):
            merge(old, [dict(self.entry(), genus=3)])

    def test_invalid_baseline_count_fails_closed(self):
        with self.assertRaises(ValueError):
            merge(dict(component_count=2, components=[self.entry()]), [])


class SnapshotTests(unittest.TestCase):
    def test_running_campaign_requires_separate_frozen_index(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'STATE.json').write_text(json.dumps(dict(status='searching')))
            live = root / 'NEW_COMPONENTS.json'
            live.write_text('[]')
            previous = dict(component_count=0, components=[])
            with self.assertRaisesRegex(ValueError, 'inactive'):
                extend(root, previous, root / 'out.json')
            with self.assertRaisesRegex(ValueError, 'separate'):
                extend(root, previous, root / 'out.json', live)
            snapshot = root / 'frozen.json'
            snapshot.write_text('[]')
            # The running index may evolve; only the explicitly frozen input is read.
            live.write_text('not the frozen index')
            checks = extend(root, previous, root / 'out.json', snapshot)
            self.assertTrue(checks['frozen_discovery_snapshot'])
            self.assertEqual(checks['state_at_export']['status'], 'searching')
            self.assertEqual(checks['discovery_index_sha256'], sha(snapshot))
            self.assertEqual(checks['added'], 0)

    def test_changed_snapshot_fails_before_output(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'STATE.json').write_text(json.dumps(dict(status='searching')))
            snapshot = root / 'frozen.json'
            snapshot.write_text('[]')
            def changed(previous, additions):
                snapshot.write_text('[{}]')
                return merge(previous, additions)
            with patch('extend_public_smooth_catalog.merge', side_effect=changed):
                with self.assertRaisesRegex(ValueError, 'index changed'):
                    extend(root, dict(component_count=0, components=[]), root / 'out.json', snapshot)
            self.assertFalse((root / 'out.json').exists())

    def test_certificate_recheck_detects_changed_input(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'certificate'
            path.write_text('original')
            reference = dict(path=str(path), sha256=sha(path))
            bindings = [dict(state_table=reference, summary=reference, proof=reference)]
            verify_bindings(bindings)
            path.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'certificate input changed'):
                verify_bindings(bindings)


if __name__ == '__main__':
    unittest.main()
