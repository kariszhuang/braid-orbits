try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
import tempfile
import unittest

from compact_smooth_outputs import archive_plan, compact, restore, sha


class CompactOutputsTests(unittest.TestCase):
    def test_lossless_round_trip_and_nonoverwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp); root = base / 'root'; root.mkdir()
            folder = root / 'repeat'; folder.mkdir()
            a = folder / 'STATES.tsv'; a.write_bytes(b'table\n' * 100)
            b = folder / '._STATES.tsv'; b.write_bytes(b'resource fork bytes')
            old = {p: sha(p) for p in [a, b]}
            out = base / 'archive'
            result = compact(root, out, sorted(old), [folder])
            self.assertEqual(result['removed_files'], 2)
            self.assertFalse(folder.exists())
            restore(out / 'REPEATED_OUTPUTS.tar.gz', root)
            self.assertEqual(old, {p: sha(p) for p in old})
            a.write_bytes(b'user change')
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                restore(out / 'REPEATED_OUTPUTS.tar.gz', root)
            self.assertEqual(a.read_bytes(), b'user change')

    def test_unique_unverified_and_protected_certificates_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve(); first = root / 'guided_0'; first.mkdir()
            second = root / 'guided_1'; second.mkdir()
            third = root / 'guided_2'; third.mkdir()
            item = dict(degree=3, key=[[0, 1, 2]], folder='component',
                        summary=dict(status='closed_smooth'))
            for job in [first, second, third]:
                folder = job / 'component'; folder.mkdir()
                (folder / 'STATES.tsv').write_text('table\n')
                (folder / 'SUMMARY.tsv').write_text('summary\n')
                (job / 'CHECKS.json').write_text(json.dumps(dict(closed_identities=[item])))
            for job in [first, second]: (job / 'VERIFIED').write_text('yes')
            (root / 'NEW_COMPONENTS.json').write_text(json.dumps([dict(item, job=first.name)]))
            (root / 'PROFILE_INPUT.json').write_text(json.dumps(dict(components=[])))
            files, folders, summary = archive_plan(root)
            self.assertEqual(folders, [second / 'component'])
            self.assertEqual(summary['unique_new_certificates_preserved'], 1)
            self.assertTrue(all(p.is_relative_to(second) for p in files))
            files, folders, _ = archive_plan(root, [second / 'component/STATES.tsv'])
            self.assertEqual(files, [])
            self.assertEqual(folders, [])

    def test_finished_profile_copies_must_agree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'NEW_COMPONENTS.json').write_text('[]')
            (root / 'PROFILE_INPUT.json').write_text('{"components":[]}')
            p = root / 'profile_000000'; p.mkdir()
            (p / 'PROFILES.json').write_text('{"start":0,"next":1,"components":[{"genus":2}]}')
            (p / 'profile_000000.json').write_text('{"genus":3}')
            with self.assertRaisesRegex(ValueError, 'disagree'):
                archive_plan(root)
            (p / 'profile_000000.json').write_text('{"genus":2}')
            self.assertEqual(archive_plan(root)[2]['finished_profile_directories'], 1)

    def test_corrupt_archive_retains_restoration_targets(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp); root = base / 'root'; root.mkdir()
            source = root / 'table'; source.write_bytes(b'source')
            archive = base / 'archive'
            compact(root, archive, [source], [])
            (archive / 'REPEATED_OUTPUTS.tar.gz').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                restore(archive / 'REPEATED_OUTPUTS.tar.gz', root)
            self.assertFalse(source.exists())

    def test_restore_checks_existing_targets_before_any_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp); root = base / 'root'; root.mkdir()
            a = root / 'a'; a.write_bytes(b'first')
            b = root / 'b'; b.write_bytes(b'second')
            out = base / 'archive'
            compact(root, out, [a, b], [])
            b.write_bytes(b'user change')
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                restore(out / 'REPEATED_OUTPUTS.tar.gz', root)
            self.assertFalse(a.exists())
            self.assertEqual(b.read_bytes(), b'user change')


if __name__ == '__main__': unittest.main()
