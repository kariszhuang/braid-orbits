try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import tempfile
import unittest

from recover_forest_interruption import ensure_inactive, negative_prefix, excluded_ids, save, sha, recovered_epochs,preserve_original_state
from test_guided_native import fixture, BINARY


class RecoveryTests(unittest.TestCase):
    def test_original_state_snapshot_cannot_be_overwritten_by_retry(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);save(root/'STATE.json',dict(status='failed',workers=0))
            original=(root/'STATE.json').read_bytes();preserve_original_state(root)
            (root/'STATE.json').write_text('{"status":"stopped"}\n')
            with self.assertRaises(FileExistsError):preserve_original_state(root)
            self.assertEqual((root/'STATE_BEFORE_RECOVERY.json').read_bytes(),original)

    def test_disk_guard_state_recovers_missing_epochs_from_job_names(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'guided_0_0000042_17_forest').mkdir()
            (root/'guided_2_0000070_12_mutate').mkdir()
            save(root/'PROFILE_SUMMARY.json',dict(checked=5595))
            self.assertEqual(recovered_epochs(root,dict(workers=6,mechanism_search=True),dict(status='failed',workers=0)),[43,0,71,0,0,5595])

    def test_active_controller_is_never_finalized(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(RuntimeError,'PID still exists'):ensure_inactive(Path(temp),os.getpid())

    def test_empty_unflushed_resolution_journal_retires_no_starts(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);folder=root/'job';folder.mkdir();fixture(folder)
            (folder/'RESOLUTIONS.tsv').write_bytes(b'')
            original=(folder/'LEADS.tsv').read_bytes()
            record=negative_prefix(folder,root/'negative',9,BINARY)
            self.assertEqual(record['records'],0)
            self.assertTrue(record['trailing_incomplete_resolution'])
            self.assertEqual((folder/'LEADS.tsv').read_bytes(),original)
            save(root/'RECOVERY.json',dict(negative_prefixes={folder.name:record}))
            self.assertEqual(excluded_ids(root,folder),set())

    def test_recover_only_verified_negative_prefix_and_bind_sources(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);folder=root/'job';folder.mkdir();fixture(folder)
            lead=(folder/'LEADS.tsv').read_text();seed=lead.split('\t',1)[1]
            with (folder/'LEADS.tsv').open('a') as f:f.write('2\t'+seed)
            with (folder/'RESOLUTIONS.tsv').open('a') as f:f.write('2\tnonsm')
            record=negative_prefix(folder,root/'negative',9,BINARY)
            self.assertEqual(record['records'],1);self.assertTrue(record['trailing_incomplete_resolution'])
            save(root/'RECOVERY.json',dict(negative_prefixes={folder.name:record}))
            self.assertEqual(excluded_ids(root,folder),{'1'})
            p=folder/'LEADS.tsv';p.write_text(p.read_text()+'3\t'+seed)
            with self.assertRaisesRegex(ValueError,'journal changed'):excluded_ids(root,folder)

    def test_false_witness_and_reordered_seed_fail_closed(self):
        for false_word in [True,False]:
            with tempfile.TemporaryDirectory() as temp:
                root=Path(temp);folder=root/'job';folder.mkdir();fixture(folder)
                if false_word:
                    p=folder/'RESOLUTIONS.tsv';p.write_text(p.read_text().replace('\tnonsmooth\t4\t','\tnonsmooth\t0,3\t'))
                    with self.assertRaises(Exception):negative_prefix(folder,root/'negative',9,BINARY)
                else:
                    p=folder/'LEADS.tsv';p.write_text(p.read_text().replace('1\t','2\t',1))
                    with self.assertRaisesRegex(ValueError,'mismatch'):negative_prefix(folder,root/'negative',9,BINARY)

    def test_audited_prefix_cannot_change_after_verification(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);folder=root/'job';folder.mkdir();fixture(folder)
            record=negative_prefix(folder,root/'negative',9,BINARY)
            save(root/'RECOVERY.json',dict(negative_prefixes={folder.name:record}))
            p=root/'negative/RESOLUTIONS.tsv';p.write_text(p.read_text().replace('\tnonsmooth\t4\t','\tnonsmooth\t0,3\t'))
            with self.assertRaisesRegex(ValueError,'prefix changed'):excluded_ids(root,folder)


if __name__=='__main__':unittest.main()
