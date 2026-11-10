"""Recovery journal safety tests; no live campaign writes."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import tempfile
import json
import unittest
from pathlib import Path
from unittest.mock import patch

from recover_genus_three import HEADER, durable_journal, inactive,recover
import run_genus_three as runner


class RecoveryTests(unittest.TestCase):
    def attention_fixture(self,root):
        (root/'jobs/job_0_00000000_14_collision').mkdir(parents=True);(root/'banks').mkdir()
        (root/'jobs/job_0_00000000_14_collision/LEADS.tsv').write_bytes(b'')
        (root/'jobs/job_0_00000000_14_collision/RESOLUTIONS.tsv').write_bytes(b'')
        (root/'banks/words_00000000.txt').write_text('1,4\n')
        (root/'COMPLETED.jsonl').write_bytes(b'');(root/'smooth_negative_audit').write_bytes(b'fixture auditor')
        runner.write(root/'LAUNCH.json',dict(pid=123))
        runner.write(root/'PLAN.json',dict(workers=1,streams=1,cap=10,closure_seconds=1,verifier='macro',deadline_epoch=None))
        runner.write(root/'STATE.json',dict(status='attention',reason='audit_or_producer_failure',epochs=[1],active_jobs=[],boundary_cursor=504000))
        runner.write(root/'FINISHED.json',dict(status='attention',reason='audit_or_producer_failure',epochs=[1],completed_jobs=0))
        runner.write(root/'CHECKPOINT.json',dict(next_epochs=[1],completed_jobs=0,counts={}))
        runner.write(root/'KNOWN_IDS.json',[]);runner.write(root/'BINDINGS.json',{'KNOWN_IDS.json':runner.sha(root/'KNOWN_IDS.json')})
        runner.write(root/'INHERITED_CHECKPOINT.json',dict(allocation_origin=[0],allocation_weights=None))

    def test_attention_recovery_needs_explicit_mode_and_preserves_all_original_receipts(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.attention_fixture(root)
            original={name:(root/name).read_bytes()for name in('STATE.json','FINISHED.json','CHECKPOINT.json')}
            with patch('recover_genus_three.inactive'):
                with self.assertRaisesRegex(ValueError,'finalized'):recover(root)
            with patch('recover_genus_three.inactive'),patch('recover_genus_three.audit_guided_native',return_value={'independent_native':{'nonsmooth':0}}):
                record=recover(root,finalized_attention=True)
            self.assertTrue(record['passed']);self.assertEqual(record['pending'],0)
            for name,data in original.items():
                self.assertEqual((root/name).read_bytes(),data)
                self.assertEqual((root/'RECOVERY'/('ORIGINAL_'+name)).read_bytes(),data)
            adapter=Path(record['admissible_predecessor']);self.assertEqual(adapter,(root/'RECOVERY/ADMISSIBLE').resolve())
            self.assertEqual(runner.finalized_predecessor(adapter)[0]['status'],'stopped')
            self.assertEqual(json.loads((adapter/'CHECKPOINT.json').read_text())['next_epochs'],[1])
            (root/'FINISHED.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'recovery'):runner.finalized_predecessor(adapter)

    def test_attention_mode_cannot_accept_running_or_wrong_final_reason(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.attention_fixture(root)
            runner.write(root/'STATE.json',dict(status='running'))
            with patch('recover_genus_three.inactive'):
                with self.assertRaisesRegex(ValueError,'attention'):recover(root,finalized_attention=True)
    def fixture(self, lead, resolution):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        folder = Path(tmp.name)
        (folder/'LEADS.tsv').write_text(lead)
        (folder/'RESOLUTIONS.tsv').write_text(resolution)
        return folder

    def test_partial_lines_are_preserved_not_admitted(self):
        folder = self.fixture('3\tseed\n4\tpartial', HEADER+'3\tnonsmooth\t4\t\tseed\n4\tpart')
        self.assertEqual(durable_journal(folder), ([(3, 'seed')], [['3', 'nonsmooth', '4', '', 'seed']]))

    def test_empty_and_partial_header_retire_nothing(self):
        for value in ['', 'candidate\tsta']:
            self.assertEqual(durable_journal(self.fixture('', value)), ([], []))

    def test_malformed_complete_header_fails(self):
        with self.assertRaises(ValueError):
            durable_journal(self.fixture('', 'wrong\n'))

    def test_wrong_seed_or_unlogged_resolution_fails(self):
        for lead in ['', '3\tother\n']:
            with self.assertRaises(ValueError):
                durable_journal(self.fixture(lead, HEADER+'3\tnonsmooth\t4\t\tseed\n'))

    def test_nonmonotone_leads_fail(self):
        with self.assertRaises(ValueError):
            durable_journal(self.fixture('3\tseed\n3\tseed\n', HEADER))

    def test_existing_controller_blocks_without_signal(self):
        with patch('recover_genus_three.os.kill') as kill:
            with self.assertRaisesRegex(RuntimeError, 'controller PID still exists'):
                inactive(Path('/safe/fixture'), 123)
            kill.assert_called_once_with(123, 0)


if __name__ == '__main__':
    unittest.main()
