try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from audit_guided_native import audit_guided_native
from test_guided_native import fixture,seed_text

PROJECT=Path(__file__).resolve().parents[2]
PRODUCER=Path(os.environ.get('DESIGN_BINARY',PROJECT/'target/discovery-design/release/smooth_forest_fast'))
AUDITOR=Path(os.environ.get('GUIDED_AUDIT_BINARY',PROJECT/'target/discovery-design/release/smooth_negative_audit'))

class DiscoveryDesignTests(unittest.TestCase):
    def test_negative_bank_literal_corruption_fails_in_both_engines(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);bank=root/'bank.tsv';bank.write_text(seed_text()+'\t4\n')
            valid=subprocess.run([str(AUDITOR),'--negative-bank',str(bank),'9'],capture_output=True,text=True)
            self.assertEqual(valid.returncode,0,valid.stderr)
            for word in ['0,3','6','',','.join(['4']*65)]:
                bank.write_text(seed_text()+'\t'+word+'\n')
                result=subprocess.run([str(AUDITOR),'--negative-bank',str(bank),'9'],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                result=subprocess.run([str(PRODUCER),'--d','9','--negative-bank',str(bank),
                    '--out',str(root/('bad_'+str(len(word))))],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)

    def test_known_alias_requires_exact_state_and_bound_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            # macOS commonly has a case-insensitive filesystem; this must not
            # alias the separate KNOWN.tsv candidate journal.
            root=Path(tmp);fixture(root);bank=root/'trusted_states.tsv'
            bank.write_text('9\t'+seed_text()+'\n')
            (root/'INPUT_BANKS.tsv').write_text('--known-bank\t'+str(bank)+'\n')
            (root/'KNOWN.tsv').write_text('2\t'+seed_text()+'\n')
            path=root/'COUNTS.tsv';header,row=path.read_text().splitlines();fields=row.split('\t')
            fields[3]='2';fields[8]='1';path.write_text(header+'\n'+'\t'.join(fields)+'\n')
            (root/'DESIGN.tsv').write_text('known_hits\n1\n')
            self.assertEqual(audit_guided_native(root,AUDITOR)['known_aliases'],1)
            bank.write_text('9\t[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]\n')
            with self.assertRaisesRegex(ValueError,'false known alias'):audit_guided_native(root,AUDITOR)
            (root/'INPUT_BANKS.tsv').write_text('')
            with self.assertRaisesRegex(ValueError,'known alias input'):audit_guided_native(root,AUDITOR)

    def test_design_launch_requires_native_and_bound_validation(self):
        result=subprocess.run(['python3',str(PROJECT/'scripts/smooth/run_forest_fast.py'),'--design-search'],capture_output=True,text=True,check=True)
        plan=json.loads(result.stdout)
        self.assertTrue(plan['design_search']);self.assertIn('solved15',plan['scope'])
        self.assertEqual(plan['workers'],6);self.assertFalse(plan['gpu'])

    def test_measured_allocation_keeps_forest_and_joint_and_composite_shapes(self):
        from run_forest_fast import design_job
        self.assertEqual(sum(design_job(1,e)[1]=='joint' for e in range(32)),4)
        self.assertEqual({design_job(0,e)[0] for e in range(10)},{11,12,13,14,16})
        self.assertEqual({design_job(2,e)[0] for e in range(6)},{15,18,20,21,24,25})
        for worker in range(3,6):
            self.assertEqual({design_job(worker,e)[0] for e in range(8)},{12,15,16,18,20,21,24,25})
            self.assertTrue(all(design_job(worker,e)[1:] == ('quotient',1) for e in range(8)))

    def test_prepared_activation_is_idempotent_and_rejects_changed_inputs(self):
        from run_forest_fast import activate_prepared
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'PLAN.json').write_text(json.dumps(dict(supervisor='launchd',design_search=True)))
            (root/'LAUNCH.json').write_text('{}')
            activate_prepared(root) # No service call and no second launch.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);old=root/'predecessor';old.mkdir();(old/'FINISHED.json').write_text('{"status":"complete"}')
            (root/'STATE.json').write_text('{"status":"prepared"}')
            (root/'source').mkdir();(root/'source/test.txt').write_text('changed')
            (root/'PLAN.json').write_text(json.dumps(dict(supervisor='launchd',design_search=True,
                predecessor=str(old),source_hashes={'test.txt':'not-the-digest'})))
            with self.assertRaisesRegex(ValueError,'prepared source changed'):activate_prepared(root)

if __name__=='__main__':unittest.main()
