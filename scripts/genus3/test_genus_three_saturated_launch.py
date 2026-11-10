"""C237 proof admission is separate from both native closures and labels."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_genus_three as runner


class AdmissionTests(unittest.TestCase):
    def test_exact_parent_scope_and_all_nine_bindings_required(self):
        proof=runner.PROJECT/'evidence/restriction_admission_20261024/ADMISSION.json'
        value=runner.admit_saturated_passports(proof);self.assertEqual(value['claim'],'C237')
        with tempfile.TemporaryDirectory()as tmp:
            path=Path(tmp)/'proof.json'
            for change in('wrong_genus','partner','missing_binding','changed_source','unresolved','extra_rule'):
                bad=copy.deepcopy(value)
                if change=='wrong_genus':bad['production_target_genus']=2
                elif change=='partner':bad['scope']['designated_partner']='at least k cycles'
                elif change=='missing_binding':bad['bindings'].pop()
                elif change=='changed_source':bad['bindings'][0]['sha256']='0'*64
                elif change=='unresolved':bad['checks']['parent_unresolved']=1
                else:bad['rules'].append('invented stronger filter')
                runner.write(path,bad)
                with self.assertRaises(ValueError,msg=change):runner.admit_saturated_passports(path)

    def test_proof_flag_and_separate_arithmetic_audit_reach_each_job(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);(root/'logs').mkdir();bank=root/'words';bank.write_text('');runner.write(root/'KNOWN_IDS.json',[])
            runner.write(root/'SATURATED_PASSPORT_THEOREM.json',{'test':'binding only'})
            folder=root/'jobs/job_2_00000000_12_collision';folder.mkdir(parents=True)
            (folder/'COUNTS.tsv').write_text('degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n12\tcollision\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n')
            (folder/'RESOLUTIONS.tsv').write_text('candidate\tstatus\tword\tfolder\tseed\n')
            (folder/'SATURATED_SUMMARY.tsv').write_text('enabled\tchecked_scope_occurrences\ttheorem_excluded_occurrences\tsampled_occurrences\ttrace_records\tstage\ntrue\t0\t0\t0\t0\tlate\n')
            (folder/'SATURATED_THEOREMS.tsv').write_text('method\tk\tprime\tmodulus\treason\texcluded_occurrences\tsampled_occurrences\n')
            (folder/'SATURATED_EXCLUSIONS.tsv').write_text('candidate\tstatus\tspecial\tpartner\tk\tprime\tmodulus\tmethod\treason\tviolator\tseed\n')
            plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',saturated_passport_reduction={'passed':True},saturated_passport_stage='late')
            with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                runner.run_job(root,plan,2,0,bank)
            command=call.call_args.args[0];self.assertIn('--saturated-passport-cuts',command);self.assertEqual(command[command.index('--saturated-stage')+1],'late')
            check=json.loads((folder/'CHECKS.json').read_text())['saturated_passport_exclusions']
            self.assertTrue(check['passed']);self.assertEqual(check['theorem_excluded_occurrences'],0)
            self.assertEqual(check['proof_receipt_sha256'],runner.sha(root/'SATURATED_PASSPORT_THEOREM.json'))

    def test_old_native_validation_cannot_silently_ignore_the_new_flag(self):
        current=Path(json.loads((runner.PROJECT/'.run-control/genus-three-current.json').read_text())['root'])
        while json.loads((current/'VALIDATION.json').read_text()).get('features',{}).get('saturated_passport_filters'):
            current=Path(json.loads((current/'PLAN.json').read_text())['inheritance']['predecessor'])
        command=[sys.executable,str(runner.PROJECT/'scripts/genus3/run_genus_three.py'),'--run','--workers','3',
            '--binary',str(current/'smooth_genus_three'),'--auditor',str(current/'smooth_negative_audit'),'--source',str(current/'source'),
            '--checks',str(current/'VALIDATION.json'),'--catalog',str(runner.PROJECT/'exports/latest-discoveries-20261021-v11/release-v11/all-smooth.json'),
            '--saturated-passport-proof',str(runner.PROJECT/'evidence/restriction_admission_20261024/ADMISSION.json')]
        result=subprocess.run(command,text=True,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('saturated passport filters require their completed source-bound validation',result.stderr)


if __name__=='__main__':unittest.main()
