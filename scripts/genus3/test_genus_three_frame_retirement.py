try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import copy
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import genus_three_frame_retirement as frames
import run_genus_three as runner


def write(path,value):path.write_text(json.dumps(value))


def fixture(root):
    """Tiny source-binding protocol fixture, not mathematical evidence."""
    source=root/'proof.txt';source.write_text('unit-test bytes, not a theorem')
    template=json.loads((runner.PROJECT/'evidence/restriction_admission_20261024/C251_PROPOSAL_ADMISSION.json').read_text())
    template['bindings']=[]
    for claim in('C239','C241','C243','C244','C247'):
        child=root/(claim+'_ADMISSION.json')
        write(child,dict(claim=claim,passed=True,bindings=[dict(path=str(source),sha256=frames.sha(source))]))
        template['bindings'].append(dict(path=str(child),sha256=frames.sha(child)))
    proof=root/'combined.json';write(proof,template)
    return proof,source


class FrameTests(unittest.TestCase):
    def test_audit_and_producer_same_named_plans_resolve_only_by_hash(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);proof,_=fixture(root);producer=root/'producer';producer.mkdir()
            write(root/'PLAN.json',dict(producer_root=str(producer),marker='audit plan'))
            write(producer/'PLAN.json',dict(marker='producer plan'))
            expected=frames.sha(producer/'PLAN.json')
            self.assertNotEqual(expected,frames.sha(root/'PLAN.json'))
            child=root/'C239_ADMISSION.json';record=json.loads(child.read_text())
            record['bindings'].append(dict(path=str(root/'PLAN.json'),sha256=frames.sha(root/'PLAN.json')))
            record['audit_finish']=dict(raw_bindings={'PLAN.json':expected})
            write(child,record);combined=json.loads(proof.read_text());combined['bindings'][0]['sha256']=frames.sha(child);write(proof,combined)
            checked={item['path']:item['sha256']for item in frames.admit(proof,runner.PROJECT)['recursively_checked_files']}
            self.assertEqual(checked[str((producer/'PLAN.json').resolve())],expected)
            self.assertEqual(checked[str((root/'PLAN.json').resolve())],frames.sha(root/'PLAN.json'))
            write(producer/'PLAN.json',dict(marker='changed producer bytes'))
            with self.assertRaises(ValueError):frames.admit(proof,runner.PROJECT)

    def test_recursive_byte_guards_and_exact_goal_scope(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);proof,source=fixture(root);result=frames.admit(proof,runner.PROJECT)
            self.assertEqual(len(result['recursively_checked_files']),6)
            for change in('target','missing','broader','objective'):
                value=json.loads(proof.read_text())
                if change=='target':value['production_target']['nonisotrivial']=False
                elif change=='missing':value['bindings'].pop()
                elif change=='broader':value['absolute_negative_rules'][0]['single_cycle']='arbitrary index2'
                else:value['objective_only_rule']['conclusion']='every frame is nonsmooth'
                altered=root/(change+'.json');write(altered,value)
                with self.assertRaises(ValueError):frames.admit(altered,runner.PROJECT)
            source.write_text('changed nested bytes')
            with self.assertRaises(ValueError):frames.admit(proof,runner.PROJECT)

    def test_multicycle_special_patterns_are_not_retired(self):
        n=17;a=list(range(n));a[0],a[1],a[2],a[3]=1,0,3,2
        b=list(range(n))
        for x in range(10):b[x]=(x+1)%10
        self.assertIsNone(frames.decision(n,(tuple(a),tuple(b)),[2,9,13,14]))

    def test_native_old_validation_cannot_ignore_new_frame_flag(self):
        with tempfile.TemporaryDirectory()as tmp:
            proof,_=fixture(Path(tmp))
            current=Path(json.loads((runner.PROJECT/'.run-control/genus-three-current.json').read_text())['root'])
            while json.loads((current/'VALIDATION.json').read_text()).get('features',{}).get('proposal_frame_retirement'):
                current=Path(json.loads((current/'PLAN.json').read_text())['inheritance']['predecessor'])
            command=[sys.executable,str(runner.PROJECT/'scripts/genus3/run_genus_three.py'),'--run','--workers','8',
                '--binary',str(current/'smooth_genus_three'),'--auditor',str(current/'smooth_negative_audit'),
                '--source',str(current/'source'),'--checks',str(current/'VALIDATION.json'),
                '--catalog',str(runner.PROJECT/'exports/latest-discoveries-20261021-v11/release-v11/all-smooth.json'),
                '--frame-retirement-proof',str(proof)]
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('proposal FRAME retirement requires its completed source-bound validation',result.stderr)

    def test_flag_routes_to_job_and_frame_audit_not_tuple_accounting(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);(root/'logs').mkdir();bank=root/'words';bank.write_text('')
            write(root/'KNOWN_IDS.json',[]);write(root/'FRAME_RETIREMENT_THEOREM.json',{'fixture':'only'})
            folder=root/'jobs/job_2_00000000_12_collision';folder.mkdir(parents=True)
            (folder/'COUNTS.tsv').write_text('degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n12\tcollision\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n')
            (folder/'RESOLUTIONS.tsv').write_text('candidate\tstatus\tword\tfolder\tseed\n')
            (folder/'FRAME_RETIREMENT_SUMMARY.tsv').write_text('constructor\tenabled\tchecked_new_frames\ttarget_excluded_frames\tsampled_frames\ncollision\ttrue\t0\t0\t0\n')
            (folder/'FRAME_RETIREMENT_REASONS.tsv').write_text('reason\ttarget_excluded_frames\n'+''.join(reason+'\t0\n'for reason in frames.REASONS))
            (folder/'FRAME_RETIREMENT_SAMPLES.tsv').write_text('sample\tk\tprofile\treason\tpartial_pair_identity_padded_NOT_A_CANDIDATE\n')
            (folder/'FRAME_RETIREMENT_SCOPE.txt').write_text('incomplete frames, no literal negatives')
            plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',frame_retirement_reduction={'passed':True})
            with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                result=runner.run_job(root,plan,2,0,bank)
            self.assertIn('--retire-proven-frames',call.call_args.args[0])
            audit=json.loads((folder/'CHECKS.json').read_text())['proposal_frame_retirement']
            self.assertTrue(audit['passed']);self.assertEqual(audit['target_excluded_frames'],0)
            self.assertEqual(int(result['counts']['closed']),0)


if __name__=='__main__':unittest.main()
