"""Finite cursor, independent construction, and scheduler regression controls."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import run_genus_three as runner
from verify_genus_three_boundary import candidates, audit_boundary, TOTAL
from verify_genus_three_finite import profile, mul


class BoundaryTests(unittest.TestCase):
    def test_intervals_and_both_path_frames(self):
        for start, stop in ((0,5040), (322500,327600), (TOTAL-5040,TOTAL)):
            whole = list(candidates(start,stop))
            half = (start+stop)//2
            self.assertEqual(whole, list(candidates(start,half))+list(candidates(half,stop)))
            self.assertTrue(whole)
            for cursor, t in whole:
                self.assertEqual([16-len(profile(p)) for p in t], [7,7,7,15])
                self.assertEqual(profile(mul(t[0],t[1])), (8,8) if cursor<322560 else (7,9))
                self.assertEqual(mul(mul(t[0],t[1]),mul(t[2],t[3])), tuple(range(16)))
        self.assertEqual(list(candidates(TOTAL,TOTAL)), [])
        with self.assertRaises(ValueError): list(candidates(0,TOTAL+1))

    def test_finite_cursor_is_forwarded_and_completed_population_is_not_repeated(self):
        for cursor, lane, n in ((123,'involution-boundary',16),(TOTAL,'collision',12)):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);(root/'logs').mkdir();folder=root/f'jobs/job_9_00000001_{n}_{lane}';folder.mkdir(parents=True)
                (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n')
                (folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
                (root/'KNOWN_IDS.json').write_text('[]');bank=root/'bank';bank.write_text('')
                boundary=dict(passed=True,start=cursor,next=cursor+1,total=TOTAL)
                with patch.object(runner.subprocess,'run') as call, patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])), patch('verify_genus_three_boundary.audit_boundary',return_value=boundary) as audit:
                    result=runner.run_job(root,dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',boundary_cursor=cursor),9,1,bank)
                cmd=call.call_args.args[0]
                if lane=='involution-boundary':
                    self.assertEqual(cmd[cmd.index('--cursor')+1],str(cursor))
                    self.assertNotIn('--index-profile',cmd)
                    self.assertEqual(result['boundary'],boundary)
                    audit.assert_called_once_with(folder,expected_start=cursor)
                else:
                    self.assertNotIn('--cursor',cmd)
                    self.assertEqual(cmd[cmd.index('--index-profile')+1],'1,9,9,9')
                    audit.assert_not_called()

    @unittest.skipUnless(os.environ.get('G3_CYCLIC_VALIDATION'), 'requires frozen validated producer')
    def test_actual_finite_journals_and_corruption(self):
        from audit_guided_native import audit_guided_native
        base=Path(os.environ['G3_CYCLIC_VALIDATION'])
        with tempfile.TemporaryDirectory() as tmp:
            for cursor in (0,322500,TOTAL-5040):
                folder=Path(tmp)/str(cursor)
                subprocess.run([str(base/'bin/compact'),'--lane','involution-boundary','--d','16',
                    '--cursor',str(cursor),'--attempts','5040','--seconds','60',
                    '--closure-seconds','5','--cap','10000','--out',str(folder)],check=True,capture_output=True)
                self.assertTrue(audit_boundary(folder,expected_start=cursor)['passed'])
                native=audit_guided_native(folder,base/'bin/auditor',native_components=True)
                self.assertEqual(native['independent_native']['records'],audit_boundary(folder)['qualified'])
                path=folder/'BOUNDARY.tsv';before=path.read_text()
                rows=list(csv.DictReader(before.splitlines(),delimiter='\t'));rows[0]['qualified']=str(int(rows[0]['qualified'])+1)
                path.write_text('\t'.join(rows[0])+'\n'+'\t'.join(rows[0].values())+'\n')
                with self.assertRaisesRegex(ValueError,'qualified mass'):audit_boundary(folder,expected_start=cursor)
                path.write_text(before)
                with self.assertRaisesRegex(ValueError,'cursor binding'):audit_boundary(folder,expected_start=cursor+1)

    def test_controller_preserves_finite_cursor_across_graceful_handoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'banks').mkdir();(root/'banks/words_00000000.txt').write_text('')
            runner.write(root/'PLAN.json',dict(campaign_lock=str(root/'lock'),workers=1,streams=10,hours=0))
            runner.write(root/'BINDINGS.json',{});runner.write(root/'KNOWN_IDS.json',[])
            runner.write(root/'INHERITED_CHECKPOINT.json',dict(next_epochs=[1]*10,boundary_cursor=5040))
            boundaries=[]
            def job(root,plan,slot,epoch,bank):
                runner.write(root/'CHECKS.json',dict(passed=True))
                result=dict(slot=slot,epoch=epoch,folder=str(root),counts=dict(constructed=1,scope=0,transitive=0,closed=0,unresolved=0),components=[],words=[])
                if slot==9 and epoch%2:
                    start=plan['boundary_cursor'];boundaries.append(start)
                    result['boundary']=dict(passed=True,start=start,next=start+100,total=TOTAL)
                    if len(boundaries)==2:(root/'STOP').touch()
                return result
            with patch.object(runner,'storage_available',return_value=True),patch.object(runner.signal,'signal'),patch.object(runner,'run_job',side_effect=job):
                runner.controller(root)
            saved=json.loads((root/'CHECKPOINT.json').read_text())
            self.assertEqual(boundaries,[5040,5140])
            self.assertEqual(saved['boundary_cursor'],5240)
            self.assertEqual(json.loads((root/'FINISHED.json').read_text())['boundary_cursor'],5240)


if __name__=='__main__': unittest.main()
