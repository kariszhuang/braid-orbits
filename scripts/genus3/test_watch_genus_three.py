"""Read-only dashboard tests; no search subprocess, scientific data or T7 writes."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from datetime import datetime, timezone
import json
import re
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import watch_genus_three as dashboard


class DashboardTests(unittest.TestCase):
    def test_recovery_adapter_preserves_ancestral_counts_without_reclassifying(self):
        a=self.entry('inherited',True);self.write(self.old/'NEW_COMPONENTS.json',[a])
        adapter=self.project/'adapter';receipt=self.project/'recovery.json'
        self.write(adapter/'PLAN.json',dict(recovery_adapter=True,inheritance=dict(predecessor=str(self.old))))
        self.write(adapter/'ADAPTER_PROVENANCE.json',dict(original_root=str(self.old),recovery_receipt=str(receipt)))
        self.write(receipt,dict(passed=True,finalized_attention=True,all_durable_leads_resolved=True,pending=0,admissible_predecessor=str(adapter),root=str(self.old)))
        self.write(self.root/'INHERITED_CHECKPOINT.json',dict(predecessor=str(adapter)))
        before={str(p):p.read_bytes()for p in self.project.rglob('*.json')}
        value=self.monitor.sample()
        self.assertEqual(value['new_genus3']['total'],1);self.assertEqual(value['new_genus3']['isotrivial'],1)
        self.assertEqual(value['new_genus3']['non_isotrivial'],0);self.assertEqual(value['known_genus3_total'],2)
        self.assertEqual(before,{str(p):p.read_bytes()for p in self.project.rglob('*.json')})
        self.write(receipt,dict(passed=True,root=str(self.old),pending=1))
        self.assertIn('Incomplete recovery ancestry',self.monitor.sample()['warnings'][0])

    def test_compact_worker_is_counted_without_matching_unowned_prefixes(self):
        from types import SimpleNamespace
        root=Path('/research/run')
        children=SimpleNamespace(returncode=0,stdout='1001\n1002\n1003\n')
        result=SimpleNamespace(returncode=0,stdout=(
            '999 1 0.1 python run_genus_three.py --worker /research/run\n'
            '1001 999 90.0 /research/run/smooth_genus_three --d 63\n'
            '1002 999 80.0 /research/run/smooth_genus_three_compact --d 21\n'
            '1003 999 20.0 /research/run/smooth_genus_three_compact_other --d 21\n'))
        with patch.object(dashboard.subprocess,'run',side_effect=[children,result]):
            value=dashboard.process_snapshot(root,999)
        self.assertEqual(value['search_workers'],2)
        self.assertEqual(value['worker_cpu_percent'],170.)

    def test_finite_progress_is_construction_count_not_new_examples(self):
        state=json.loads((self.root/'STATE.json').read_text())
        state.update(boundary_cursor=12345,boundary_total=504000)
        self.write(self.root/'STATE.json',state)
        value=self.monitor.sample()
        self.assertEqual(value['boundary'],dict(next=12345,total=504000))
        text=dashboard.render(value,5)
        self.assertIn('12,345 / 504,000 constructed',text)
        self.assertEqual(value['new_genus3']['total'],0)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.project = Path(self.tmp.name)
        self.root = self.project/'run'; self.root.mkdir()
        self.old = self.project/'old'; self.old.mkdir()
        self.clock = 100.
        self.processes = dict(controller_alive=True,search_workers=6,audit_workers=0,worker_cpu_percent=590.)
        self.monitor = dashboard.Monitor(self.project,processes=lambda r,p:dict(self.processes),clock=lambda:self.clock)
        self.write(self.project/'.run-control/genus-three-current.json',dict(root=str(self.root)))
        self.catalog = self.project/'published.json'
        self.write(self.catalog,dict(dataset_version=11,components=[dict(id='published',genus=3),dict(id='other',genus=4)]))
        self.write(self.root/'PLAN.json',dict(catalog=dict(path=str(self.catalog))))
        self.write(self.root/'LAUNCH.json',dict(updated=self.iso(0)))
        self.state()
        self.write(self.root/'INHERITED_CHECKPOINT.json',dict(predecessor=str(self.old)))

    def iso(self, seconds):
        return datetime.fromtimestamp(seconds,timezone.utc).isoformat()

    def write(self, path, value):
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(value))

    def state(self, status='running', draws=10000):
        self.write(self.root/'STATE.json',dict(status=status,pid=999,updated=self.iso(self.clock),
            completed_jobs=6,counts=dict(constructed=draws,scope=100,unresolved=0)))

    def entry(self, cid, status=None, genus=3, verified=True):
        folder = self.project/'certs'/cid/'component'
        folder.mkdir(parents=True,exist_ok=True)
        if verified:
            (folder.parent/'VERIFIED').write_text('admitted')
        return dict(id=cid,genus=genus,degree=21,raw_folder=str(folder),
                    original_isotrivial=status,prospective_witness_count=9)

    def test_distinct_ancestral_results_exclude_published_and_other_genera(self):
        a=self.entry('a',True); b=self.entry('b',False); u=self.entry('u',None)
        self.write(self.root/'NEW_COMPONENTS.json',[a,b,u,self.entry('published',False),self.entry('g4',False,genus=4)])
        self.write(self.old/'NEW_COMPONENTS.json',[a])
        value=self.monitor.sample()
        self.assertEqual(value['new_genus3'],dict(total=3,isotrivial=1,non_isotrivial=1,unknown=1,by_degree={21:3}))
        self.assertEqual(value['known_genus3_total'],4)
        self.assertEqual(value['average_draws_per_second'],100.)

    def test_unknown_and_prospective_witnesses_are_not_nonisotrivial(self):
        self.write(self.root/'NEW_COMPONENTS.json',[self.entry('pending-proof')])
        value=self.monitor.sample()
        self.assertEqual(value['new_genus3']['unknown'],1)
        self.assertEqual(value['new_genus3']['non_isotrivial'],0)

    def test_missing_marker_not_counted_then_retried(self):
        item=self.entry('pending',True,verified=False)
        self.write(self.root/'NEW_COMPONENTS.json',[item])
        self.assertEqual(self.monitor.sample()['new_genus3']['total'],0)
        (Path(item['raw_folder']).parent/'VERIFIED').write_text('admitted')
        self.assertEqual(self.monitor.sample()['new_genus3']['total'],1)

    def test_running_file_is_not_evidence_of_a_live_process(self):
        self.processes['controller_alive']=False
        value=self.monitor.sample()
        self.assertIn('NOT RUNNING',value['status'])
        self.assertIn('historical',dashboard.render(value,5))

    def test_stopped_average_uses_finished_duration_not_time_since_launch(self):
        self.state(status='stopped')
        self.write(self.root/'FINISHED.json',dict(seconds=50))
        self.clock=10000.
        value=self.monitor.sample()
        self.assertEqual(value['average_draws_per_second'],200.)
        self.assertIsNone(value['recent_draws_per_second'])

    def test_rolling_speed_is_completed_work(self):
        first=self.monitor.sample();self.assertIsNone(first['recent_draws_per_second'])
        self.clock=160.;self.state(draws=16000)
        self.assertEqual(self.monitor.sample()['recent_draws_per_second'],100.)
        self.clock=165.;self.state(draws=16000)
        self.assertAlmostEqual(self.monitor.sample()['recent_draws_per_second'],6000/65)

    def test_root_switch_resets_speed_history_and_follows_the_new_pointer(self):
        self.monitor.sample();self.clock=160.;self.state(draws=16000)
        self.assertIsNotNone(self.monitor.sample()['recent_draws_per_second'])
        new=self.project/'successor';new.mkdir()
        self.write(new/'PLAN.json',dict(catalog=dict(path=str(self.catalog))))
        self.write(new/'LAUNCH.json',dict(updated=self.iso(160)))
        self.write(new/'STATE.json',dict(status='running',pid=999,updated=self.iso(160),counts=dict(constructed=0)))
        self.write(new/'INHERITED_CHECKPOINT.json',dict(predecessor=str(self.root)))
        self.write(self.project/'.run-control/genus-three-current.json',dict(root=str(new)))
        value=self.monitor.sample()
        self.assertIsNone(value['recent_draws_per_second'])
        self.assertEqual(value['root'],str(new))

    def test_changed_publication_baseline_excludes_newly_published_ids(self):
        self.write(self.root/'NEW_COMPONENTS.json',[self.entry('a',True)])
        self.assertEqual(self.monitor.sample()['new_genus3']['total'],1)
        self.write(self.catalog,dict(dataset_version=12,components=[dict(id='a',genus=3)]))
        value=self.monitor.sample()
        self.assertEqual(value['published_baseline_version'],12)
        self.assertEqual(value['new_genus3']['total'],0)

    def test_monitor_does_not_write_or_modify_search_files(self):
        self.write(self.root/'NEW_COMPONENTS.json',[self.entry('a',True)])
        before={str(p):p.stat().st_mtime_ns for p in self.project.rglob('*')}
        self.monitor.sample();self.monitor.sample()
        after={str(p):p.stat().st_mtime_ns for p in self.project.rglob('*')}
        self.assertEqual(before,after)

    def test_unchanged_json_files_are_not_reparsed(self):
        self.write(self.root/'NEW_COMPONENTS.json',[self.entry('a',True)])
        self.monitor.sample();reads=self.monitor.files.reads
        with patch.object(Path,'read_text',side_effect=AssertionError('unneeded parse')):
            again=self.monitor.sample()
        self.assertEqual(self.monitor.files.reads,reads)
        self.assertEqual(again['new_genus3']['total'],1)

    def test_index_changes_and_conflicting_classifications_are_visible(self):
        a=self.entry('a',True);self.write(self.root/'NEW_COMPONENTS.json',[a])
        self.assertEqual(self.monitor.sample()['new_genus3']['total'],1)
        self.write(self.root/'NEW_COMPONENTS.json',[a,self.entry('b',False)])
        self.assertEqual(self.monitor.sample()['new_genus3']['total'],2)
        self.write(self.old/'NEW_COMPONENTS.json',[dict(a,original_isotrivial=False)])
        value=self.monitor.sample()
        self.assertEqual(value['new_genus3']['unknown'],1)
        self.assertTrue(value['warnings'])

    def test_missing_history_warns_instead_of_silently_claiming_complete_totals(self):
        self.write(self.root/'INHERITED_CHECKPOINT.json',dict(predecessor=str(self.project/'absent')))
        self.assertIn('Missing predecessor',self.monitor.sample()['warnings'][0])

    def test_process_filter_ignores_shell_command_strings_and_reused_pid(self):
        output='999 1 0.0 python run_genus_three.py --worker '+str(self.root)+'\n'
        output+='1000 999 98.0 '+str(self.root/'smooth_genus_three')+' --d 21\n'
        output+='1001 999 9.0 '+str(self.root/'smooth_negative_audit')+' --guided x\n'
        output+='1002 1 10.0 zsh -c echo '+str(self.root/'smooth_genus_three')+' --d 21\n'
        response=type('Result',(),dict(stdout=output,returncode=0))()
        children=type('Result',(),dict(stdout='1000\n1001\n',returncode=0))()
        with patch.object(dashboard.subprocess,'run',side_effect=[children,response]) as run:
            value=dashboard.process_snapshot(self.root,999)
        self.assertIn('999,1000,1001',run.call_args.args[0])
        self.assertEqual(value,dict(controller_alive=True,search_workers=1,audit_workers=1,worker_cpu_percent=107.))

    def test_color_is_only_presentation_and_plain_output_has_no_escapes(self):
        value=self.monitor.sample()
        plain=dashboard.render(value,5)
        colored=dashboard.render(value,5,color=True)
        self.assertNotIn('\033[',plain)
        self.assertIn('\033[1;36m',colored)
        self.assertEqual(re.sub(r'\x1b\[[0-9;]*m','',colored),plain)
        self.assertNotIn('Ctrl-C',plain)
        self.assertNotIn('Counts are',plain)
        self.assertIn('draws/s',plain)

    def test_nonisotrivial_success_and_unknown_have_different_colors(self):
        self.write(self.root/'NEW_COMPONENTS.json',[self.entry('found',False),self.entry('review',None)])
        text=dashboard.render(self.monitor.sample(),5,color=True)
        self.assertIn('Non-isotrivial     \033[1;32m1',text)
        self.assertIn('Unknown            \033[1;33m1',text)


if __name__=='__main__':unittest.main()
