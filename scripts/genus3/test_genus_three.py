"""Small launcher/certificate boundaries; no launchd or real research mutation."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_genus_three as runner


class GenusThreeTests(unittest.TestCase):
    def test_complete_thirteen_receipt_retires_only_thirteen_and_rotates_open_fourteen(self):
        import copy
        proof=runner.PROJECT/'evidence/resumed_genus3_20261023/FULL_D13_CHECKS.json'
        receipt=runner.admit_complete_degree_thirteen(proof)
        plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',degree12_reduction={'passed':True},degree13_reduction={'passed':True},complete_degree13_reduction=receipt,boundary_cursor=504000)
        self.assertNotIn(13,{runner.scheduled_task(s,e,plan)[0]for s in range(10)for e in range(150)})
        for epoch,profile in enumerate(('2,8,11,11','2,9,10,11','2,10,10,10','3,9,10,10')):
            with tempfile.TemporaryDirectory()as tmp:
                root=Path(tmp);(root/'logs').mkdir();(root/'KNOWN_IDS.json').write_text('[]');bank=root/'words';bank.write_text('')
                n,lane=runner.scheduled_task(9,epoch,plan);self.assertEqual((n,lane),(14,'collision'))
                folder=root/f'jobs/job_9_{epoch:08}_14_collision';folder.mkdir(parents=True)
                (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n');(folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
                with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                    runner.run_job(root,plan,9,epoch,bank)
                command=call.call_args.args[0];self.assertEqual(command[command.index('--index-profile')+1],profile)
                self.assertEqual(command.count('--index-profile'),1)
        with tempfile.TemporaryDirectory()as tmp:
            path=Path(tmp)/'receipt.json';bad=copy.deepcopy(receipt);bad['double_transposition']['unresolved']=1;runner.write(path,bad)
            with self.assertRaisesRegex(ValueError,'incomplete'):runner.admit_complete_degree_thirteen(path)
            bad=copy.deepcopy(receipt);bad['bindings'][0]['sha256']='0'*64;runner.write(path,bad)
            with self.assertRaisesRegex(ValueError,'binding'):runner.admit_complete_degree_thirteen(path)

    def expired_predecessor(self,root):
        runner.write(root/'PLAN.json',dict(deadline_epoch=1000))
        runner.write(root/'STATE.json',dict(status='complete',active_jobs=[]))
        runner.write(root/'FINISHED.json',dict(status='complete',reason='absolute_deadline_drain',epochs=[5]*10,completed_jobs=50))
        runner.write(root/'CHECKPOINT.json',dict(deadline_epoch=1000,next_epochs=[5]*10,completed_jobs=50))

    def test_expired_deadline_default_remains_inherited(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.expired_predecessor(root)
            self.assertEqual(runner.successor_deadline(None,root,timestamp=1100),(1000,None))
            self.assertEqual(runner.successor_deadline(2000,root,timestamp=1100),(1000,None))

    def test_new_human_authority_can_reset_only_completed_expired_mission(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.expired_predecessor(root)
            original={p.name:p.read_bytes()for p in root.iterdir()}
            deadline,receipt=runner.successor_deadline(None,root,'User: dont stop yo',1100)
            self.assertIsNone(deadline);self.assertEqual(receipt['old_deadline_epoch'],1000)
            self.assertEqual(receipt['reason'],'User: dont stop yo')
            self.assertEqual(receipt['finished_sha256'],runner.sha(root/'FINISHED.json'))
            self.assertEqual(original,{p.name:p.read_bytes()for p in root.iterdir()})
            self.assertEqual(runner.successor_deadline(2000,root,'New user deadline',1100)[0],2000)
            with self.assertRaisesRegex(ValueError,'expired'):runner.successor_deadline(None,root,'User resumes',999)
            with self.assertRaisesRegex(ValueError,'reason'):runner.successor_deadline(None,root,' ',1100)

    def test_deadline_reset_rejects_live_failed_and_inconsistent_predecessors(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp)
            for status in('running','failed','attention','stopped'):
                self.expired_predecessor(root);runner.write(root/'STATE.json',dict(status=status,active_jobs=[]))
                with self.assertRaises(ValueError):runner.successor_deadline(None,root,'User resumes',1100)
            self.expired_predecessor(root);runner.write(root/'STATE.json',dict(status='complete',active_jobs=['raw']))
            with self.assertRaisesRegex(ValueError,'active'):runner.successor_deadline(None,root,'User resumes',1100)
            self.expired_predecessor(root);runner.write(root/'CHECKPOINT.json',dict(deadline_epoch=1000,next_epochs=[4]*10,completed_jobs=50))
            with self.assertRaisesRegex(ValueError,'checkpoint'):runner.successor_deadline(None,root,'User resumes',1100)
        with self.assertRaisesRegex(ValueError,'predecessor'):runner.successor_deadline(None,None,'User resumes',1100)

    def test_degree_twelve_proof_scope_and_focused_thirteen_preserve_open_shapes(self):
        import copy
        proof=runner.PROJECT/'evidence/overnight_genus3_20261023/actions/D12_CHECKS.json'
        r=runner.admit_degree_twelve(proof);self.assertTrue(r['passed'])
        with tempfile.TemporaryDirectory()as tmp:
            path=Path(tmp)/'receipt.json';bad=copy.deepcopy(r);bad['double_transposition']['unresolved']=1;runner.write(path,bad)
            with self.assertRaisesRegex(ValueError,'incomplete'):runner.admit_degree_twelve(path)
            bad=copy.deepcopy(r);bad['bindings'][0]['sha256']='0'*64;runner.write(path,bad)
            with self.assertRaisesRegex(ValueError,'binding'):runner.admit_degree_twelve(path)
        plan=dict(degree12_reduction=r,search_policy='overnight-refactor',boundary_cursor=504000)
        self.assertNotIn(12,{runner.scheduled_task(s,e,plan)[0]for s in range(10)for e in range(150)})
        for epoch,profile in((0,'2,8,10,10'),(2,'2,9,9,10'),(4,'3,9,9,9')):
            with tempfile.TemporaryDirectory()as tmp:
                root=Path(tmp);(root/'logs').mkdir();(root/'KNOWN_IDS.json').write_text('[]');bank=root/'words';bank.write_text('')
                folder=root/f'jobs/job_9_{epoch:08}_13_collision';folder.mkdir(parents=True)
                (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n');(folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
                config=dict(plan,seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro')
                with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                    runner.run_job(root,config,9,epoch,bank)
                command=call.call_args.args[0];self.assertEqual(command[command.index('--index-profile')+1],profile)
                self.assertNotIn('--exclude-single-three',command)

    def test_saturated_thirteen_proof_preserves_middle_single_three(self):
        r=runner.admit_degree_thirteen(runner.PROJECT/'evidence/overnight_genus3_20261023/actions/D13_SATURATED_CHECKS.json')
        self.assertEqual(r['remaining_special_types'],['3','2,2'])
        for slot,epoch in((9,0),(9,2),(9,4),(3,0)):
            with tempfile.TemporaryDirectory()as tmp:
                root=Path(tmp);(root/'logs').mkdir();(root/'KNOWN_IDS.json').write_text('[]');bank=root/'words';bank.write_text('')
                plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',degree12_reduction={'passed':True},degree13_reduction=r,boundary_cursor=504000)
                n,lane=runner.scheduled_task(slot,epoch,plan)
                if n!=13:continue
                folder=root/f'jobs/job_{slot}_{epoch:08}_{n}_{lane}';folder.mkdir(parents=True)
                (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n');(folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
                with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                    runner.run_job(root,plan,slot,epoch,bank)
                command=call.call_args.args[0];self.assertEqual(command[command.index('--index-profile')+1],'2,9,9,10')
                self.assertEqual(command.count('--index-profile'),1)
                self.assertNotIn('--exclude-single-three',command)
    def test_joined_refactor_alternates_and_does_not_replace_old_shapes(self):
        for epoch in range(12):
            with tempfile.TemporaryDirectory()as tmp:
                root=Path(tmp);(root/'logs').mkdir();(root/'KNOWN_IDS.json').write_text('[]');bank=root/'words';bank.write_text('')
                n,lane=runner.task(8,epoch,'overnight-refactor')
                if lane!='refactor':continue
                folder=root/f'jobs/job_8_{epoch:08}_{n}_{lane}';folder.mkdir(parents=True)
                (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n');(folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
                plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',search_policy='overnight-refactor',joined_refactor=True)
                with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                    runner.run_job(root,plan,8,epoch,bank)
                self.assertEqual('--refactor-joined'in call.call_args.args[0],(epoch//3)%2==1)
    def test_positive_profile_subprocess_obeys_remaining_absolute_deadline(self):
        import subprocess
        with patch('audit_guided_native.time.time',return_value=100),patch.object(runner.subprocess,'check_output',return_value='{"id":"case"}')as child:
            self.assertEqual(runner.bounded_profile_component(Path('/tmp/cert'),'case',113),{'id':'case'})
            self.assertEqual(child.call_args.kwargs['timeout'],13)
            self.assertIn('--profile',child.call_args.args[0])
        with patch('audit_guided_native.time.time',return_value=113),patch.object(runner.subprocess,'check_output')as child:
            with self.assertRaises(TimeoutError):runner.bounded_profile_component(Path('/tmp/cert'),'case',113)
            child.assert_not_called()
        with patch.object(runner,'profile_component',return_value={'id':'legacy'})as inline:
            self.assertEqual(runner.bounded_profile_component(Path('/tmp/cert'),'legacy'),{'id':'legacy'})
            inline.assert_called_once()
        with patch('audit_guided_native.time.time',return_value=100),patch.object(runner.subprocess,'check_output',side_effect=subprocess.TimeoutExpired('profile',13)):
            with self.assertRaises(subprocess.TimeoutExpired):runner.bounded_profile_component(Path('/tmp/cert'),'case',113)

    def test_transposition_proof_admission_is_scope_and_hash_bound(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);math=root/'proof.md';paper=root/'paper.pdf';math.write_text('reviewed proof');paper.write_bytes(b'primary source control')
            receipt=dict(schema=1,passed=True,theorem='no-positive-genus-all-smooth-transposition',source_genus='positive',transitive=True,product_one=True,production_target_genus=3,conclusion='not_all_smooth',reviewed_by_parent_and_geometry_against_primary=True,
                mathematics=dict(path=str(math),sha256=runner.sha(math)),primary_source=dict(path=str(paper),sha256=runner.sha(paper)))
            path=root/'receipt.json';runner.write(path,receipt);self.assertTrue(runner.admit_transposition_reduction(path)['passed'])
            runner.write(path,dict(receipt,production_target_genus=0))
            with self.assertRaisesRegex(ValueError,'scope'):runner.admit_transposition_reduction(path)
            runner.write(path,receipt);math.write_text('changed source')
            with self.assertRaisesRegex(ValueError,'binding'):runner.admit_transposition_reduction(path)
    def test_refactor_policy_adds_law_without_retiring_saved_streams(self):
        rows=[runner.task(8,e,'overnight-refactor')for e in range(90)]
        self.assertEqual({d for d,lane in rows if lane=='refactor'},{18,20,21,24,25})
        self.assertIn('general-fiber',{lane for _,lane in rows})
        self.assertEqual(sum(runner.task(2,e,'overnight-refactor')[1]=='refactor'for e in range(100)),50)
        self.assertEqual(sum(runner.task(2,e,'overnight-refactor')[1]=='collision'for e in range(100)),50)
        for s in range(10):
            if s not in (2,8):
                self.assertTrue(all(runner.task(s,e)==runner.task(s,e,'overnight-refactor')for e in range(100)))
        scheduled=[runner.task(s,e,'overnight-refactor')for s in range(10)for e in range(150)]
        self.assertEqual({n for n,lane in scheduled if lane=='collision'},set(range(12,64)))
        with self.assertRaisesRegex(ValueError,'policy'):runner.task(8,0,'bad')

    def test_adaptive_proposal_templates_only_enter_one_of_three_refactor_rounds(self):
        for epoch,enabled in ((0,True),(1,True),(3,False),(4,False),(6,False),(7,False),(9,True),(10,True)):
            with tempfile.TemporaryDirectory()as tmp:
                root=Path(tmp);(root/'logs').mkdir();(root/'KNOWN_IDS.json').write_text('[]');bank=root/'words';bank.write_text('')
                n,lane=runner.task(8,epoch,'overnight-refactor');folder=root/f'jobs/job_8_{epoch:08}_{n}_{lane}';folder.mkdir(parents=True)
                (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n');(folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
                plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=10,verifier='macro',search_policy='overnight-refactor',adaptive_refactor=256)
                with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                    runner.run_job(root,plan,8,epoch,bank)
                self.assertEqual('--adaptive-refactor'in call.call_args.args[0],enabled)
    def test_absolute_deadline_has_explicit_zone_and_cannot_extend(self):
        target=runner.absolute_deadline('2026-10-23T08:30:00-05:00')
        self.assertEqual(target,runner.absolute_deadline('2026-10-23T13:30:00Z'))
        with self.assertRaisesRegex(ValueError,'offset'):runner.absolute_deadline('2026-10-23T08:30:00')
        self.assertEqual(runner.inherited_deadline(target+3600,{'deadline_epoch':target}),target)
        self.assertEqual(runner.inherited_deadline(None,{'deadline_epoch':target}),target)
        self.assertEqual(runner.inherited_deadline(target-1,{'deadline_epoch':target}),target-1)
        plan=dict(deadline_epoch=target,job_seconds=45,closure_seconds=60)
        self.assertFalse(runner.deadline_dispatch_closed(plan,target-226))
        self.assertTrue(runner.deadline_dispatch_closed(plan,target-225))
        self.assertFalse(runner.deadline_dispatch_closed({},target))
        from audit_guided_native import remaining_timeout
        with patch('audit_guided_native.time.time',return_value=target-13):
            self.assertEqual(remaining_timeout(target),13)
            self.assertIsNone(remaining_timeout(None))
        with patch('audit_guided_native.time.time',return_value=target):
            with self.assertRaisesRegex(TimeoutError,'deadline'):remaining_timeout(target)

    def test_deadline_drain_starts_no_new_job_and_preserves_all_cursors(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);(root/'banks').mkdir();(root/'banks/words_00000000.txt').write_text('1,4\n')
            runner.write(root/'PLAN.json',dict(campaign_lock=str(root/'lock'),workers=8,streams=10,hours=0,deadline_epoch=1000))
            runner.write(root/'BINDINGS.json',{});runner.write(root/'KNOWN_IDS.json',[])
            runner.write(root/'INHERITED_CHECKPOINT.json',dict(next_epochs=[123]*10,boundary_cursor=504000))
            with patch.object(runner,'storage_available',return_value=True),patch.object(runner.signal,'signal'),patch.object(runner.time,'time',return_value=1000),patch.object(runner,'run_job')as job:
                runner.controller(root)
            job.assert_not_called()
            final=json.loads((root/'FINISHED.json').read_text());checkpoint=json.loads((root/'CHECKPOINT.json').read_text())
            self.assertEqual(final['reason'],'absolute_deadline_drain')
            self.assertEqual(checkpoint['next_epochs'],[123]*10)
            self.assertEqual(checkpoint['boundary_cursor'],504000)
            self.assertEqual(checkpoint['deadline_epoch'],1000)
    def test_elementary_exclusion_requires_both_bound_complete_slices(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);checks=root/'checks.json';math=root/'math.md'
            checks.write_text('{"passed":true}');math.write_text('elementary proofs')
            audits=[]
            for i,profile in enumerate(((1,9,8,8),(2,8,8,8))):
                path=root/f'audit{i}.json';runner.write(path,dict(passed=True,status='complete_negative_slice',coverage=dict(degree=11,pending=0,ordered_indices=profile)))
                audits.append(dict(path=str(path),sha256=runner.sha(path)))
            receipt=dict(schema=1,source_genus=3,identity_free=True,excluded_degrees=[11],checks=dict(path=str(checks),sha256=runner.sha(checks)),mathematics=dict(path=str(math),sha256=runner.sha(math)),finite_degree11_audits=audits)
            file=root/'sealed.json';runner.write(file,receipt)
            self.assertEqual(runner.admit_elementary_reductions(file)['excluded_degrees'],[11])
            changed=dict(receipt,finite_degree11_audits=audits[:1]);runner.write(file,changed)
            with self.assertRaisesRegex(ValueError,'missing finite'):runner.admit_elementary_reductions(file)
            runner.write(file,receipt);math.write_text('changed proof')
            with self.assertRaisesRegex(ValueError,'binding changed'):runner.admit_elementary_reductions(file)

    def test_weighted_allocation_retains_exploration_and_does_not_duplicate_busy_slots(self):
        counts=[0]*10;cursor=0;weights=[2,2,3,3,1,2,3,3,2,1]
        for _ in range(2200):
            ready,cursor=runner.available_slots(10,set(),1,cursor,weights,counts)
            self.assertEqual(len(ready),1);counts[ready[0]]+=1
        self.assertEqual(counts,[100*w for w in weights])
        ready,_=runner.available_slots(10,{0,1,2,3,4},6,0,weights,counts)
        self.assertEqual(len(ready),1);self.assertNotIn(ready[0],{0,1,2,3,4})
        with self.assertRaises(ValueError):runner.available_slots(10,set(),6,0,[0]*10,counts)
    def test_six_workers_rotate_all_ten_saved_streams(self):
        slots,cursor=runner.available_slots(10,set(),6,0)
        self.assertEqual(slots,list(range(6)));self.assertEqual(cursor,6)
        following,cursor=runner.available_slots(10,set(slots[2:]),6,cursor)
        self.assertEqual(following,[6,7]);self.assertEqual(cursor,8)
        self.assertEqual(runner.available_slots(10,set(range(6)),6,8),([],8))
        seen=set(slots+following)
        for _ in range(12):
            batch,cursor=runner.available_slots(10,set(),6,cursor)
            self.assertLessEqual(len(batch),6);self.assertEqual(len(set(batch)),len(batch))
            seen.update(batch)
        self.assertEqual(seen,set(range(10)))
        for args in [(5,set(),6,0),(10,{10},6,0),(10,set(),6,10)]:
            with self.assertRaises(ValueError):runner.available_slots(*args)

    def test_controller_limits_concurrency_without_losing_streams(self):
        import threading
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'banks').mkdir();(root/'banks/words_00000000.txt').write_text('1,4\n')
            runner.write(root/'PLAN.json',dict(campaign_lock=str(root/'lock'),workers=6,streams=10,hours=0))
            runner.write(root/'BINDINGS.json',{});runner.write(root/'KNOWN_IDS.json',[])
            old=[711,710,711,711,711,711,303,303,303,303]
            runner.write(root/'INHERITED_CHECKPOINT.json',dict(next_epochs=old))
            guard=threading.Lock(); active=0; peak=0; seen=[]
            def job(root,plan,slot,epoch,bank):
                nonlocal active,peak
                with guard:
                    active+=1;peak=max(peak,active);seen.append((slot,epoch))
                    if len(seen)>=16:(root/'STOP').touch()
                    active-=1
                return dict(slot=slot,epoch=epoch,counts=dict(constructed=1,scope=0,transitive=0,closed=0,unresolved=0),components=[],words=[])
            with patch.object(runner,'storage_available',return_value=True),patch.object(runner.signal,'signal'),patch.object(runner,'run_job',side_effect=job):
                runner.controller(root)
            checkpoint=json.loads((root/'CHECKPOINT.json').read_text())
            self.assertLessEqual(peak,6);self.assertEqual({s for s,e in seen},set(range(10)))
            self.assertEqual(len(checkpoint['next_epochs']),10)
            self.assertEqual(len(set(seen)),len(seen))
            self.assertTrue(all(e>=old[s] for s,e in seen))
            self.assertEqual(Path(checkpoint['word_bank']).read_text(),'1,4\n')

    def test_ten_worker_expansion_preserves_old_cursors_and_seed_streams(self):
        old=[117,98,103,101,70,89]
        self.assertEqual(runner.expanded_epochs(old,10),old+[0]*4)
        self.assertEqual(old,[117,98,103,101,70,89])
        with self.assertRaises(ValueError):runner.expanded_epochs(old,4)
        with self.assertRaises(ValueError):runner.expanded_epochs([-1],10)
        for slot in range(6):
            for epoch in [0,1,100,2**30]:
                expected=(20261021+0x9E3779B97F4A7C15*(1+slot+6*epoch))%2**64
                self.assertEqual(runner.job_seed(20261021,slot,epoch),expected)
        seeds={runner.job_seed(20261021,slot,epoch) for slot in range(10) for epoch in range(1000)}
        self.assertEqual(len(seeds),10000)
        self.assertTrue(all(12<=runner.task(slot,e)[0]<=63 for slot in range(10) for e in range(140)))
    def test_schedule_has_independent_and_structured_lanes_and_all_degrees(self):
        scheduled = [runner.task(s, e) for s in range(6) for e in range(136)]
        self.assertEqual({lane for _, lane in scheduled}, {'collision','legacy-types','elliptic','rational','mutation','full-fiber','partial-fiber','general-fiber','four-quotient'})
        self.assertEqual({n for n, lane in scheduled if lane == 'collision'}, set(range(12,64)))
        self.assertTrue(all(n % r == 0 for n, lane in scheduled if lane == 'elliptic'
                            for r in [next(r for r in range(3,n//3+1) if n%r==0)]))
        self.assertTrue({36,48,63}<={n for n,lane in scheduled if lane=='full-fiber'})
        self.assertFalse({15,18,20}&{n for n,lane in scheduled if lane=='full-fiber'})
        self.assertEqual({runner.task(9,e)for e in range(50)},{(12,'collision'),(16,'involution-boundary')})

    def test_focused_stream_passes_scope_to_actual_producer(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);(root/'logs').mkdir();folder=root/'jobs/job_9_00000000_12_collision';folder.mkdir(parents=True)
            (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n');(folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
            (root/'KNOWN_IDS.json').write_text('[]');bank=root/'words';bank.write_text('')
            with patch.object(runner.subprocess,'run')as call,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                runner.run_job(root,dict(seed=1,job_seconds=1,closure_seconds=1,cap=100,verifier='macro'),9,0,bank)
            cmd=call.call_args.args[0];self.assertEqual(cmd[cmd.index('--index-profile')+1],'1,9,9,9')
            self.assertIn('--literal-trie',cmd)

    def test_zero_new_successor_retains_ancestral_known_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'KNOWN_IDS.json').write_text('["ancestral", "catalog"]')
            self.assertEqual(runner.inherited_known(root), ['ancestral','catalog'])
            runner.write(root/'BINDINGS.json', {'KNOWN_IDS.json': runner.sha(root/'KNOWN_IDS.json')})
            self.assertEqual(runner.inherited_known(root), ['ancestral','catalog'])
            (root/'KNOWN_IDS.json').write_text('[]')
            with self.assertRaisesRegex(ValueError,'binding changed'):
                runner.inherited_known(root)

    def test_words_are_bounded_ranked_and_do_not_claim_smoothness(self):
        words = runner.learned_words([[1,4,2]]*4 + [[0,1]]*2 + [[9],[],[0]*25])
        self.assertEqual(words, [(1,4,2),(0,1)])
        self.assertEqual(len(runner.learned_words([[x//6,x%6] for x in range(36)], 5)), 5)

    def test_empty_jobs_cannot_erase_inherited_word_bank(self):
        inherited=[[1,4,2],[0,1]]
        self.assertEqual(runner.next_word_bank([],inherited),[(1,4,2),(0,1)])
        self.assertEqual(runner.next_word_bank([[0,1]]*3,inherited),[(0,1),(1,4,2)])

    def test_component_identity_is_same_compact_json_convention(self):
        from export_research_catalog import identifier
        self.assertEqual(runner.identity(4, ((1,2,3,0),)*3), identifier(4, ((1,2,3,0),)*3))

    def test_atomic_state_has_no_half_pointer(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'STATE.json'; runner.write(p, {'status':'running'})
            self.assertEqual(json.loads(p.read_text()), {'status':'running'})
            runner.write(p, {'status':'stopped'})
            self.assertFalse(p.with_suffix('.json.tmp').exists())

    def test_storage_guard_stops_before_any_producer_and_retains_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'banks').mkdir();(root/'banks/words_00000000.txt').write_text('')
            runner.write(root/'PLAN.json',dict(campaign_lock=str(root/'lock'),workers=2,hours=0))
            runner.write(root/'BINDINGS.json',{});runner.write(root/'KNOWN_IDS.json',[])
            with patch.object(runner,'storage_available',return_value=False),patch.object(runner.signal,'signal'),patch.object(runner,'run_job') as run:
                runner.controller(root)
            run.assert_not_called()
            self.assertEqual(json.loads((root/'FINISHED.json').read_text())['reason'],'disk_guard_10GiB')
            self.assertEqual(json.loads((root/'CHECKPOINT.json').read_text())['next_epochs'],[0,0])

    def test_only_positive_genus_three_seeds_enter_mutation_input(self):
        entry=dict(degree=4,genus=3,tuple=['()']*4,block_quotients=[dict(blocks=[[1,2],[3,4]])])
        text=runner.mutation_seed_text(dict(components=[entry,dict(genus=2)]))
        self.assertEqual(text,'4\t[(), (), (), ()]\t0,0,1,1\n')

    def test_job_only_admits_after_full_independent_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'logs').mkdir(); folder=root/'jobs/job_2_00000000_12_collision'
            folder.mkdir(parents=True); bank=root/'words.txt'; bank.write_text('')
            (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n')
            (folder/'RESOLUTIONS.tsv').write_text('status\tword\nnonsmooth\t1,4\n')
            (root/'KNOWN_IDS.json').write_text('[]')
            plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=100,verifier='macro')
            checks=dict(closed_identities=[])
            with patch.object(runner.subprocess,'run'), patch.object(runner,'audit_guided_native',return_value=checks):
                result=runner.run_job(root,plan,2,0,bank)
            self.assertEqual(result['words'], [[1,4]])
            self.assertTrue((folder/'VERIFIED').exists())
            (folder/'VERIFIED').unlink()
            with patch.object(runner.subprocess,'run'), patch.object(runner,'audit_guided_native',side_effect=ValueError('bad edge')):
                with self.assertRaisesRegex(ValueError,'bad edge'):runner.run_job(root,plan,2,0,bank)
            self.assertFalse((folder/'VERIFIED').exists())

    def test_compact_variant_is_used_only_with_an_explicit_validated_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'logs').mkdir();bank=root/'words';bank.write_text('')
            (root/'KNOWN_IDS.json').write_text('[]')
            folder=root/'jobs/job_2_00000000_12_collision';folder.mkdir(parents=True)
            (folder/'COUNTS.tsv').write_text('closed\tunresolved\n0\t0\n')
            (folder/'RESOLUTIONS.tsv').write_text('status\tword\n')
            plan=dict(seed=1,job_seconds=1,closure_seconds=1,cap=100,verifier='macro',compact_binary_sha256='bound')
            with patch.object(runner.subprocess,'run')as launch,patch.object(runner,'audit_guided_native',return_value=dict(closed_identities=[])):
                runner.run_job(root,plan,2,0,bank)
            self.assertEqual(launch.call_args.args[0][0],str(root/'smooth_genus_three_compact'))


if __name__=='__main__': unittest.main()
