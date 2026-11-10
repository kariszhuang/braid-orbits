"""Small stdlib-only checks for the long-running research scheduler."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import datetime as dt
from pathlib import Path
import tempfile
import json
import itertools
from unittest.mock import patch
from collections import Counter

from run_campaign import candidates_from_files, certificate_header, report, research_jobs, seed_profile, class_size, exact_slices, genus_zero_profiles, exact_priority, CORE_TIMEZONE, scheduled_threads, seconds_until_thread_change, census_progress, format_status_bar
from run_campaign import degree11_involution_slices, braid_fixed_jobs, rigid_glue_jobs, certificate_identity, block_degree, s18_family_jobs, mechanism_jobs, structure_jobs, pair_coverage, counterexample_pair_jobs, small_symmetric_slices, exact_scan_complete
from run_campaign import structured_slices, structured_followup, exact_coverage
from run_campaign import discovery_plan, discovery_command, discovery_pair_profile, rigid_coverage
from run_campaign import deduplicate_representations
from run_campaign import completed_rigid_profiles
from run_campaign import exploration_passport, max_power_index, inherit_discovery, discovery_coverage, weak_symmetry_rows
from run_campaign import hybrid_discovery_rows, discovery_worker_cost, discovery_quotas, hybrid_discovery_summary, discovery_prefixes, unfinished_hybrid_rows
from run_campaign import defect_discovery_row, defect_discovery_summary
from run_campaign import DiscoveryProofCache, coupled_profiles
from run_campaign import small_passport_plan, stability_passport_plan, useful_passport_plan, inherited_passport_exclusions, PassportQueue, reservation_allows, small_degree_scope_audit, covered_construction_targets

SEED = '[(1 7 5)(2 9 6)(3 8 4), (1 2 3)(4 6 5)(7 8 9), (1 8 4)(2 3 7)(5 6 9), (1 5 6)(2 8 7)(3 9 4)]'


def check():
    full=discovery_quotas(10,True,2,stability_mission=True)
    assert full==dict(exact=2,exploration=1,small=1,mechanism=2,hybrid=4)
    assert sum(full.values())==10
    for workers in (3,4,8):
        try:discovery_quotas(workers,True,2,stability_mission=True)
        except ValueError:pass
        else:raise AssertionError('Unsafe hybrid stability budget accepted')
    population=stability_passport_plan({})
    useful_population=useful_passport_plan({})
    assert {r['genus'] for r in useful_population['rows'] if (r['group'],r['d'])==('S',10)}==set(range(6))
    assert {r['name'] for r in population['rows']}<={r['name'] for r in useful_population['rows']}
    audit=small_degree_scope_audit({})
    assert not covered_construction_targets(audit)
    complete=dict(audit,strata=[dict(r,missing_passports=0) for r in audit['strata']])
    assert covered_construction_targets(complete)=={('A',9),('S',9),('A',10)}
    complete['strata'][0]['missing_passports']=1
    assert ('A',9) not in covered_construction_targets(complete)
    counts=Counter((r['group'],r['d'],r['genus']) for r in population['rows'])
    assert counts[('S',9,2)]==1365 and counts[('S',9,3)]==588
    assert [counts[('A',10,g)] for g in range(4)]==[765,1005,825,495]
    assert set((g,d) for g,d,_ in counts)=={('S',9),('A',10),('S',10),('A',11),('S',11)}
    assert len({r['name'] for r in population['rows']})==len(population['rows'])
    controls=weak_symmetry_rows(0,Path('known'),include_joint=True)
    joint=[r for r in controls if r['family']=='joint-pure']
    assert len(joint)==36 and {r['joint_repair_rounds'] for r in joint}=={0,8}
    assert all(r['pure_power']==1 and r['pure_braid_start'] for r in joint)
    command=discovery_command(joint[0],Path('binary'),Path('attempt'),5,Path('known'),Path('source'),None,123)
    assert '--joint-repair-rounds' in command and '--no-braid-fixed' not in command
    useful=weak_symmetry_rows(0,Path('known'),include_joint=True,include_useful=True)
    upgraded=[r for r in useful if r['family']=='joint-pure']
    assert upgraded and all(r['joint_admission_aware'] and r['joint_repair_rounds']==12 for r in upgraded)
    assert {r['d'] for r in upgraded}==set(range(9,18))|{25}
    hidden=next(r for r in upgraded if r.get('joint_hidden_fiber'))
    command=discovery_command(hidden,Path('binary'),Path('attempt'),5,Path('known'),Path('source'),None,123)
    assert all(flag in command for flag in ('--joint-admission-aware','--joint-pure-basis','--joint-hidden-fiber','--construction-witness-limit'))
    assert sum(discovery_quotas(7,True,1,stability_mission=True,useful_mission=True).values())==7
    assert discovery_quotas(7,True,1,stability_mission=True,useful_mission=True)==dict(exact=1,exploration=1,small=1,mechanism=2,hybrid=2)
    audit=small_degree_scope_audit({})
    assert not audit['necessary_filter_scope_complete'] and not audit['whole_degree_nonexistence_claim']
    assert {(r['group'],r['degree'],r['genus']) for r in audit['strata']}=={('A',9,g) for g in range(5)}|{('S',9,g) for g in range(4)}|{('A',10,g) for g in range(4)}
    with tempfile.TemporaryDirectory(prefix='hurwitz-exclusion-') as temporary:
        root=Path(temporary);folder=root/'jobs'/'passport';attempt=folder/'attempt_000';attempt.mkdir(parents=True)
        row=dict(name='passport',kind='passport',group='S',d=9,passport='3/3/3/3')
        (attempt/'VERIFIED').touch()
        (attempt/'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='test',passport=row['passport'],enumeration='centralizer-v1')))
        summary=attempt/'SUMMARY.txt'
        with patch('run_campaign.finite_passport_exclusions',return_value={}),patch('run_campaign.discovery_chain',return_value=[(root,dict(binary_sha256='test'),dict(lanes={'exact':[row]}))]):
            summary.write_text('scan start = 0\nscan next = 20\nscan total candidate pairs = 20\nscan enumeration = centralizer-v1\nhard unsolved candidates = 1\n')
            assert inherited_passport_exclusions(root)=={}
            summary.write_text(summary.read_text().replace('candidates = 1','candidates = 0'))
            assert ('S',9,('3','3','3','3')) in inherited_passport_exclusions(root)
    quotas=discovery_quotas(3,False,passport_mission=True)
    assert quotas==dict(exact=2,discovery=1)
    assert not reservation_allows('discovery',dict(exact=1,discovery=1),{'exact','discovery'},quotas)
    assert reservation_allows('exact',dict(exact=1,discovery=1),{'exact','discovery'},quotas)
    assert not reservation_allows('discovery',dict(exact=1,discovery=1),{'discovery'},quotas)
    assert reservation_allows('discovery',dict(exact=0,discovery=2),{'discovery'},quotas)
    for workers,hybrid in [(2,False),(5,True)]:
        try:discovery_quotas(workers,hybrid,passport_mission=True)
        except ValueError:pass
        else:raise AssertionError('invalid passport resource plan accepted')
    with tempfile.TemporaryDirectory(prefix='hurwitz-proof-cache-') as temporary:
        root=Path(temporary);(root/'inputs').mkdir();cache=DiscoveryProofCache(root)
        original=cache.snapshot(9);attempt=root/'jobs/test/attempt_000';attempt.mkdir(parents=True)
        (attempt/'POSITIVE_WITNESSES.tsv').write_text(f'degree\tmoves\tseed\n9\t1,0\t{SEED}\n')
        (attempt/'SUMMARY.txt').write_text('representation key = example\nentry order divides = 0\nsearch complete = false\n')
        cache.register(attempt);assert cache.snapshot(9)==original and not cache.representations
        (attempt/'VERIFIED').touch();cache.register(attempt)
        current=cache.snapshot(9);assert current!=original and SEED in Path(current).read_text()
        assert SEED not in Path(original).read_text()
        (attempt/'SUMMARY.txt').write_text('representation key = example\nentry order divides = 0\nsearch complete = true\nhard unsolved candidates = 0\n')
        cache.register(attempt);assert ('example',0) in cache.representations
        restored=DiscoveryProofCache(root);assert restored.words==cache.words and restored.representations==cache.representations
    stability=weak_symmetry_rows(0,Path('evidence/certificates'),include_extensions=True,include_stability=True)
    coupled=[r for r in stability if r['family']=='coupled-rigid']
    assert {9,10,11}<={r['d'] for r in coupled}
    assert len({r['name'] for r in stability})==len(stability)
    for r in coupled:
        assert sum(int(b.split(':')[0]) for b in r['representation_blocks'])==r['d']
        cmd=discovery_command(r,'binary','out',1,'known','source',seed=19)
        assert cmd[1]=='rigid-glue' and '--certificate' not in cmd and '--representation-seed' in cmd
    broad=[r for r in stability if r.get('comparison_block')]
    assert len(broad)==48
    for block in {r['comparison_block'] for r in broad}:
        arms=[r for r in broad if r['comparison_block']==block]
        assert len(arms)==4 and len({r['seed_key'] for r in arms})==1
        assert {r['local_frame_cap'] for r in arms}=={512,2048}
    assert all(r.get('retire_boundary_plans') for r in stability if r['family'] in ('independent-glue','double-coset-glue'))
    assert not any(r.get('retire_boundary_plans') for r in stability if r['family']=='coordinated-glue')
    assert {(defect_discovery_row(i)['group'],defect_discovery_row(i)['d']) for i in range(22)}=={(g,d) for d in (*range(9,18),18,25) for g in ('A','S')}
    row=defect_discovery_row(0)
    cmd=discovery_command(row,'binary','out',30,'known','source',seed=7)
    assert cmd[1]=='defect-mission' and cmd[cmd.index('--threads')+1]=='1' and '--weak-only' not in cmd
    with tempfile.TemporaryDirectory(prefix='hurwitz-defect-batch-') as temporary:
        folder=Path(temporary);row=dict(row,draws=1)
        draw=folder/'draw_000000';draw.mkdir();(draw/'DONE').write_text('walk_positive_defect')
        (draw/'SUMMARY.tsv').write_text('group\tdegree\tgenus\tidentity_entry\tstatus\nA\t9\t1\tfalse\twalk_positive_defect\n')
        defect_discovery_summary(folder,row,0.1)
        assert 'accepted generating seeds = 0' in (folder/'SUMMARY.txt').read_text()
        assert 'defect eligible starts = 1' in (folder/'SUMMARY.txt').read_text()
        (draw/'SUMMARY.tsv').write_text('group\tdegree\tgenus\tidentity_entry\tstatus\nA\t9\t1\tfalse\tclosed\n')
        try:defect_discovery_summary(folder,row,0.1)
        except RuntimeError:pass
        else:raise AssertionError('uncertified all-zero claim accepted')
        (draw/'SUMMARY.tsv').write_text('group\tdegree\tgenus\tidentity_entry\tstatus\nA\t9\t1\tfalse\tunknown\n')
        try:defect_discovery_summary(folder,row,0.1)
        except RuntimeError:pass
        else:raise AssertionError('unknown defect completion accepted')
    from compute_spectra import component_files
    import compute_spectra as spectra
    with tempfile.TemporaryDirectory(prefix='hurwitz-s8-catalog-') as temporary:
        root=Path(temporary)
        with patch.object(spectra,'PROJECT',root),patch.object(spectra,'COMPONENTS',root/'legacy'):
            assert component_files('S',8,8)==([],['S8'])
            for i in range(16):
                folder=root/f'runs/exact/S8_shards/shard_{i}';folder.mkdir(parents=True)
                (folder/'generating_components.csv').touch()
                (folder/'SUMMARY.txt').write_text(f'group = S_8\nsignature shard = {i}/16\n')
            paths,missing=component_files('S',8,8)
            assert len(paths)==16 and not missing
            (root/'runs/exact/S8_shards/shard_0/SUMMARY.txt').unlink()
            assert component_files('S',8,8)==([],['S8'])
    from run_campaign import candidate_paths, recover_interrupted_candidates
    with tempfile.TemporaryDirectory(prefix='hurwitz-live-survivor-') as temporary:
        root=Path(temporary);(root/'inputs').mkdir()
        active=root/'jobs/live/attempt_000';active.mkdir(parents=True)
        survivor=active/'survivors.txt';survivor.write_text('group=')
        # A worker may be paused halfway through any record, not just a line.
        assert candidate_paths(root)==[]
        assert candidates_from_files(candidate_paths(root),weak_only=True)=={}
        assert candidate_paths(root,verified=False)==[survivor]
        survivor.write_text(f'group=S_18 seed={SEED}\n')
        assert recover_interrupted_candidates(root)=={'S_18':{SEED}}
        recovered=root/'inputs/RECOVERED_PENDING.txt'
        assert recovered in candidate_paths(root)
        assert candidates_from_files(candidate_paths(root),weak_only=True)=={'S_18':{SEED}}
        (active/'VERIFIED').touch()
        (active/'RESULTS.tsv').write_text(f'status\tseed\nweak_found\t{SEED}\n')
        assert candidates_from_files(candidate_paths(root),weak_only=True)=={}
        survivor.write_text('group=')
        try:candidates_from_files(candidate_paths(root),weak_only=True)
        except ValueError:pass
        else:raise AssertionError('malformed closed candidate silently accepted')
    with tempfile.TemporaryDirectory(prefix='hurwitz-report-sidecar-') as temporary:
        folder=Path(temporary);profiles=folder/'known_profiles';profiles.mkdir()
        (profiles/'._CERTIFIED_example.tsv').write_bytes(bytes.fromhex('00051607')+b'\xb0metadata')
        (profiles/'CERTIFIED_example.tsv').write_text('group\tdegree\tpassport\nS\t18\t2/2/2/2\n')
        report(folder)
        assert (folder/'REPORT.txt').exists()
        from run_campaign import run_event
        run_event(folder,'job_started',seed='18446744073709551615')
        run_event(folder,'job_verified',coverage=dict(traversal_exhausted=False))
        events=[json.loads(l) for l in (folder/'EVENTS.jsonl').read_text().splitlines()]
        assert [e['event'] for e in events]==['job_started','job_verified']
        assert events[0]['seed']=='18446744073709551615' and events[0]['time_utc'].endswith('+00:00')
    # Saved launch metadata keeps its old absolute path after project relocation.
    # Exercise command generation with no compatibility symlink.
    import run_campaign as runner
    with tempfile.TemporaryDirectory(prefix='hurwitz-data-') as temporary:
        data=Path(temporary)/'data';link=Path(temporary)/'runs';link.symlink_to(data/'Archive/runs')
        with patch.object(runner,'RESULTS',data):
            runner.require_results_volume(runner.PROJECT/'internal-fixture')
            try:runner.require_results_volume(data/'campaign_new')
            except RuntimeError:pass
            else:raise AssertionError('missing data/ folder silently accepted results')
            try:runner.require_results_volume(link/'campaigns/old')
            except RuntimeError:pass
            else:raise AssertionError('compatibility link bypassed the data/ guard')
            data.mkdir();runner.require_results_volume(data/'campaign_new')
    with tempfile.TemporaryDirectory(prefix='hurwitz-relocation-') as temporary:
        moved=Path(temporary);executable=moved/'runs/old/bin/pinned'
        executable.parent.mkdir(parents=True);executable.touch()
        legacy=Path('/missing/hurwitz-project')
        row=dict(kind='sample',group='S',d=9)
        with patch.object(runner,'PROJECT',moved),patch.object(runner,'LEGACY_PROJECT',legacy,create=True):
            cmd=discovery_command(row,str(legacy/'runs/old/bin/pinned'),'attempt',30,'known','source')
            assert cmd[0]==str(executable) and Path(cmd[0]).is_file()
            assert runner.coverage_prefix_list([str(legacy/'runs/old')])==[str(moved/'runs/old')]
            assert runner.relocated_path('')=='' and runner.relocated_path('/unrelated/path')=='/unrelated/path'
            try:runner.coverage_prefix_list([str(legacy/'runs/old'),str(moved/'runs/old')])
            except RuntimeError:pass
            else:raise AssertionError('relocation concealed duplicate coverage')
            # Never-started jobs legitimately have no prefix directory. Missing
            # recorded work, however, must not silently restart at zero.
            assert runner.relocated_path(legacy/'missing')==str(moved/'missing')
            inherited=dict(kind='pair',pair_blocks=[2],pair_paths=[],even_paths=[],third_order=3,group='S',
                           coverage_prefixes=[str(legacy/'missing')],prior_attempts=1)
            try:discovery_coverage(inherited,moved/'current','sha')
            except RuntimeError:pass
            else:raise AssertionError('missing verified work reset coverage')
            inherited['prior_attempts']=0
            assert discovery_coverage(inherited,moved/'current','sha')['verified_attempts']==0
    gpu_rows=hybrid_discovery_rows(0,262144,3)
    assert len(gpu_rows)==21 and len({r['name'] for r in gpu_rows})==21
    assert {r['d'] for r in gpu_rows}==set(range(9,19))|{25}
    assert not {r['name'] for r in gpu_rows}&{r['name'] for r in hybrid_discovery_rows(1,262144,3)}
    assert all(discovery_worker_cost(r)==2 for r in gpu_rows)
    for workers in range(5,11):
        quotas=discovery_quotas(workers,True)
        assert sum(quotas.values())==workers and quotas['hybrid']==2 and all(n>=1 for n in quotas.values())
        assert sum(discovery_worker_cost(r) for r in [gpu_rows[0]]+[dict(kind='pair')]*(workers-2))==workers
    for workers in range(7,11):
        quotas=discovery_quotas(workers,True,2)
        assert sum(quotas.values())==workers and quotas['hybrid']==4
        assert sum(discovery_worker_cost(r) for r in gpu_rows[:2]+[dict(kind='pair')]*(workers-4))==workers
    cmd=discovery_command(dict(gpu_rows[0],runtime='runtime'),'binary','attempt',30,'known','source',seed=123)
    assert cmd[0]=='runtime/hybrid_pilot' and cmd[5]=='gpu' and cmd[-1]=='2' and cmd[10]=='2'
    compact=discovery_command(dict(gpu_rows[0],runtime='runtime',compact_inputs=True),'binary','attempt',30,'known','source',seed=123)
    assert compact==cmd+['1']
    with tempfile.TemporaryDirectory(prefix='hurwitz-hybrid-report-') as temporary:
        folder=Path(temporary)
        # An old epoch can still have a job in flight when the next epoch starts.
        # Resume recovers it instead of skipping directly to the newer epoch.
        recovery=folder/'jobs'/gpu_rows[0]['name'];saved=recovery/'attempt_000';saved.mkdir(parents=True)
        (saved/'JOB_INPUT.json').write_text(json.dumps(gpu_rows[0]))
        assert unfinished_hybrid_rows(folder)==[gpu_rows[0]]
        (recovery/'DONE').touch();assert unfinished_hybrid_rows(folder)==[]
        summary=dict(backend='gpu',screen_mode=2,shadow_mismatches=0,degree=18,family=2,pipeline_seconds=1.0,shadow_audit_seconds=0.5,rows=[dict(stats=[0,0,0,0,0,7,3,2,0,2,0,0])])
        (folder/'SUMMARY.json').write_text(json.dumps(summary));hybrid_discovery_summary(folder)
        assert 'hard unsolved candidates = 2' in (folder/'SUMMARY.txt').read_text()
        assert 'search complete = false' in (folder/'SUMMARY.txt').read_text()
        # Pending inputs survive interruption; completed witness rows suppress
        # replay, while capped inputs remain eligible for deeper exact closure.
        (folder/'survivors.txt').write_text(f'group=S_9 seed={SEED}\ngroup=S_9 seed=pending\n')
        (folder/'RESULTS.tsv').write_text(f'status\tseed\nweak_found\t{SEED}\n')
        assert candidates_from_files([folder/'survivors.txt'],weak_only=True)=={'S_9':{'pending'}}
        summary['shadow_mismatches']=1;(folder/'SUMMARY.json').write_text(json.dumps(summary))
        try:hybrid_discovery_summary(folder)
        except RuntimeError:pass
        else:raise AssertionError('unaudited hybrid job accepted')
    from hybrid_pilot import schedule
    hybrid=schedule([6,8,10],[65536,262144],2)
    assert hybrid==schedule([6,8,10],[65536,262144],2) and len(hybrid)==216
    for i in range(0,len(hybrid),2):
        a,b=hybrid[i:i+2]
        assert {a['backend'],b['backend']}=={'cpu','gpu'}
        assert all(a[k]==b[k] for k in ('degree','family','workers','batch','repeat','seed','block'))
        producers=2 if a['workers']>=8 else 1
        assert producers+(a['workers']-producers)==a['workers']
    frontier=discovery_plan({},frontier=True)
    small=frontier['lanes']['small']
    assert {r['d'] for r in small}=={9,10,11}
    assert {r['third_order'] for r in small}==set(range(2,7))
    assert small[0]['group']=='S'
    assert all(r['group'] in ('S','A') for r in small)
    assert any(r['third_order']==5 for r in small)
    def power_partitions(n,q,minimum=1):
        if not n:yield ()
        for k in range(minimum,n+1):
            if q%k==0:
                for rest in power_partitions(n-k,q,k):yield (k,)+rest
    for d in range(4,12):
        for q in range(2,7):
            assert max_power_index(d,q)==max(d-len(p) for p in power_partitions(d,q))
    for d in (9,10,11):
        for group in ('S','A'):
            for genus in range(d-4):
                passport=exploration_passport(d,group,genus,20260929)
                # S9 genus four is impossible: all four maximal indices six
                # would be even, and hence generate a subgroup of A9.
                if passport is None:
                    assert (group=='S' and d%2==1 and genus==d-5) or (group=='A' and d%2==0 and genus>d-7)
                    continue
                indices=[sum(int(n)-1 for n in p.split('.')) for p in passport.split('/')]
                assert sum(indices)==2*d-2+2*genus and all(i<=d-3 for i in indices)
                assert all(i%2==0 for i in indices) if group=='A' else any(i%2 for i in indices)
    with tempfile.TemporaryDirectory(prefix='hurwitz-prefix-inheritance-') as temporary:
        parent=Path(temporary)/'old';parent.mkdir();new=Path(temporary)/'new';new.mkdir()
        row=dict(small[0]);row['name']='partial';profile=discovery_pair_profile(row)
        old_folder=parent/'jobs'/'partial';old=old_folder/'attempt_000';old.mkdir(parents=True)
        new_folder=new/'partial';part=new_folder/'attempt_000';part.mkdir(parents=True)
        for folder,start,end in [(old,0,5),(part,5,10)]:
            (folder/'VERIFIED').touch();(folder/'PAIR_INPUT.json').write_text(json.dumps(dict(binary_sha256='old-hash',**profile)))
            (folder/'SUMMARY.txt').write_text(f'scan start h = 0\nscan start third = {start}\nscan next h = 0\nscan next third = {end}\nscan range resolved = true\ntraversal exhausted = false\n')
        c=pair_coverage(new_folder,'old-hash',profile,old_folder)
        assert c['next_third']==10 and c['verified_attempts']==2
        (parent/'STATE.json').write_text(json.dumps(dict(plan='discovery',status='stopped',binary_sha256='old-hash')))
        (parent/'inputs').mkdir();(parent/'inputs/DISCOVERY_PLAN.json').write_text(json.dumps(dict(lanes={'small':[row]})))
        schedule=inherit_discovery(dict(lanes={'small':[dict(row)]},excluded=[]),parent)
        inherited=schedule['lanes']['small'][0]
        assert inherited['enumeration_sha256']=='old-hash' and inherited['prior_attempts']==1
        assert discovery_coverage(inherited,new_folder,'new-hash')['next_third']==10
        middle=Path(temporary)/'middle';middle.mkdir();(middle/'inputs').mkdir()
        (middle/'STATE.json').write_text(json.dumps(dict(plan='discovery-frontier',status='stopped',binary_sha256='middle-hash')))
        (middle/'inputs/DISCOVERY_PLAN.json').write_text(json.dumps(dict(lanes={'small':[inherited]})))
        middle_job=middle/'jobs'/'partial'/'attempt_000';middle_job.mkdir(parents=True)
        for name in ('VERIFIED','PAIR_INPUT.json','SUMMARY.txt'):
            (middle_job/name).write_bytes((part/name).read_bytes())
        chained=inherit_discovery(dict(lanes={'small':[dict(row)]},excluded=[]),middle)['lanes']['small'][0]
        assert discovery_prefixes(chained)==[str(old_folder),str(middle_job.parent)]
        assert chained['enumeration_sha256']=='old-hash' and chained['enumeration_binary']==inherited['enumeration_binary']
        assert chained['prior_attempts']==2
        assert discovery_coverage(chained,Path(temporary)/'empty-job','newest-hash')['next_third']==10
        try:discovery_coverage(dict(chained,coverage_prefixes=[str(old_folder)]*2),new_folder,'newest-hash')
        except RuntimeError:pass
        else:raise AssertionError('duplicate prefix was accepted')
        bad=part/'PAIR_INPUT.json';bad.write_text(json.dumps(dict(binary_sha256='new-hash',**profile)))
        try:discovery_coverage(inherited,new_folder,'new-hash')
        except RuntimeError:pass
        else:raise AssertionError('changed binary was allowed to inherit a prefix')
        summary=old/'SUMMARY.txt';summary.write_text(summary.read_text().replace('traversal exhausted = false','traversal exhausted = true'))
        complete=inherit_discovery(dict(lanes={'small':[dict(row)]},excluded=[]),parent)
        assert not complete['lanes']['small'] and len(complete['excluded'])==1
        # A completed row is absent from middle's lanes; it must stay absent
        # when a fresh third-generation plan is built.
        (middle/'inputs/DISCOVERY_PLAN.json').write_text(json.dumps(dict(lanes={'small':[]},predecessor=str(parent))))
        grandchild=inherit_discovery(dict(lanes={'small':[dict(row)]},excluded=[]),middle)
        assert not grandchild['lanes']['small'] and len(grandchild['excluded'])==1
        assert grandchild['predecessor']==str(middle)
        rigid_profile=dict(blocks=['b8:2:0'],order=3,enumeration='root-branches-v1')
        for folder,start,end,done in [(old,0,2,False),(part,2,4,True)]:
            (folder/'RIGID_INPUT.json').write_text(json.dumps(dict(binary_sha256='old-hash',**rigid_profile)))
            (folder/'SUMMARY.txt').write_text(f'rigid start branch = {start}\nrigid next branch = {end}\nrigid total branches = 4\nrigid traversal exhausted = {str(done).lower()}\n')
        assert rigid_coverage(new_folder,'old-hash',rigid_profile,old_folder)['traversal_exhausted']
        assert rigid_coverage(Path(temporary)/'empty-job','old-hash',rigid_profile,[str(old_folder),str(new_folder)])['traversal_exhausted']
    discovery=discovery_plan({})
    rows=list(itertools.chain.from_iterable(discovery['lanes'].values()))
    assert len(rows)==len({r['name'] for r in rows})
    assert {r['d'] for r in discovery['lanes']['middle']}==set(range(12,18))
    assert {r['third_order'] for r in discovery['lanes']['powers']}=={3,4}
    assert any(not r['pair_paths'] and not r['even_paths'] and r['d']==14 for r in discovery['lanes']['middle'])
    assert any((r['d']-len(r['pair_paths'])-2*len(r['even_paths']))//2%2==0 and r['third_pairs']%2==1 and r['group']=='S' for r in discovery['lanes']['middle'])
    for row in rows:
        if row['kind']=='pair':
            assert sum(row['pair_blocks'])+sum(row['pair_paths'])+2*sum(row['even_paths'])==row['d']
            cmd=discovery_command(row,'binary','out',30,'known','source',dict(next_h=0,next_third=0))
            assert cmd[cmd.index('--group')+1]==row['group']
            if row['third_order']!=2:assert '--third-pairs' not in cmd
        else:
            assert any(str(b).startswith('b') for b in row['rigid_blocks'])
            assert sum(map(block_degree,row['rigid_blocks']))==row['d']<=34
    assert discovery_command(dict(kind='sample',group='S',d=12),'binary','out',30,'known','source')[1:3]==['smart-hunt','--group']
    fresh=weak_symmetry_rows(0,Path('evidence/certificates'))
    assert len([r for r in fresh if r['family']=='square-symmetry'])==40
    assert len([r for r in fresh if r['family']=='triangle-symmetry'])==40
    old=weak_symmetry_rows(0,Path('evidence/certificates'),include_triangle=False)
    assert not any(r['family']=='triangle-symmetry' for r in old)
    assert len([r for r in old if r['family']=='square-symmetry'])==40
    pure=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True)
    assert len([r for r in pure if r['family']=='pure-symmetry'])==120
    assert len({r['name'] for r in pure})==len(pure)
    assert {r['pair_model'] for r in pure if r['family']=='pure-symmetry'}=={'A4','S4','A5'}
    next_rows=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True,include_next=True)
    adaptive=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True,include_next=True,include_frontier=True)
    finite=[r for r in adaptive if r.get('finite_construction')]
    assert len(finite)==6 and {(r['group'],r['d']) for r in finite}=={(g,d) for g in ('A','S') for d in (9,10,11)}
    assert {r['name'] for r in finite}=={r['name'] for r in weak_symmetry_rows(3,Path('evidence/certificates'),include_frontier=True) if r.get('finite_construction')}
    assert len({r['name'] for r in adaptive})==len(adaptive)
    assert adaptive[0]['source'].endswith('CERTIFIED_A25_XXYY_genus4_weak.txt')
    patterns=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True,include_next=True,include_frontier=True,include_patterns=True)
    assert patterns[0]['family']=='marked-symmetry' and patterns[0]['d']==9
    assert len({r['name'] for r in patterns})==len(patterns)
    marked=[r for r in patterns if r['family']=='marked-symmetry']
    assert {r['pair_model'] for r in marked}=={'A4','S4','A5','D5','D7','D9','D11','D13'}
    assert any(r['pair_model']=='D9' and r['pure_power']==4 and r['novel_power'] for r in marked)
    repaired=[r for r in patterns if r['family']=='resistant-glue']
    assert {3,5,12}<={r['glue_genus'] for r in repaired}
    assert not any(r['group']=='A' and r['d']==9 and r['glue_genus']==0 for r in repaired)
    for r in marked+repaired:
        cmd=discovery_command(r,'binary','out',30,'known','source')
        assert '--no-braid-fixed' not in cmd
        if r['family']=='resistant-glue':
            assert '--expand-pair-actions' in cmd and cmd[cmd.index('--glue-repair-steps')+1]==str(r['glue_repair_steps'])
        if r.get('novel_power'):assert '--novel-power' in cmd
    assert all(r['glue_repair_steps']==0 for r in repaired if r['d']<=11 or (r['d']==18 and r['glue_genus']==1))
    assert all(r['glue_repair_steps']==16 for r in repaired if r['d']==25 and r['glue_genus'] in (4,12))
    fibers=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True,include_next=True,
        include_frontier=True,include_patterns=True,include_fibers=True)
    assert fibers[0]['family']=='pair-fiber' and fibers[0]['d']==18 and fibers[0]['group']=='A'
    assert len({r['name'] for r in fibers})==len(fibers)
    refined=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True,include_next=True,
        include_frontier=True,include_patterns=True,include_fibers=True,include_refinement=True)
    refined_next=weak_symmetry_rows(1,Path('evidence/certificates'),include_pure=True,include_next=True,
        include_frontier=True,include_patterns=True,include_fibers=True,include_refinement=True)
    finite_fibers=[r for r in refined if r['family']=='pair-fiber-exact']
    assert len(finite_fibers)==15
    assert {r['name'] for r in finite_fibers}=={r['name'] for r in refined_next if r['family']=='pair-fiber-exact'}
    assert not any(r['family']=='pair-fiber' and r['d']==18 and r['group']=='A' and r['fiber_radius']==0 for r in refined)
    assert all(r.get('post_cut_base_frames') and r['name'].startswith('fiber_glue_v2_') for r in refined if r['family']=='cut-join-glue')
    assert not any(r.get('post_cut_base_frames') for r in fibers)
    extended=weak_symmetry_rows(0,Path('evidence/certificates'),include_pure=True,include_next=True,
        include_frontier=True,include_patterns=True,include_fibers=True,include_refinement=True,include_extensions=True)
    assert len({r['name'] for r in extended})==len(extended)
    assert {r['name'] for r in refined}<={r['name'] for r in extended}
    new=[r for r in extended if r['family'] in {'independent-symmetry','independent-glue','double-coset-glue','coordinated-glue','coordinated-fiber'}]
    assert new[0]['d']==9 and new[0]['family']=='independent-symmetry'
    assert {r['d'] for r in new if r['family']=='independent-symmetry'}==set(range(9,18))
    assert not any(r['group']=='A' and r['d']==9 and r.get('glue_genus')==0 for r in new)
    assert {r['allow_partition_change'] for r in new if r['family']=='coordinated-fiber'}=={False,True}
    for r in new:
        cmd=discovery_command(r,'binary','out',30,'known','source')
        assert '--no-braid-fixed' not in cmd
        if r['family']=='independent-symmetry':assert '--local-frame-cap' in cmd
        if r.get('allow_partition_change'):assert '--allow-partition-change' in cmd
        assert not r.get('finite_construction')
    for row in finite_fibers:
        cmd=discovery_command(row,'binary','out',120,'known','source')
        assert '--fiber-state-start' in cmd and '--fiber-state-count' in cmd and '--no-braid-fixed' not in cmd
    assert len({r['name'] for r in refined})==len(refined)
    from run_campaign import discovery_round_complete
    assert discovery_round_complete(False,False,['small','hybrid'],{'small','hybrid'},{},[],[])
    for stop,expired,dead,active,replay,hot in [(True,False,{'small','hybrid'},{},[],[]),
        (False,True,{'small','hybrid'},{},[],[]),(False,False,{'small'},{},[],[]),
        (False,False,{'small','hybrid'},{1:True},[],[]),(False,False,{'small','hybrid'},{},[1],[]),
        (False,False,{'small','hybrid'},{},[],[1])]:
        assert not discovery_round_complete(stop,expired,['small','hybrid'],dead,active,replay,hot)
    sources=[r for r in fibers if r['family']=='pair-fiber']
    assert any(r['d']==18 and r['source_group']=='S' and r['group']=='A' for r in sources)
    assert any(r['d']==25 and r['source_group']=='A' and r['group']=='S' for r in sources)
    assert {r['fiber_radius'] for r in sources}=={0,2,4}
    cutglue=[r for r in fibers if r['family']=='cut-join-glue']
    assert not any(r['group']=='A' and r['d']==9 and r['glue_genus']==0 for r in cutglue)
    assert any(r['group']=='S' and r['d']==9 and r['glue_genus']==0 for r in cutglue)
    assert all(r.get('screen_steps')==64 and r.get('skip_known_starts') for r in fibers if r['family']!='s3-exact')
    for row in sources+cutglue:
        cmd=discovery_command(row,'binary','out',30,'known','source')
        assert '--no-braid-fixed' not in cmd and '--skip-known-starts' in cmd
        assert cmd[cmd.index('--fiber-radius')+1]==str(row['fiber_radius'])
        assert cmd[cmd.index('--screen-steps')+1]=='64'
    from run_campaign import VerifiedDiscoveryIndex
    with tempfile.TemporaryDirectory(prefix='hurwitz-verified-index-') as temporary:
        root=Path(temporary);attempt=root/'jobs/test/attempt_000';attempt.mkdir(parents=True)
        (attempt/'JOB_INPUT.json').write_text(json.dumps(dict(lane='symmetry')))
        (attempt/'SUMMARY.txt').write_text('accepted generating seeds = 5\nsolved seeds = 4\nhard unsolved candidates = 1\nelapsed seconds = 2.5\n')
        (attempt/'hard_candidates.txt').write_text(f'group=A_9 seed={SEED}\n')
        index=VerifiedDiscoveryIndex(root);assert not index.paths and index.totals==[0,0,0]
        (attempt/'VERIFIED').touch();index.register(attempt);index.register(attempt)
        assert index.totals==[5,4,1] and index.metrics['symmetry']==[1,5,1,2.5]
        assert index.paths==set(candidate_paths(root))
        rebuilt=VerifiedDiscoveryIndex(root);assert rebuilt.metrics==index.metrics
        index.live_report();assert json.loads((root/'LIVE_REPORT.json').read_text())['verified_attempts']==1
    for row in finite:
        cmd=discovery_command(dict(row,excluded_passports='verified.tsv'),'binary','out',120,'known','source')
        assert '--no-braid-fixed' not in cmd and cmd[cmd.index('--excluded-passports')+1]=='verified.tsv'
    glue=next(r for r in adaptive if r['family']=='pair-glue' and r['d']==25)
    cmd=discovery_command(glue,'binary','out',30,'known','source')
    assert '--no-braid-fixed' not in cmd and cmd[cmd.index('--pair-genus-max')+1]=='2'
    from run_campaign import completed_s3_populations
    with tempfile.TemporaryDirectory(prefix='hurwitz-finite-s3-') as temporary:
        root=Path(temporary);attempt=root/'jobs/s3_exact_v1_A9/attempt_000';attempt.mkdir(parents=True)
        (attempt/'SUMMARY.txt').write_text('family = s3-exact\nfinite construction complete = true\nfinite population resolved = false\n')
        assert completed_s3_populations(root)==set()
        (attempt/'VERIFIED').touch();assert completed_s3_populations(root)=={'s3_exact_v1_A9'}
        construction=attempt/'construction_survivors.txt';construction.write_text(f'group=A_9 seed={SEED}\n')
        assert construction in candidate_paths(root)
        shard=root/'jobs/fiber_exact_v1_fixture_A18_s000/attempt_000';shard.mkdir(parents=True)
        (shard/'SUMMARY.txt').write_text('family = pair-fiber-exact\nfinite construction complete = false\n')
        (shard/'VERIFIED').touch();assert shard.parent.name not in completed_s3_populations(root)
        (shard/'SUMMARY.txt').write_text('family = pair-fiber-exact\nfinite construction complete = true\nfinite population resolved = false\n')
        assert shard.parent.name in completed_s3_populations(root)
    from run_campaign import completed_fiber_validations
    with tempfile.TemporaryDirectory(prefix='hurwitz-fiber-ledger-') as temporary:
        root=Path(temporary);source=root/'source.txt';source.write_text('frozen fixture')
        shard=root/'fiber_exact_v1_fixture_A18_s000';shard.mkdir();summary=shard/'SUMMARY.txt'
        summary.write_text(f'family = pair-fiber-exact\nfinite construction complete = true\nfinite population resolved = true\nsource certificate = {source}\n')
        import hashlib
        index=root/'index.json';index.write_text(json.dumps(dict(shards=[dict(name=shard.name,summary=str(summary),
            summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())])))
        assert completed_fiber_validations(index)=={shard.name}
        summary.write_text('tampered')
        try:completed_fiber_validations(index)
        except RuntimeError:pass
        else:raise AssertionError('Changed finite evidence was silently accepted')
    assert len({r['name'] for r in next_rows})==len(next_rows)
    assert {r['pure_power'] for r in next_rows if r['family']=='pure-symmetry'}=={1,2,3}
    assert {r['pair_model'] for r in next_rows if r['family']=='pure-symmetry'}=={'A4','S4','A5','PSL27'}
    assert {r['d'] for r in next_rows if r.get('pair_model')=='PSL27'}==set(range(9,19))
    glues=[r for r in next_rows if r['family']=='dessin-glue']
    assert not any(r['group']=='A' and r['d']==9 and r['glue_genus']==0 for r in glues)
    for row in next_rows:
        cmd=discovery_command(row,'binary','out',30,'known','source')
        if row.get('pure_power'):assert cmd[cmd.index('--pure-power')+1]==str(row['pure_power'])
        if row['family']=='dessin-glue':
            assert '--no-braid-fixed' not in cmd
            assert cmd[cmd.index('--glue-genus')+1]==str(row['glue_genus'])
    for row in pure:
        if not row['family'].startswith('pure-'):continue
        cmd=discovery_command(row,'binary','out',30,'known','source')
        assert '--no-braid-fixed' not in cmd
        if row['family']=='pure-symmetry':assert cmd[cmd.index('--pair-model')+1]==row['pair_model']
        else:assert '--certificate' in cmd
    assert len({r['name'] for r in fresh})==len(fresh)
    assert {r['d'] for r in fresh if r['lane']=='small-symmetry'}=={9,10,11}
    assert not {r['name'] for r in fresh}&{r['name'] for r in weak_symmetry_rows(1,Path('evidence/certificates'))}
    for row in fresh:
        cmd=discovery_command(row,'binary','out',30,'known','source')
        assert '--no-braid-fixed' not in cmd and '--weak-only' in cmd
        if row['family']=='symmetry-centralizer':assert '--certificate' in cmd
    with tempfile.TemporaryDirectory(prefix='hurwitz-new-fixed-sources-') as temporary:
        root=Path(temporary);known=root/'known';resolved=root/'resolved';profiles=root/'classifications'
        for folder in (known,resolved,profiles):folder.mkdir()
        source=resolved/'CERTIFIED_fresh_weak.txt'
        source.write_bytes(Path('evidence/certificates/CERTIFIED_S18_XXXX_genus1_weak.txt').read_bytes())
        profile=profiles/'fresh.tsv'
        baseline=len(weak_symmetry_rows(0,known))
        profile.write_text('any_braid_fixed_states\n0\n')
        assert len(weak_symmetry_rows(0,known,resolved))==baseline
        pure_baseline=len(weak_symmetry_rows(0,known,include_pure=True))
        assert len(weak_symmetry_rows(0,known,resolved,include_pure=True))==pure_baseline+1
        profile.write_text('any_braid_fixed_states\n2\n')
        assert len(weak_symmetry_rows(0,known,resolved))==baseline+1
        (known/'CERTIFIED_duplicate_weak.txt').write_bytes(source.read_bytes())
        assert len(weak_symmetry_rows(0,known,resolved))==baseline+1
    with tempfile.TemporaryDirectory(prefix='hurwitz-inherited-rigid-alias-') as temporary:
        root=Path(temporary);parent=root/'parent';middle=root/'middle'
        old_row=dict(name='original',kind='rigid',rigid_blocks=['old-model'],order=5)
        for folder,plan in ((parent,dict(lanes={'rigid':[old_row]})),
                            (middle,dict(lanes={'rigid':[]},predecessor=str(parent)))):
            (folder/'inputs').mkdir(parents=True)
            (folder/'STATE.json').write_text(json.dumps(dict(plan='discovery-frontier',status='stopped',binary_sha256='original-hash')))
            (folder/'inputs/DISCOVERY_PLAN.json').write_text(json.dumps(plan))
        attempt=parent/'jobs/original/attempt_000';attempt.mkdir(parents=True)
        (attempt/'VERIFIED').touch()
        (attempt/'RIGID_INPUT.json').write_text(json.dumps(dict(binary_sha256='original-hash',blocks=['old-model'],order=5,enumeration='root-branches-v1')))
        (attempt/'SUMMARY.txt').write_text('rigid start branch = 0\nrigid next branch = 1\nrigid total branches = 1\nrigid traversal exhausted = true\n')
        previous=[p for p in completed_rigid_profiles(middle) if p['blocks']==['old-model']]
        assert len(previous)==1
        fresh=dict(lanes={'rigid':[dict(rigid_blocks=['alias-model'],order=5,name='renamed')]},excluded=[])
        deduplicate_representations(fresh,{'old-model':'same-action','alias-model':'same-action'},previous)
        assert not fresh['lanes']['rigid'] and len(fresh['excluded'])==1
    duplicate=dict(lanes={'rigid':[dict(name='a',rigid_blocks=['b8:2:0'],order=3),dict(name='b',rigid_blocks=['b8:2:1'],order=3),dict(name='c',rigid_blocks=['b8:2:2'],order=0)]},excluded=[])
    keys={f'b8:2:{i}':'same exact representation' for i in range(3)}
    result=deduplicate_representations(duplicate,keys,[dict(blocks=['b8:2:2'],order=6,evidence='completed')])
    assert [r['name'] for r in result['lanes']['rigid']]==['c']  # order six covers order three, NOT unrestricted
    with tempfile.TemporaryDirectory(prefix='hurwitz-new-cursors-') as temporary:
        folder=Path(temporary);profile=dict(blocks=['b8:2:0'],order=3,enumeration='root-branches-v1')
        for i,(start,end,done) in enumerate([(0,2,False),(2,2,False),(2,4,True)]):
            attempt=folder/f'attempt_{i:03d}';attempt.mkdir();(attempt/'VERIFIED').touch()
            (attempt/'RIGID_INPUT.json').write_text(json.dumps(dict(binary_sha256='test',**profile)))
            (attempt/'SUMMARY.txt').write_text(f'rigid start branch = {start}\nrigid next branch = {end}\nrigid total branches = 4\nrigid traversal exhausted = {str(done).lower()}\n')
        assert rigid_coverage(folder,'test',profile)['traversal_exhausted']
        try:rigid_coverage(folder,'changed',profile)
        except RuntimeError:pass
        else:raise AssertionError('changed representation binary accepted')
        path=folder/'attempt_002/SUMMARY.txt';path.write_text(path.read_text().replace('rigid start branch = 2','rigid start branch = 3'))
        try:rigid_coverage(folder,'test',profile)
        except RuntimeError:pass
        else:raise AssertionError('rigid branch gap accepted')
    followup=structured_followup({('A',10,tuple(sorted(['2.2.2.2']*3+['3.3.3']))):'isolated test exclusion'})
    assert len({r['name'] for r in followup['exact']})==len(followup['exact'])
    excluded={(r['group'],r['d'],tuple(sorted(r['passport'].split('/')))) for r in followup['excluded']}
    def partitions(n,minimum=1):
        if not n:yield ()
        for first in range(minimum,n+1):
            for rest in partitions(n-first,first):yield (first,)+rest
    for group,d,genera in [('A',9,(1,)),('A',10,(0,1)),('A',11,(0,1)),('S',11,(0,1))]:
        types=[p for p in partitions(d) if len(p)<d and (group=='S' or (d-len(p))%2==0)]
        expected=set()
        for four in itertools.combinations_with_replacement(types,4):
            if sorted(Counter(four).values()) not in ([4],[2,2],[1,3]):continue
            if group=='S' and all((d-len(p))%2==0 for p in four):continue
            if 1-d+sum(d-len(p) for p in four)//2 not in genera:continue
            if sum(d-len(p) for p in four)%2:continue
            expected.add(tuple(sorted('.'.join(str(x) for x in p if x>1) for p in four)))
        generated=structured_slices(group,(d,),genera,quotient=True)
        assert {tuple(sorted(r['passport'].split('/'))) for r in generated}==expected
        planned={tuple(sorted(r['passport'].split('/'))) for r in followup['exact'] if (r['group'],r['d'])==(group,d)}
        assert planned=={p for p in expected if (group,d,p) not in excluded}
    # Independently enumerate arbitrary type multisets, including XXYZ/XYZW.
    mission=small_passport_plan({});rows=mission['rows'];all_types=list(partitions(9))
    expected=set()
    for four in itertools.combinations_with_replacement(all_types,4):
        if any(not 3<=len(p)<9 for p in four):continue
        indices=[9-len(p) for p in four]
        if sum(indices)%2:continue
        group='S' if any(i%2 for i in indices) else 'A'
        genus=1-9+sum(indices)//2
        if genus not in ((1,2,3,4) if group=='A' else (0,1)):continue
        expected.add((group,tuple(sorted('.'.join(str(x) for x in p if x>1) for p in four))))
    assert {(r['group'],tuple(sorted(r['passport'].split('/')))) for r in rows}==expected
    target=next(r for r in rows if r['group']=='A' and r['genus']==1 and
                tuple(sorted(r['passport'].split('/')))==tuple(sorted(['3.3.3','3','2.2.2.2','2.2.5'])))
    assert target['phase']==0
    key=(target['group'],9,tuple(sorted(target['passport'].split('/'))))
    reduced=small_passport_plan({key:'verified_fixture'})
    assert len(reduced['rows'])==len(rows)-1 and reduced['excluded'][0]['evidence']=='verified_fixture'
    assert len({r['name'] for r in rows})==len(rows)
    queue=PassportQueue();first,second=rows[:2]
    queue.push(first);queue.push(second)
    assert queue.pop()['name']==first['name']
    queue.push(first,10);assert queue.pop()['name']==second['name']
    assert queue.pop()['name']==first['name'] and queue.pop() is None
    cmd=discovery_command(target,'binary','out',15,'known','source',dict(next=128),seed=19)
    assert cmd[1]=='passport-exact' and cmd[cmd.index('--start')+1]=='128'
    assert '--centralizer-rows' in cmd and '--weak-only' in cmd and '--no-braid-fixed' not in cmd
    assert cmd[cmd.index('--probe-witness-limit')+1]=='256'
    assert all(r['centralizer'] and not (r['d']==9 and r['genus']==0) for r in followup['exact'])
    assert {r['d'] for r in followup['mechanisms']}=={15,16,17}
    assert all(r['group']=='A' and r['resumable'] for r in followup['mechanisms'])
    for r in followup['mechanisms']:
        m=(r['d']-len(r.get('pair_paths',[]))-2*len(r.get('even_paths',[])))//2
        assert r['third_pairs']!=m and r['third_pairs']%2==0 and r['third_pairs']<=r['d']//2
        assert 2*m+r['third_pairs']>=r['d']+1
        components=len(r['pair_blocks'])+len(r.get('pair_paths',[]))+2*len(r.get('even_paths',[]))
        assert 1<components<=r['third_pairs']+1
    assert all(r['reason'] for r in followup['mechanism_excluded'])
    with tempfile.TemporaryDirectory(prefix='hurwitz-exact-coverage-') as temporary:
        folder=Path(temporary)
        assert not exact_coverage(folder,'test','3/3/3/3')['enumeration_complete']
        for i,(start,end,hard) in enumerate([(0,10,1),(10,20,0)]):
            attempt=folder/f'attempt_{i:03d}';attempt.mkdir();(attempt/'VERIFIED').touch()
            (attempt/'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='test',passport='3/3/3/3',enumeration='centralizer-v1')))
            (attempt/'SUMMARY.txt').write_text(f'scan start = {start}\nscan next = {end}\nscan total candidate pairs = 20\nscan enumeration = centralizer-v1\nhard unsolved candidates = {hard}\n')
        c=exact_coverage(folder,'test','3/3/3/3')
        assert c['enumeration_complete'] and not c['search_complete']
        import run_campaign as runner
        memo={}
        with patch.object(runner,'extend_exact_coverage',wraps=runner.extend_exact_coverage) as validate:
            cached=exact_coverage(folder,'test','3/3/3/3',memo=memo)
            assert cached==c and validate.call_count==2
            cached['search_complete']=True
            assert exact_coverage(folder,'test','3/3/3/3',memo=memo)==c
            assert validate.call_count==2
            try:exact_coverage(folder,'changed','3/3/3/3',memo=memo)
            except RuntimeError:pass
            else:raise AssertionError('cached prefix hid a changed binary')
        for digest,passport in [('changed','3/3/3/3'),('test','2.2/2.2/2.2/2.2')]:
            try:exact_coverage(folder,digest,passport)
            except RuntimeError:pass
            else:raise AssertionError('unsafe exact cursor reuse accepted')
        last=folder/'attempt_001/SUMMARY.txt'
        last.write_text(last.read_text().replace('scan start = 10','scan start = 11'))
        try:exact_coverage(folder,'test','3/3/3/3')
        except RuntimeError:pass
        else:raise AssertionError('exact coverage gap accepted')
    with tempfile.TemporaryDirectory(prefix='hurwitz-passport-inheritance-') as temporary:
        parent=Path(temporary)/'old';current=Path(temporary)/'new'
        for root in [parent,current]:(root/'jobs').mkdir(parents=True)
        row=dict(target);old=parent/'jobs'/row['name'];old.mkdir()
        for folder,start,end in [(old,0,128),(current/'jobs'/row['name'],128,256)]:
            attempt=folder/'attempt_000';attempt.mkdir(parents=True);(attempt/'VERIFIED').touch()
            (attempt/'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='old-hash',passport=row['passport'],enumeration='centralizer-v1')))
            (attempt/'SUMMARY.txt').write_text(f'scan start = {start}\nscan next = {end}\nscan total candidate pairs = 256\nscan enumeration = centralizer-v1\nhard unsolved candidates = 0\n')
        (parent/'STATE.json').write_text(json.dumps(dict(status='stopped',plan='discovery-frontier',binary_sha256='old-hash')))
        (parent/'inputs').mkdir();(parent/'inputs/DISCOVERY_PLAN.json').write_text(json.dumps(dict(lanes={'exact':[row]},excluded=[])))
        inherited=inherit_discovery(dict(lanes={'exact':[dict(row)]},excluded=[]),parent)['lanes']['exact'][0]
        assert inherited['enumeration_sha256']=='old-hash'
        c=discovery_coverage(inherited,current/'jobs'/row['name'],'different-new-hash')
        assert c['next']==256 and c['traversal_exhausted'] and c['search_complete']
        try:discovery_coverage(dict(inherited,enumeration_sha256='wrong'),current/'jobs'/row['name'],'different-new-hash')
        except RuntimeError:pass
        else:raise AssertionError('changed inherited passport binary accepted')
        empty=Path(temporary)/'zero';attempt=empty/'attempt_000';attempt.mkdir(parents=True);(attempt/'VERIFIED').touch()
        (attempt/'SCAN_INPUT.json').write_text(json.dumps(dict(binary_sha256='zero',passport='3/3/3/3',enumeration='centralizer-v1')))
        (attempt/'SUMMARY.txt').write_text('scan start = 0\nscan next = 0\nscan total candidate pairs = 0\nscan enumeration = centralizer-v1\nhard unsolved candidates = 0\n')
        assert exact_coverage(empty,'zero','3/3/3/3')['search_complete']
    small=small_symmetric_slices()
    assert len(small)==320 and sum(r['total'] for r in small)==97250587651
    assert sum(r['d']==9 for r in small)==118 and sum(r['d']==10 for r in small)==202
    assert max(r['total'] for r in small)==3657830400
    pilot=small_symmetric_slices(5_000_000)
    assert len(pilot)==122 and sum(r['total'] for r in pilot)==101235393
    reduced={'scan next':'60','scan raw class pairs':'100','scan total candidate pairs':'60','scan range resolved':'true'}
    assert exact_scan_complete(reduced,(0,100)) and not exact_scan_complete(reduced,(20,80))
    assert not exact_scan_complete(dict(reduced,**{'scan next':'100','scan range resolved':'false'}),(0,100))
    assert all(r['group']=='S' and r['genus'] in (0,1) for r in small)
    assert any(r['passport']=='2.2.4/2.2.2.2.2/2.2.2.2.2/2.2.4' and r['tier']==0 for r in small)
    assert any('a25-three-one' in r['mechanisms'] for r in small)
    assert len({(r['d'],tuple(sorted(r['passport'].split('/')))) for r in small})==len(small)
    new_pairs=counterexample_pair_jobs()
    assert len(new_pairs)==len({j['name'] for j in new_pairs})==18
    assert {j['d'] for j in new_pairs}=={20,22,24,26,28,30,32}
    assert all(j['resumable'] and sum(j['pair_blocks'])==j['d'] and
               (j['group']=='A')==(j['d']%4==0) for j in new_pairs)
    structure=structure_jobs()
    assert len({j['name'] for j in structure})==len(structure)
    assert {j['operation'] for j in structure if 'operation' in j}=={'braid-symmetry','independent','research','rayleigh'}
    assert {j['d'] for j in structure if j.get('even_paths')}==set(range(12,18))
    for j in structure:
        if j.get('even_paths'):
            assert j['resumable']
            assert sum(j['pair_blocks'])+sum(j['pair_paths'])+2*sum(j['even_paths'])==j['d']
            fixed=len(j['pair_paths'])+2*len(j['even_paths'])
            m=(j['d']-fixed)//2
            assert (j['group']=='S')==(m%2==1) and 3*m>=j['d']-1
    with tempfile.TemporaryDirectory(prefix='hurwitz-pair-coverage-') as temporary:
        job=Path(temporary)
        assert not pair_coverage(job,'test')['search_complete']
        for i,(start,end,exhausted) in enumerate([((0,0),(3,7),False),((3,7),(8,0),True)]):
            attempt=job/f'attempt_{i:03d}';attempt.mkdir()
            (attempt/'VERIFIED').touch()
            (attempt/'PAIR_INPUT.json').write_text(json.dumps(dict(binary_sha256='test',cycles=[2],paths=[],even_paths=[2])))
            (attempt/'SUMMARY.txt').write_text(f'scan start h = {start[0]}\nscan start third = {start[1]}\nscan next h = {end[0]}\nscan next third = {end[1]}\nscan range resolved = true\ntraversal exhausted = {str(exhausted).lower()}\n')
        assert pair_coverage(job,'test')['search_complete']
        for digest,profile in [('changed',None),('test',dict(cycles=[],paths=[],even_paths=[2]))]:
            try:pair_coverage(job,digest,profile)
            except RuntimeError:pass
            else:raise AssertionError('unsafe cursor reuse accepted')
        try:pair_coverage(job,'test',dict(cycles=[2],paths=[],even_paths=[2],third_pairs=0))
        except RuntimeError:pass
        else:raise AssertionError('changed third type accepted')
        last=job/'attempt_001/SUMMARY.txt'
        last.write_text(last.read_text().replace('scan start third = 7','scan start third = 8'))
        try:pair_coverage(job,'test')
        except RuntimeError:pass
        else:raise AssertionError('coverage gap accepted')
    mechanisms=mechanism_jobs()
    assert len({j['name'] for j in mechanisms})==len(mechanisms)
    assert mechanisms[0]['position_spectrum'] and mechanisms[1]['order']==0
    assert all(4<=j['d']<=34 for j in mechanisms)
    assert {j['d'] for j in mechanisms if j.get('pair_paths')}==set(range(12,18))
    for j in mechanisms:
        if j.get('pair_paths'):
            assert j['d']==sum(j['pair_blocks'])+sum(j['pair_paths'])
            assert all(x%2==1 for x in j['pair_paths'])
            assert all(x%2==0 for x in j['pair_blocks'])
        if j.get('rigid_blocks'):
            assert sum(map(block_degree,j['rigid_blocks']))==j['d']
    for group,d,m,exceptional in s18_family_jobs():
        index=sum(int(l)-1 for l in exceptional.split('.'))
        assert 3*m+index==2*d and 2*m==d and (group=='S')==(m%2==1)
    glued=rigid_glue_jobs()
    assert len(glued)==len(set(glued))
    assert all(4<=sum(b)<=34 and set(b)<={4,8,16,*range(3,35,3)} for b in glued)
    assert {(15,16),(15,15),(18,),(12,16)}<=set(glued)
    assert {(4,9),(3,16),(8,16),(9,9,16),(3,3,3,3,16)}<=set(glued)
    general=rigid_glue_jobs(3)
    assert set(glued)<set(general) and len(general)==len(set(general))
    assert {(4,),(6,),(8,),(3,3)}<=set(general)
    assert all(4<=sum(b)<=34 for b in general)
    for order in (2,3):
        lifted=rigid_glue_jobs(order,True)
        assert len(lifted)==len(set(lifted))
        assert all(4<=sum(map(block_degree,b))<=34 and any(str(x).startswith('t') for x in b) for b in lifted)
    structural=braid_fixed_jobs()
    assert structural[0]==('A',25,10,0)
    assert len(structural)==len(set(structural))
    assert any(g==1 and d<12 for _,d,_,g in structural)
    for group,d,m,g in structural:
        assert 4<=d<=25 and 1<=m<=d//2 and g in (0,1)
        assert 3<=3*m-d+2-2*g<=d and (group=='S')==(m%2==1)
        assert g>0 or d>11
    small = degree11_involution_slices()
    assert len(small) == 15 and len(set(small)) == 15
    assert sum(g == 'A' for g, _, _ in small) == 8
    assert sum(g == 'S' for g, _, _ in small) == 7
    for group, profile, total in small:
        types = profile.split('/')
        assert types[0] == types[1] == types[2]
        indices = [sum(int(n)-1 for n in t.split('.')) for t in types]
        assert sum(indices) == 20
        assert all(i % 2 == (group == 'S') for i in indices)
        assert total == class_size(types[1], 11)*class_size(types[2], 11)
        assert types[3] not in ('3.4.4', '3.3.5') or group != 'A'
    with tempfile.TemporaryDirectory(prefix='hurwitz-checkpoint-') as temporary:
        progress = Path(temporary) / 'CENSUS_PROGRESS.tsv'
        progress.write_text('passport\ttotal_pairs\tverified_scanned_prefix\tverified_resolved_prefix\tcomplete\n'
                            'p0\t10\t10\t10\tTrue\n'
                            'p1\t20\t5\t5\tFalse\n')
        assert census_progress(Path(temporary)) == (1, 2, 15, 30)
        line = format_status_bar({'root': temporary, 'forever': True, 'phase': 'A11 deform strong',
                                  'completed_jobs': 3, 'threads': 8, 'job': 'r000'}, False)
        assert 'A9 exact' in line and 'phase A11 deform strong' in line and 'unbounded campaign' in line
    assert scheduled_threads(dt.datetime(2026,9,7,0,59,tzinfo=CORE_TIMEZONE)) == 8
    assert scheduled_threads(dt.datetime(2026,9,7,1,0,tzinfo=CORE_TIMEZONE)) == 4
    assert scheduled_threads(dt.datetime(2026,9,7,15,59,tzinfo=CORE_TIMEZONE)) == 4
    assert scheduled_threads(dt.datetime(2026,9,7,16,0,tzinfo=CORE_TIMEZONE)) == 8
    assert seconds_until_thread_change(dt.datetime(2026,9,7,0,59,tzinfo=CORE_TIMEZONE)) == 60
    assert seconds_until_thread_change(dt.datetime(2026,9,7,15,59,tzinfo=CORE_TIMEZONE)) == 60
    passport, genus = seed_profile(SEED,9)
    assert passport == '/'.join(['3.3.3']*4) and genus == 4
    assert class_size('3.3',9)==3360 and class_size('3.3.3',9)==2240
    known = [dict(group='A',d='9',seed=SEED)]
    jobs = research_jobs([9,10,11,8,16],known,0)
    assert jobs[0] == ('A',10,'deform',False,None)
    assert jobs.index(('A',11,'deform',False,None)) < jobs.index(('A',9,'passport',False,passport))
    assert all(row[1]<=11 for row in jobs[:5])
    assert all(not (group=='A' and d==8) for group,d,*_ in jobs)
    assert ('A',9,'passport',True,passport) in jobs
    assert ('A',11,'deform',False,None) in jobs
    a13 = research_jobs([13],known,0)
    assert ('A',13,'passport',False,'/'.join(['3.3.3']*4)) in a13
    assert not any(row[1]==16 and row[2]=='deform' for row in jobs)
    assert research_jobs([9],known,0) == research_jobs([9],known,0)
    assert set(research_jobs([9],known,1)) == set(research_jobs([9],known,2))
    slices=exact_slices(known)
    assert ('3.3/3.3/3.3/3.3',11289600) in slices
    assert ('3.3.3/3.3.3/3.3.3/3.3.3',5017600) in slices
    resolved = [p for p,_ in slices]
    filtered = research_jobs([9,10,11,13,16],known,0,resolved)
    assert not any(g=='A' and d==9 and f=='passport' for g,d,f,*_ in filtered)
    assert ('A',9,'deform',False,None) in filtered
    assert ('A',10,'passport',False,'3.3/3.3/3.3/3.3.3') in filtered
    assert filtered.index(('A',13,'passport',False,'/'.join(['3.3.3']*4))) < next(i for i,r in enumerate(filtered) if r[1]==16)
    # Sorting execution preserves stored IDs and the complete workload on resume.
    planned=exact_slices(known,True)
    reordered=sorted(enumerate(planned),key=exact_priority)
    assert dict(reordered)==dict(enumerate(planned))
    assert all(exact_priority(item)[0]==False for item in reordered[:8])
    assert sum(total for _,total in planned)==978564937 and len(planned)==380
    all_profiles=genus_zero_profiles(9)
    scheduled={tuple(sorted(p.split('/'))) for p,_ in exact_slices(known,True)}
    assert {tuple(sorted(p.split('/'))) for p in all_profiles} <= scheduled
    for profile in all_profiles:
        assert sum(sum(int(n)-1 for n in p.split('.')) for p in profile.split('/'))==16
    assert max(total for _,total in exact_slices(known,True))<=60_000_000
    assert len(exact_slices(known,True))==len({p for p,_ in exact_slices(known,True)})
    with tempfile.TemporaryDirectory(prefix='hurwitz-runner-check-') as temporary:
        root = Path(temporary)
        earlier, later = root/'old',root/'new'
        earlier.mkdir(); later.mkdir()
        hard = earlier/'hard_candidates.txt'
        hard.write_text(f'group=A_9 status=capped seed={SEED}\n')
        new = later/'hard_candidates.txt'
        new.touch()
        results = later/'RESULTS.tsv'
        def resolved(status):
            results.write_text(f'status\tquotient_states\tweak\tseed\n{status}\t1\ttrue\t{SEED}\n')
        resolved('weak_found')
        assert candidates_from_files([hard,new]) == {'A_9':{SEED}}
        assert candidates_from_files([hard,new],weak_only=True) == {}
        resolved('strong_found')
        assert candidates_from_files([hard,new]) == {}
        resolved('timeout')
        assert candidates_from_files([hard,new]) == {'A_9':{SEED}}
        # Header reads stay cheap; component identity must ignore filenames.
        catalog = root/'known'
        catalog.mkdir()
        certificate = catalog/'CERTIFIED_A9_fixture.txt'
        certificate.write_text(f'HURWITZ_SPARSE_CERTIFICATE_V1\ngroup=A\nd=9\nseed={SEED}\nstates=36\nstatus=certified_strong_failure\nweak_witness={SEED}\nORDER\t504\t6\nSTATE\t1\t2\t3\n')
        copy=root/'jobs'/'rediscovery'/'attempt_000'/'CERTIFIED_different_name.txt'
        copy.parent.mkdir(parents=True)
        copy.write_text(certificate.read_text())
        assert certificate_identity(copy)==certificate_identity(certificate)
        assert certificate_header(certificate)['max_pair_order']==504
        report(root)
        with (root/'CATALOG.tsv').open() as stream:
            row = next(csv.DictReader(stream,delimiter='\t'))
        assert row['genus']=='4' and row['quotient_states']=='36'
        assert 'New components beyond starting catalog: 0' in (root/'REPORT.txt').read_text()
    print('runner checks passed')


if __name__ == '__main__':
    check()
