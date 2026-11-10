"""Small source-bound handoff records; no discovery or external publication."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
from run_smooth_degree9 import sha


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','checks','replay','lattice','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--replace-unpublished',action='store_true')
    p.add_argument('--rust-checks',action='store_true')
    args=p.parse_args();args.out.mkdir(exist_ok=True)
    if (args.out/'CHECKS.json').exists():
        project=Path(__file__).resolve().parents[2];manifest=json.loads((project/'evidence/manifest.json').read_text())
        if not args.replace_unpublished or any((project/e['path']).resolve().is_relative_to(args.out.resolve()) for e in manifest['sources']):raise ValueError('do not overwrite published handoff')
    validations=json.loads(args.checks.read_text());plan=json.loads((args.root/'PLAN.json').read_text())
    if not validations['passed'] or plan['validation_sha256']!=sha(args.checks):raise ValueError('launch validation binding')
    bindings=[args.checks,args.checks.parent/'SOURCE_HASHES.json',args.checks.parent/'BENCHMARK.json',
        args.root/'PLAN.json',args.root/'LAUNCH.json',args.root/'NOVELTY_INPUT.json',args.lattice,
        args.replay/'RESULTS.json',args.replay/'INPUTS.json']
    replay=json.loads((args.replay/'RESULTS.json').read_text())
    replay_summary={k:v for k,v in replay.items() if k not in ['results']}
    replay_summary['partial_overlap_groups']=replay_summary.pop('deduplicated_components',replay_summary.get('partial_overlap_groups'))
    replay_summary['interpretation']='No overlap of saved partial graphs was observed. These139 negative occurrences need not be globally different braid components.'
    save(args.out/'REPLAY.json',replay_summary)
    lattice=json.loads(args.lattice.read_text());entries=lattice['components']
    if lattice['queue_sha256']!=plan['profile_input_sha256'] or not lattice['complete']:raise ValueError('lattice queue binding')
    hist=Counter();core=Counter();genus0=[];idcounts=Counter();core_id=Counter();all_id=Counter();identity_free_components=0
    for e in entries:
        b=e['block_lattice']
        if not b['complete']:raise ValueError('incomplete block lattice')
        hist.update(s['genus'] for s in b['proper_systems'])
        core.update(s['genus'] for i,s in enumerate(b['proper_systems']) if i in b['primitive_quotient_indices'])
        if any(s['genus']==0 for s in b['proper_systems']):genus0.append(e['queue_index'])
        idcounts.update(len(s['identity_entries']) for s in b['proper_systems'] if s['genus']==0)
        all_id.update(len(s['identity_entries']) for s in b['proper_systems'])
        core_id.update((s['genus'],len(s['identity_entries'])) for i,s in enumerate(b['proper_systems']) if i in b['primitive_quotient_indices'])
        identity_free_components+=any(not s['identity_entries'] for s in b['proper_systems'])
    save(args.out/'LATTICE.json',dict(components=len(entries),complete_lattices=len(entries),primitive=sum(e['block_lattice']['primitive'] for e in entries),
        proper_block_systems=sum(hist.values()),quotient_genus_system_counts=dict(hist),primitive_quotient_genus_counts=dict(core),
        components_with_genus_zero_quotient=len(genus0),genus_zero_identity_entry_counts=dict(idcounts),
        proper_quotient_identity_counts=dict(all_id),components_with_identity_free_intermediate_quotient=identity_free_components,
        primitive_core_genus_identity_counts={f'{h},{identities}':count for (h,identities),count in core_id.items()},
        example=entries[genus0[0]],report=dict(path=str(args.lattice),sha256=sha(args.lattice)),
        scope='All invariant partitions of these7169 transitive representatives. Quotient genera are not component counts. No geometric classification.'))
    save(args.out/'VALIDATION.json',validations)
    save(args.out/'BENCHMARK.json',json.loads((args.checks.parent/'BENCHMARK.json').read_text()))
    affine=args.root/'AFFINE_COVERAGE_agl16_16.json'
    if affine.exists():
        report=json.loads(affine.read_text())
        for item in report['proofs']:
            path=Path(item['path']);part=json.loads(path.read_text())
            if sha(path)!=item['sha256'] or not part['passed']:raise ValueError('affine independent proof changed')
            if any(sha(path.parent/name)!=digest for name,digest in part['sources'].items()):raise ValueError('affine population source changed')
        save(args.out/'AFFINE.json',report);bindings.append(affine)
    if args.rust_checks:
        log=args.checks.parent/'RUST_TESTS.log'
        if log.exists():raise ValueError('do not replace a prior Rust test log')
        project=Path(__file__).resolve().parents[2]
        result=subprocess.run(['cargo','test','--release','--bin','smooth_forest_fast'],cwd=project,capture_output=True,text=True)
        log.write_text(result.stdout+result.stderr)
        match=re.search(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',log.read_text())
        if result.returncode or not match:raise ValueError('Rust regressions failed')
        save(args.out/'RUST_CHECKS.json',dict(passed=True,tests_passed=int(match[1]),failed=int(match[2]),ignored=int(match[3]),
            command=['cargo','test','--release','--bin','smooth_forest_fast'],log=dict(path=str(log),sha256=sha(log)),
            source_hashes=plan['source_hashes']))
        bindings.append(log)
    state=json.loads((args.root/'STATE.json').read_text());launch=json.loads((args.root/'LAUNCH.json').read_text())
    native=subprocess.check_output(['/bin/launchctl','list',launch['launchd_label']],text=True)
    if state['status']!='searching' or f'"PID" = {launch["pid"]};' not in native or state['pid']!=launch['pid']:raise ValueError('controller/native/searching mismatch')
    jobs=sorted(p for p in args.root.glob('guided_*') if (p/'VERIFIED').exists())
    jobs+=sorted(p for p in args.root.glob('replay_*') if (p/'VERIFIED').exists())
    if not jobs:raise ValueError('first independently verified job not finished')
    first=jobs[0];r=json.loads((first/'CHECKS.json').read_text());bindings.append(first/'CHECKS.json')
    save(args.out/'LAUNCH.json',dict(root=str(args.root),state=state,launch=launch,native_ownership=native,
        first_verified_job=str(first),first_checks_sha256=sha(first/'CHECKS.json'),first_construction=r['construction'],
        scope='Startup snapshot, not a final result tally. Live STATE/VERIFIED authoritative.'))
    save(args.out/'CHECKS.json',dict(passed=True,human_scientific_review='pending',no_tex_edit=True,
        bindings=[dict(path=str(path),sha256=sha(path)) for path in bindings],
        limits='Biased proposals are not censuses; named affine omissions independently checked. No weak/primitive/geometry novelty inferred.'))
    print('RECORDED='+str(args.out),flush=True)


if __name__=='__main__':main()
