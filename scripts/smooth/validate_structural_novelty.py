"""Bounded constructor, omission, cache and lifecycle validation before launch.

Discovery draws are biased proposals, not uniform samples or complete passports.
The independent finite-action audit reconstructs omissions as well as journal
witnesses. Benchmark pairs use identical rank intervals and literal outputs.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import sys
import time

from audit_guided_native import audit_guided_native
from run_smooth_degree9 import sha,save

PROJECT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['binary','auditor','inputs','out']:p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();root=args.out.resolve();root.mkdir()
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binary=args.binary.resolve();auditor=args.auditor.resolve();inputs=args.inputs.resolve()
    env=os.environ.copy();env.update(DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    hashes={};source=root/'source';source.mkdir()
    for path in [*sorted((PROJECT/'src').rglob('*.rs')),*sorted((PROJECT/'scripts').rglob('*.py'))]:
        rel=path.relative_to(PROJECT);target=source/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);hashes[str(rel)]=sha(target)
    save(root/'SOURCE_HASHES.json',hashes)
    shutil.copy2(binary,root/'smooth_forest_fast');shutil.copy2(auditor,root/'smooth_negative_audit')
    tests=[]
    for name in ['test_structural_novelty','test_smooth_block_lattice','test_saved_forest_replay',
                 'test_mechanism_breaker','test_native_components','test_theory_frontier','test_discovery_design',
                 'test_guided_native','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        result=subprocess.run([sys.executable,str(_scripts_path.find(f'{name}.py'))],env=env,capture_output=True,text=True)
        log=root/f'{name}.log';log.write_text(result.stdout+result.stderr)
        if result.returncode:raise RuntimeError(name+'\n'+log.read_text())
        tests.append(dict(name=name,log_sha256=sha(log)));print('PASS',name,flush=True)
    control=subprocess.check_output([str(binary),'--check-controls',str(inputs/'controls.tsv')],text=True)
    (root/'CONTROLS.log').write_text(control);print(control.strip(),flush=True)
    exclusions=json.loads((Path('data/smooth_fast_20261013_113427_284452')/'PRIMITIVE_EXCLUSIONS.json').read_text())
    for item in exclusions.values():
        proof=Path(item['path']);r=json.loads(proof.read_text())
        if sha(proof)!=item['sha256'] or not r['passed'] or not r['independent_complete_population'] or not r['no_unresolved']:raise ValueError('prior finite proof binding')
        if any(sha(proof.parent/name)!=digest for name,digest in r['sources'].items()):raise ValueError('prior finite source changed')
    novelty=json.loads((inputs/'NOVELTY_INPUT.json').read_text())
    for item in novelty['negative_frontier']:
        folder=Path(item['folder']);path=folder/'CHECKS.json';r=json.loads(path.read_text())
        if sha(path)!=item['checks_sha256'] or r['independent_native']['closed'] or r['independent_native']['unresolved']:raise ValueError('negative frontier not bound')
        audit_guided_native(folder,auditor,native_components=True)
    pilots=[]
    def job(name,n,lane,extra=(),attempts=10000):
        folder=root/name;began=time.monotonic()
        command=[str(binary),'--d',str(n),'--lane',lane,'--frame-reuse','1','--seed',str(20261013000+n),
                 '--seconds','60','--attempts',str(attempts),'--cap','10000','--closure-seconds','2',
                 '--full-closure-budget','--out',str(folder),*map(str,extra)]
        subprocess.run(command,check=True,capture_output=True)
        r=audit_guided_native(folder,auditor,native_components=True)
        save(folder/'CHECKS.json',r);(folder/'VERIFIED').write_text('Independent bounded validation\n')
        entry=dict(name=name,degree=n,lane=lane,construction=r['construction'],checks_sha256=sha(folder/'CHECKS.json'),wall_seconds=time.monotonic()-began)
        if lane=='primitive':
            proof=json.loads(subprocess.check_output([str(auditor),'--primitive-population',str(folder)],text=True))
            proof['sources']={f:sha(folder/f) for f in ['GROUP.tsv','PAIR_FRAMES.tsv','PRIMITIVE.tsv','LEADS.tsv','RESOLUTIONS.tsv','COUNTS.tsv']}
            save(folder/'EXACT_POPULATION.json',proof);entry['population']=proof
        pilots.append(entry);print('PILOT',json.dumps(entry),flush=True);return folder,entry
    # Complete only this named affine action, not all primitive degree16 groups.
    folder,asp=job('asp16_complete',16,'primitive',['--primitive-group','asp16'],attempts=1000000)
    r=asp['population']
    if r['independent_complete_population'] and r['no_unresolved'] and r['counts']['qualifying']==r['counts']['literal_negative']:
        exclusions['asp16_16']=dict(path=str(folder/'EXACT_POPULATION.json'),sha256=sha(folder/'EXACT_POPULATION.json'))
    # Larger action gets a bounded interval; never retire it from a prefix.
    agl_folder,_=job('agl16_prefix',16,'primitive',['--primitive-group','agl16'],attempts=5000)
    novelty['primitive_prefixes']={'agl16_16':dict(path=str(agl_folder/'EXACT_POPULATION.json'),sha256=sha(agl_folder/'EXACT_POPULATION.json'))}
    for n in range(11,19):job(f'triangle_{n}',n,'triangle')
    for n in [12,14,15,16,18,20,21,24,25]:job(f'returns_{n}',n,'returns')
    pairs=[];schedule=random.Random(20261013)
    for i in range(7):
        order=['cached','uncached'];schedule.shuffle(order);values={};folders={}
        for mode in order:
            folder=root/f'cache_{i}_{mode}';began=time.monotonic()
            command=[str(binary),'--d','16','--lane','primitive','--primitive-group','asp16','--frame-reuse','1',
                     '--primitive-cursor','150000','--attempts','50000','--seconds','60','--out',str(folder)]
            if mode=='uncached':command+=['--uncached-primitive-frames']
            subprocess.run(command,check=True,capture_output=True);wall=time.monotonic()-began
            with (folder/'COUNTS.tsv').open() as f:counts=next(csv.DictReader(f,delimiter='\t'))
            if int(counts['constructed'])!=50000:raise ValueError('benchmark did not finish fixed rank interval')
            seconds=float(counts.pop('seconds'));values[mode]=dict(producer_seconds=seconds,wall_seconds=wall,counts=counts);folders[mode]=folder
        for name in ['LEADS.tsv','RESOLUTIONS.tsv','KNOWN.tsv','GROUP.tsv','PAIR_FRAMES.tsv']:
            if sha(folders['cached']/name)!=sha(folders['uncached']/name):raise ValueError('benchmark changes mathematical output')
        if values['cached']['counts']!=values['uncached']['counts']:raise ValueError('benchmark population changed')
        pairs.append(dict(order=order,**values))
    benchmark=dict(pairs=pairs,identical_literal_outputs=True,rank_interval=[150000,200000],
        producer_speedup=statistics.median(x['uncached']['producer_seconds']/x['cached']['producer_seconds'] for x in pairs),
        whole_job_speedup=statistics.median(x['uncached']['wall_seconds']/x['cached']['wall_seconds'] for x in pairs),
        limitations='Matched affine workload only, includes identical proposal population. Not a discovery-yield or whole-search multiplier.')
    save(root/'BENCHMARK.json',benchmark);print('BENCHMARK',json.dumps({k:v for k,v in benchmark.items() if k!='pairs'}),flush=True)
    for name in ['known_states.tsv','controls.tsv','mutation_seeds.tsv','double_bases.tsv','profile_queue.json','INPUTS.json']:
        shutil.copy2(inputs/name,root/name)
    save(root/'NOVELTY_INPUT.json',novelty)
    save(root/'primitive_exclusions.json',exclusions)
    if any(sha(PROJECT/name)!=digest for name,digest in hashes.items()):raise ValueError('source changed during validation')
    result=dict(passed=True,design_passed=True,mechanism_passed=True,native_components_passed=True,novelty_passed=True,
        binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),controls=66,control_states=25192,
        known_bank_sha256=sha(root/'known_states.tsv'),mutation_seeds_sha256=sha(root/'mutation_seeds.tsv'),
        double_bases_sha256=sha(root/'double_bases.tsv'),profile_queue_sha256=sha(root/'profile_queue.json'),
        primitive_exclusions_sha256=sha(root/'primitive_exclusions.json'),novelty_input_sha256=sha(root/'NOVELTY_INPUT.json'),
        sources_sha256=sha(root/'SOURCE_HASHES.json'),tests=tests,pilots=pilots,benchmark_sha256=sha(root/'BENCHMARK.json'),
        limits='Only named finite-action intervals are censused; triangle/returns draws remain biased. No geometry, whole-degree or primitive nonexistence claim.')
    save(root/'RESULTS.json',result);print('VALIDATION='+str(root/'RESULTS.json'),flush=True)


if __name__=='__main__':main()
