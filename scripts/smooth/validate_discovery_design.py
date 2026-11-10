"""Bounded discovery-design checks; one CPU and a held campaign lock.

Paired logging measurement preserves candidate streams. Constructor pilots
are DIFFERENT biased populations; yield ratios are not runtime speedups.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import json
from pathlib import Path
import random
import statistics
import subprocess
import time

from audit_guided_native import audit_guided_native,known_keys
from benchmark_forest_second import interval
from run_smooth_degree9 import sha
from run_s18_linked_roots import state_key
from odd_order_spin import parse_seed

PROJECT=Path(__file__).resolve().parents[2]

def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--attempts',type=int,default=1000000)
    p.add_argument('--repeats',type=int,default=7)
    args=p.parse_args();root=args.out.resolve();root.mkdir()
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binary=args.binary.resolve();auditor=args.auditor.resolve()
    catalog=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components']
    controls=root/'controls.tsv';known=root/'known_states.tsv'
    controls.write_text(''.join(f"{c['degree']}\t{c['states']}\t{c['representative']}\n" for c in catalog))
    text=subprocess.check_output([str(binary),'--check-controls',str(controls),'--export-known',str(known)],text=True)
    (root/'CONTROLS.log').write_text(text);print(text.strip(),flush=True)
    expected=set()
    for c in catalog:
        folder=Path(c['raw_folder'])
        if sha(folder/'STATES.tsv')!=c['table_sha256'] or sha(folder/'SUMMARY.tsv')!=c['summary_sha256']:raise ValueError('catalog input changed')
        for row in rows(folder/'STATES.tsv'):
            expected.add((c['degree'],state_key(parse_seed(row['tuple'],c['degree']))))
    if known_keys(known)!=expected or len(expected)!=25192:raise ValueError('known export differs from independent catalog states')
    print('KNOWN_BANK exact independent set equality: 25192',flush=True)
    rng=random.Random(202610131);timings=[]
    for repeat in range(args.repeats):
        order=['dense','compact'];rng.shuffle(order);pair={};checks={};folders={}
        for variant in order:
            folder=root/f'logging_{repeat}_{variant}';folders[variant]=folder
            cmd=[binary,'--d','11','--lane','forest','--frame-reuse','1','--seed',str(11000+repeat),
                 '--seconds','120','--attempts',str(args.attempts),'--out',folder]
            if variant=='compact':cmd+=['--compact-frames']
            began=time.monotonic();subprocess.run(list(map(str,cmd)),check=True,stdout=subprocess.DEVNULL)
            producer=time.monotonic()-began;began=time.monotonic();check=audit_guided_native(folder,auditor)
            total=producer+time.monotonic()-began;checks[variant]=check
            (folder/'CHECKS.json').write_text(json.dumps(check,indent=2)+'\n')
            pair[variant]=dict(producer=producer,total=total,frame_bytes=(folder/'FRAMES.tsv').stat().st_size)
        if sha(folders['dense']/'LEADS.tsv')!=sha(folders['compact']/'LEADS.tsv') or sha(folders['dense']/'RESOLUTIONS.tsv')!=sha(folders['compact']/'RESOLUTIONS.tsv'):raise ValueError('logging changed literal outcomes')
        for k,v in checks['dense']['construction'].items():
            if k!='seconds' and v!=checks['compact']['construction'][k]:raise ValueError('logging changed '+k)
        if int(checks['dense']['construction']['constructed'])!=args.attempts:raise ValueError('benchmark hit time cap')
        sampled=rows(folders['compact']/'FRAMES.tsv');cursor=0
        with (folders['dense']/'FRAMES.tsv').open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if (int(row['frame'])-1)%4096==0:
                    if {k:sampled[cursor][k] for k in row}!=row:raise ValueError('sparse frame mismatch')
                    cursor+=1
        if cursor!=len(sampled):raise ValueError('sparse frame mass')
        if rows(folders['dense']/'DESIGN.tsv')[0]['rng_final']!=rows(folders['compact']/'DESIGN.tsv')[0]['rng_final']:raise ValueError('RNG changed')
        timings.append(pair);print('LOGGING',repeat,json.dumps(pair),flush=True)
    pilot=[];distinct={}
    for n,lane in [(11,'forest'),(11,'joint'),(14,'forest'),(14,'joint'),
                   (15,'solved'),(18,'solved'),(20,'solved'),(21,'solved'),(24,'solved'),(25,'solved')]:
        folder=root/f'pilot_{n}_{lane}';began=time.monotonic()
        cmd=[binary,'--d',str(n),'--lane',lane,'--frame-reuse','1' if lane=='solved' else '256',
             '--seed',str(202610130+n),'--seconds','8','--attempts','10000' if lane!='solved' else '300',
             '--cap','10000','--closure-seconds','1','--compact-frames','--export-negative-cache',
             '--known-bank',known,'--out',folder]
        subprocess.run(list(map(str,cmd)),check=True,stdout=subprocess.DEVNULL)
        check=audit_guided_native(folder,auditor)
        (folder/'CHECKS.json').write_text(json.dumps(check,indent=2)+'\n')
        for item in check['closed_identities']:
            key=json.dumps([item['degree'],item['key']],separators=(',',':'))
            distinct.setdefault(key,dict(**item,job=str(folder)))
        pilot.append(dict(degree=n,lane=lane,counts=check['construction'],known_aliases=check['known_aliases'],
                          verified_closed=len(check['closed_identities']),total_seconds=time.monotonic()-began))
        print('PILOT',json.dumps(pilot[-1]),flush=True)
    # Reuse an exported, independently checked proof bank on another seed.
    source=root/'pilot_11_forest'/'NEGATIVE_CACHE.tsv';folder=root/'bank_reuse'
    subprocess.run(list(map(str,[binary,'--d','11','--seed','202610148','--seconds','8','--attempts','10000',
                    '--negative-bank',source,'--export-negative-cache','--compact-frames','--out',folder])),check=True,stdout=subprocess.DEVNULL)
    reuse=audit_guided_native(folder,auditor);(folder/'CHECKS.json').write_text(json.dumps(reuse,indent=2)+'\n')
    result=dict(passed=True,design_passed=True,controls=66,control_states=25192,known_bank_sha256=sha(known),
                binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),
                logging=dict(repeats=args.repeats,attempts=args.attempts,timings=timings,
                    pipeline=interval([t['dense']['total']/t['compact']['total'] for t in timings]),
                    byte_reduction=statistics.median(t['dense']['frame_bytes']/t['compact']['frame_bytes'] for t in timings),
                    identical_promoted_seeds_words_counters_rng=True),
                pilots=pilot,distinct_new_components=list(distinct.values()),bank_reuse_passed=True,
                limitations='Bounded biased pilots, not census or discovery-rate promise. Exact new S_d components, not weak failures, full A/S or automatic geometry classification.')
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print('RESULTS='+str(root/'RESULTS.json'),flush=True)
    lock.close()

if __name__=='__main__':main()
