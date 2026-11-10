"""Bounded same-input second-pass timings, including independent verification."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import time

from audit_guided_native import audit_guided_native
from run_smooth_frontier import audit_guided
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]


def interval(ratios):
    rng=random.Random(20261012)
    draws=sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(10000))
    return dict(paired_median_speedup=statistics.median(ratios),
                bootstrap_median_interval95=[draws[249],draws[9749]],paired_ratios=ratios)


def compare(old,new):
    if old['construction'] != new['construction'] or old['components'] != new['components']:
        raise ValueError('independent audits differ')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old',type=Path,required=True);p.add_argument('--new',type=Path,required=True)
    p.add_argument('--auditor',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source-job',type=Path,required=True);p.add_argument('--attempts',type=int,default=100000)
    p.add_argument('--repeats',type=int,default=7);args=p.parse_args()
    root=args.out.resolve();root.mkdir();old=args.old.resolve();new=args.new.resolve();auditor=args.auditor.resolve()
    catalog=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components']
    controls=root/'controls.tsv';controls.write_text(''.join(f"{c['degree']}\t{c['states']}\t{c['representative']}\n" for c in catalog))
    text=subprocess.check_output([str(new),'--check-controls',str(controls)],text=True)
    (root/'CONTROLS.log').write_text(text);print('CONTROLS',text.strip(),flush=True)
    # Equal prefix of immutable negative records, not unequal-time jobs.
    sample=root/'audit_sample';sample.mkdir()
    shutil.copy2(args.source_job/'COUNTS.tsv',sample/'COUNTS.tsv')
    for name in ['LEADS.tsv','RESOLUTIONS.tsv']:
        with (args.source_job/name).open() as source,(sample/name).open('w') as target:
            for i,line in enumerate(source):
                if i>=20000+int(name=='RESOLUTIONS.tsv'):break
                target.write(line)
    audit_times=[]
    for repeat in range(args.repeats):
        pair={};reports={}
        for variant in (['python','native'] if repeat%2==0 else ['native','python']):
            began=time.monotonic()
            reports[variant]=audit_guided(sample) if variant=='python' else audit_guided_native(sample,auditor)
            pair[variant]=time.monotonic()-began
        compare(reports['python'],reports['native']);audit_times.append(pair)
    audits=interval([r['python']/r['native'] for r in audit_times]);audits['timings']=audit_times
    print('AUDIT_TIMING',json.dumps(audits),flush=True)
    # Unprofiled whole completed job. Profiling overhead is not a speed claim.
    full={};full_reports={}
    for variant in ['native','python']:
        began=time.monotonic()
        full_reports[variant]=audit_guided(args.source_job) if variant=='python' else audit_guided_native(args.source_job,auditor)
        full[variant]=time.monotonic()-began
    compare(full_reports['python'],full_reports['native'])
    print('FULL_AUDIT',json.dumps(full),flush=True)
    stages=[];streams=[]
    for n,lane in [(11,'forest'),(12,'forest'),(13,'forest'),(14,'forest'),(12,'blocks')]:
        timings=[]
        for repeat in range(args.repeats):
            pair={};reports={};folders={}
            for variant in (['old','new'] if repeat%2==0 else ['new','old']):
                folder=root/f'{lane}_{n}_{repeat}_{variant}';folders[variant]=folder
                command=[str(old if variant=='old' else new),'--d',str(n),'--lane',lane,'--frame-reuse','256',
                         '--seconds','120','--attempts',str(args.attempts),'--seed',str(7100+repeat),'--out',str(folder)]
                began=time.monotonic();subprocess.run(command,check=True,stdout=subprocess.DEVNULL)
                producer=time.monotonic()-began
                began=time.monotonic();reports[variant]=audit_guided(folder) if variant=='old' else audit_guided_native(folder,auditor)
                audit=time.monotonic()-began
                pair[variant]=dict(producer=producer,audit=audit,total=producer+audit)
                (folder/'CHECKS.json').write_text(json.dumps(reports[variant],indent=2)+'\n')
            fields=['constructed','scope','transitive','passport','short_reject','duplicate','closed','unresolved']
            if any(reports['old']['construction'][f]!=reports['new']['construction'][f] for f in fields):
                raise ValueError('candidate stream counters differ')
            for name in ['LEADS.tsv','RESOLUTIONS.tsv','FRAMES.tsv']:
                if sha(folders['old']/name)!=sha(folders['new']/name):raise ValueError('candidate/witness/frame stream differs: '+name)
            timings.append(pair)
        result=dict(degree=n,lane=lane,timings=timings,
                    producer=interval([v['old']['producer']/v['new']['producer'] for v in timings]),
                    pipeline=interval([v['old']['total']/v['new']['total'] for v in timings]))
        stages.append(result);print('PIPELINE',n,lane,json.dumps(result),flush=True)
        streams.append(dict(degree=n,lane=lane,identical_counts_leads_words_frames=True))
    result=dict(passed=True,controls=len(catalog),control_states=sum(c['states'] for c in catalog),
                binary_sha256=sha(new),audit_binary_sha256=sha(auditor),old_binary_sha256=sha(old),
                audit_sample=audits,full_job_audit=full,full_job_records=full_reports['native']['independent_native']['records'],
                stages=stages,streams=streams,repeats=args.repeats,
                limitations='Shared busy Mac; equal finite inputs, startup/I/O included; no discovery-rate guarantee. Rejection-only full-job audit; positive/partial components use general Python audit.')
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('RESULTS='+str(root/'RESULTS.json'),flush=True)


if __name__=='__main__':main()
