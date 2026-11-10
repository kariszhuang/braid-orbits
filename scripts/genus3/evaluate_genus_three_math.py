"""Paired same-stream regression plus bounded new-population pilots.

Continuous production may contend for CPU. Timings are observations, not a
discovery multiplier. Every retained negative/exception has an independent audit.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
import random
import statistics
import subprocess
import time
from audit_guided_native import audit_guided_native
from run_genus_three import sha,write
from verify_genus_three_exclusions import audit_exclusions


def run(binary,auditor,words,out,n,lane,seed,attempts=0,seconds=3):
    command=[str(binary),'--d',str(n),'--lane',lane,'--seed',str(seed),'--seconds',str(seconds),
        '--cap','1000000','--closure-seconds','10','--out',str(out),'--words',str(words),'--trace-early','128']
    if attempts:command.extend(['--attempts',str(attempts)])
    start=time.perf_counter();proc=subprocess.run(command,capture_output=True,text=True,check=True)
    wall=time.perf_counter()-start
    with (out/'COUNTS.tsv').open()as f:counts=next(csv.DictReader(f,delimiter='\t'))
    with (out/'STAGES.tsv').open()as f:stages=next(csv.DictReader(f,delimiter='\t'))
    checks=audit_guided_native(out,auditor,native_components=True)
    checks['isotrivial_exclusions']=audit_exclusions(out/'ISOTRIVIAL.tsv',int(stages['isotrivial_outside_target']))
    write(out/'CHECKS.json',checks)
    return dict(wall=wall,counts=counts,stages=stages,
        journals={name:sha(out/name)for name in ['LEADS.tsv','RESOLUTIONS.tsv','ISOTRIVIAL.tsv']},
        audit=checks['independent_native'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['binary','baseline','auditor','words','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--compact-binary',type=Path)
    a=p.parse_args();a.out.mkdir();pairs=[];order=random.Random(272103)
    for n,lane in [(18,'collision'),(24,'full-fiber'),(33,'collision')]:
        for rep in range(3):
            names=['before','after'];order.shuffle(names);data={}
            for name in names:
                binary=a.baseline if name=='before'else (a.compact_binary if n<=25 and a.compact_binary else a.binary)
                data[name]=run(binary,a.auditor,a.words,a.out/f'{n}_{rep}_{name}',n,lane,92311+rep,
                    attempts=1000000,seconds=60)
            for key in ['journals','stages']:assert data['before'][key]==data['after'][key],(n,rep,key)
            cb={k:v for k,v in data['before']['counts'].items()if k!='seconds'}
            ca={k:v for k,v in data['after']['counts'].items()if k!='seconds'}
            assert cb==ca,(n,rep,'counters')
            pairs.append(dict(degree=n,lane=lane,replicate=rep,data=data,
                speed_ratio=data['before']['wall']/data['after']['wall']))
    pilots=[]
    for n,lane in [(36,'collision'),(48,'collision'),(60,'collision'),(63,'collision'),
        (36,'partial-fiber'),(42,'partial-fiber'),(63,'partial-fiber'),
        (24,'general-fiber'),(36,'general-fiber'),(48,'general-fiber'),(60,'general-fiber'),
        (28,'four-quotient'),(35,'four-quotient'),(52,'four-quotient'),(56,'four-quotient'),(60,'four-quotient'),(63,'four-quotient')]:
        result=run(a.binary,a.auditor,a.words,a.out/f'pilot_{n}_{lane}',n,lane,44011,seconds=3)
        pilots.append(dict(degree=n,lane=lane,result=result));print(n,lane,result['counts']['transitive'],flush=True)
    write(a.out/'RESULT.json',dict(passed=True,binary_sha256=sha(a.binary),baseline_sha256=sha(a.baseline),
        compact_binary_sha256=sha(a.compact_binary) if a.compact_binary else None,
        audit_binary_sha256=sha(a.auditor),words_sha256=sha(a.words),same_stream_pairs=pairs,pilots=pilots,
        median_wall_speed_ratios={str(n):statistics.median(v['speed_ratio']for v in pairs if v['degree']==n)for n in [18,24,33]},
        scope='same-stream assertions apply only to paired unchanged populations; novel pilots sampled, no degree exhaustion or discovery-speed guarantee'))


if __name__=='__main__':main()
