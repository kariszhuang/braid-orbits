#!/usr/bin/env python3
"""Read only finalized VERIFIED work; produce small diagnostic/benchmark fixtures.

Proposal occurrences, late labelled starts and braid components are different
grains. No failed sample or zero-yield lane is an exclusion of a family.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
from run_genus_three import sha,write

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();finish=json.loads((a.root/'FINISHED.json').read_text())
    if finish['status']!='stopped':raise ValueError('only finalized stopped roots')
    a.out.mkdir();raw=(a.root/'COMPLETED.jsonl').read_bytes();rows=[json.loads(l)for l in raw.splitlines()]
    lanes=defaultdict(Counter);hist=defaultdict(Counter);bound=hashlib.sha256();fixture=[];degree_samples=Counter()
    for row in rows:
        folder=Path(row['folder']);assert folder.parent==a.root/'jobs' and (folder/'VERIFIED').exists()
        lane=row['counts']['lane'];c=lanes[lane];c['jobs']+=1
        for k,v in row['counts'].items():
            if k not in ['degree','lane','seed']:c[k]+=float(v)if k=='seconds'else int(v)
        data=(folder/'STAGES.tsv').read_bytes();bound.update(folder.name.encode()+b'\0'+data)
        stage=next(csv.DictReader(data.decode().splitlines(),delimiter='\t'))
        for k,v in stage.items():
            if k!='rng_final':c[k]+=int(v)
        with (folder/'RESOLUTIONS.tsv').open() as stream:
            for result in csv.DictReader(stream,delimiter='\t'):
                assert result['status']=='nonsmooth'
                word=tuple(map(int,result['word'].split(',')));assert word
                hist[lane][len(word)]+=1
                n=int(row['counts']['degree'])
                if degree_samples[n]<10:
                    fixture.append(f"{n}\t{result['seed']}\n");degree_samples[n]+=1
    assert len(rows)==finish['completed_jobs']
    for name,key in [('constructed','constructed'),('scope','scope'),('transitive','transitive'),('closed','closed_occurrences'),('unresolved','unresolved')]:
        assert sum(c[name]for c in lanes.values())==finish['counts'][key]
    assert sum(sum(h.values())for h in hist.values())==sum(c['closure_calls']for c in lanes.values())
    total_seconds=sum(c['seconds']for c in lanes.values())
    for c in lanes.values():
        c['job_time_share']=c['seconds']/total_seconds
        c['late_calls_per_job_second']=c['closure_calls']/c['seconds']
    (a.out/'TRIE_FIXTURE.tsv').write_text(''.join(fixture))
    result=dict(root=str(a.root),finished=finish,completed_sha256=hashlib.sha256(raw).hexdigest(),
        stages_stream_sha256=bound.hexdigest(),lanes=lanes,late_negative_word_lengths=hist,
        late_calls=sum(c['closure_calls']for c in lanes.values()),fixture_rows=len(fixture),
        fixture_sha256=sha(a.out/'TRIE_FIXTURE.tsv'),
        scope='VERIFIED proposal occurrences,not distinct states/orbits; transitive counter occurs AFTER short/deep rejection; job-seconds shares approximate CPU allocation,not measured CPU time')
    write(a.out/'REVIEW.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
