"""Bind saved, independently checked partial pilot starts for later replay."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
from run_smooth_degree9 import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args();leads=[];sources=[]
    for folder in sorted(args.validation.glob('pilot_*')):
        if not (folder/'CHECKS.json').exists():raise ValueError('pilot lacks independent checks')
        check=json.loads((folder/'CHECKS.json').read_text());n=int(check['construction']['degree']);mass=0
        with (folder/'RESOLUTIONS.tsv').open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if row['status'] not in ['timeout','capped']:continue
                leads.append(dict(degree=n,seed=row['seed'],source=str(folder),candidate=row['candidate']))
                mass+=1
        if mass!=int(check['construction']['unresolved']):raise ValueError('unresolved mass')
        if mass:sources.append(dict(folder=str(folder),checks_sha256=sha(folder/'CHECKS.json'),resolutions_sha256=sha(folder/'RESOLUTIONS.tsv')))
    with args.out.open('x') as f:json.dump(leads,f,indent=2);f.write('\n')
    with args.out.with_suffix('.binding.json').open('x') as f:json.dump(dict(leads=len(leads),sha256=sha(args.out),sources=sources),f,indent=2);f.write('\n')
    print('SAVED_UNRESOLVED='+str(len(leads)))

if __name__=='__main__':main()
