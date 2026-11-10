"""Freeze saved positive banks, mutation controls, and a source-bound profile queue.

No search starts here. Only a finalized STOP predecessor is admissible.
Existing complete closure reports remain the proof; no capped table is imported.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
import shutil

from run_smooth_degree9 import sha
from run_forest_fast import predecessor_finished


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--predecessor',type=Path,required=True);p.add_argument('--theory',type=Path,required=True)
    p.add_argument('--prior-queue',type=Path,help='preserve completed profile entries, then append only new identities')
    p.add_argument('--out',type=Path,required=True);args=p.parse_args();args.out.mkdir()
    prior=args.predecessor.resolve();predecessor_finished(prior,True)
    entries=json.loads(args.theory.read_text())['components'];queue=[];seeds=[];double=[]
    for e in entries:
        folder=Path(e['raw_folder'])
        if sha(folder/'STATES.tsv')!=e['table_sha256']:raise ValueError('curated source changed')
        queue.append(dict(e,summary_sha256=sha(folder/'SUMMARY.tsv'),
            group_order_provenance='prior independent SymPy group profile'))
        if e['degree'] in [9,12]:double.append(f"{e['degree']}\t{e['representative']}\n")
        if e['degree']<=10:continue # Retired complete smooth censuses, not fresh exploration.
        labels=[None]*e['degree']
        for j,b in enumerate(e['block_partition']):
            for x in b:labels[x-1]=j
        if any(x is None for x in labels):raise ValueError('mutation partition mass')
        seeds.append(f"{e['degree']}\t{e['representative']}\t"+','.join(map(str,labels))+'\n')
    (args.out/'mutation_seeds.tsv').write_text(''.join(seeds))
    (args.out/'double_bases.tsv').write_text(''.join(double))
    if args.prior_queue:queue=json.loads(args.prior_queue.read_text())['components']
    new=json.loads((prior/'NEW_COMPONENTS.json').read_text());reports={}
    for e in new:
        job=prior/e['job']
        if not (job/'VERIFIED').exists():raise ValueError('unverified profile queue job')
        if str(job) not in reports:reports[str(job)]=json.loads((job/'CHECKS.json').read_text())
        report=reports[str(job)]
        candidates=[(c,i) for i in report['closed_identities'] for c in report['components']
                    if i['folder']==e['folder'] and i['key']==e['key'] and c['summary']==i['summary']]
        if len(candidates)!=1:raise ValueError('new closed identity binding')
        c,item=candidates[0]
        if c['summary']['status']!='closed_smooth' or not c['independent_six_move_closure']:raise ValueError('not complete closure')
        folder=job/e['folder']
        if sha(folder/'STATES.tsv')!=c['table_sha256']:raise ValueError('new closure changed')
        queue.append(dict(degree=e['degree'],states=int(e['summary']['states']),representative=e['summary']['seed'],
            component_identity_key=e['key'],
            group_order=int(e['summary']['group_order']),raw_folder=str(folder),table_sha256=c['table_sha256'],
            summary_sha256=sha(folder/'SUMMARY.tsv'),checks_sha256=sha(job/'CHECKS.json'),
            group_order_provenance='producer Schreier order; independent profile not yet curated'))
    (args.out/'profile_queue.json').write_text(json.dumps(dict(components=queue,
        initial_theory_sha256=sha(args.theory),new_index_sha256=sha(prior/'NEW_COMPONENTS.json'),
        prior=str(prior),curated=len(entries),new_closed=len(new)),indent=2)+'\n')
    bank=prior/'VERIFIED_BANKS';latest={};bindings=[]
    for path in sorted(bank.glob('degree_*_v*.json')):
        value=json.loads(path.read_text());n=value['degree'];latest[n]=(path,value)
    with (args.out/'known_states.tsv').open('w') as output:
        for n,(meta,value) in sorted(latest.items()):
            path=meta.with_suffix('.tsv')
            if sha(path)!=value['sha256']:raise ValueError('saved verified bank changed')
            provenance=value['provenance']
            if provenance.get('kind')=='independently verified six-move closures':
                source=Path(provenance['job'])
                if not (source/'VERIFIED').exists() or sha(source/'CHECKS.json')!=provenance['checks_sha256']:raise ValueError('bank proof binding')
            with path.open() as stream:shutil.copyfileobj(stream,output)
            bindings.append(dict(path=str(path),sha256=value['sha256'],rows=value['rows'],degree=n))
        with (prior/'KNOWN_STATES.tsv').open() as stream:
            for line in stream:
                if int(line.split('\t',1)[0]) not in latest:output.write(line)
    shutil.copy2(prior/'CONTROLS.tsv',args.out/'controls.tsv')
    result=dict(predecessor=str(prior),predecessor_finished_sha256=sha(prior/'FINISHED.json'),
        source_theory_sha256=sha(args.theory),new_index_sha256=sha(prior/'NEW_COMPONENTS.json'),
        known_bank_sha256=sha(args.out/'known_states.tsv'),known_bank_provenance=bindings,
        mutation_seeds_sha256=sha(args.out/'mutation_seeds.tsv'),mutation_seeds=len(seeds),
        double_bases_sha256=sha(args.out/'double_bases.tsv'),double_bases=len(double),
        profile_queue_sha256=sha(args.out/'profile_queue.json'),profile_components=len(queue),
        profile_start=len(queue)-len(new) if args.prior_queue else 0,
        limitations='Bank is RAM-bounded, not every newly closed state. Profile queue includes all closed identities; group orders/geometry await stated independent checks.')
    (args.out/'INPUTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
