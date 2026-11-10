"""Freeze resolved-frontier and profile-cursor inputs for the authorized mission."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
from run_smooth_degree9 import sha,save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--replay',type=Path,required=True)
    p.add_argument('--lattice',type=Path,required=True);args=p.parse_args()
    report=json.loads((args.replay/'RESULTS.json').read_text())
    if not report['passed'] or report['original_starts_certified_negative']!=139 or report['original_occurrences']!=139:raise ValueError('frontier not completely negative')
    frontier=[];mass=0
    for proof in report['original_witness_audits']:
        path=Path(proof['path']);r=json.loads(path.read_text())
        if sha(path)!=proof['sha256'] or r['independent_native']['closed'] or r['independent_native']['unresolved']:raise ValueError('negative proof changed')
        for name,digest in r['sources'].items():
            if sha(path.parent/name)!=digest:raise ValueError('negative source changed')
        n=int(r['construction']['degree']);mass+=r['independent_native']['nonsmooth']
        frontier.append(dict(folder=str(path.parent),degree=n,checks_sha256=proof['sha256']))
    if mass!=139:raise ValueError('negative frontier mass')
    queue=args.inputs/'profile_queue.json';items=json.loads(queue.read_text())['components']
    lattice=json.loads(args.lattice.read_text())
    if lattice['queue_sha256']!=sha(queue) or lattice['start']!=0 or not lattice['complete'] or len(lattice['components'])!=len(items):raise ValueError('lattice queue not bound')
    if not all(e['block_lattice']['complete'] for e in lattice['components']):raise ValueError('lattice includes unresolved caps')
    out=args.inputs/'NOVELTY_INPUT.json'
    if out.exists():raise ValueError('refusing to overwrite frozen novelty input')
    save(out,dict(profile_start=5595,profile_total=len(items),profile_queue_sha256=sha(queue),
        negative_frontier=frontier,replay_results_sha256=sha(args.replay/'RESULTS.json'),
        complete_lattice=dict(path=str(args.lattice.resolve()),sha256=sha(args.lattice)),
        resources=dict(cpu_slots=6,gpu=False,unbounded=True,disk_guard_gib=30),
        priorities=['triangle index budgets','multi-active-block prescribed returns','proper affine primitive actions16','all block systems of new closures'],
        retired='Complete degree9/10 strata and exact named primitive actions are not rerun. Negative original139 starts stay retired.',
        assumptions='Biased discovery draws. A smooth seed is not a smooth orbit. Only complete independent six-move closure certifies a result.'))
    print('NOVELTY_INPUT='+str(out),flush=True)


if __name__=='__main__':main()
