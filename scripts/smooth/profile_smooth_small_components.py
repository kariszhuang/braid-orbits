#!/usr/bin/env python3
"""Independent exact group orders and passport descriptions of smooth closures."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import json
from pathlib import Path
import sympy
from sympy.combinatorics import Permutation, PermutationGroup

from odd_order_spin import cycle_lengths, parse_seed
from run_smooth_mission import audit_component, save, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(); results = []
    for folder in sorted(args.root.glob('degree_*/component_*')):
        audit = audit_component(folder); summary = audit['summary']; n = int(summary['degree'])
        t = parse_seed(summary['seed'], n)
        order = int(PermutationGroup([Permutation(list(q)) for q in t]).order())
        if order != int(summary['group_order']):
            raise ValueError('independent group order')
        types = sorted((sorted(cycle_lengths(q), reverse=True) for q in t))
        results.append(dict(folder=str(folder.relative_to(args.root)), degree=n, genus=int(summary['genus']),
                            states=int(summary['states']), group_order=order, full_group=summary['full_group']=='true',
                            passport=types, has_full_cycle=any(lengths == [n] for lengths in types), seed=summary['seed'],
                            independent_six_move_closure=audit['independent_six_move_closure'],
                            summary_sha256=sha(folder/'SUMMARY.tsv'), states_sha256=sha(folder/'STATES.tsv')))
    answer = dict(input_root=str(args.root), verifier_sha256=sha(Path(__file__)), sympy_version=sympy.__version__,
                  components=results, by_degree=dict(Counter(r['degree'] for r in results)),
                  by_genus=dict(Counter(r['genus'] for r in results)), all_have_full_cycle=all(r['has_full_cycle'] for r in results))
    save(args.out, answer)
    print(json.dumps({k: v for k, v in answer.items() if k!='components'}, indent=2))


if __name__=='__main__':
    main()
