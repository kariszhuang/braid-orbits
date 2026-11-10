"""Count failures of Nitiwit's sufficient criterion; reads JSON, writes nothing."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
from classify_smooth_isotrivial import classify


def count_failures(components):
    totals=dict(unknown_components=0,block_systems=0,a_failed=0,
        a_passing_candidates=0,b_failed=0,ab_passing_candidates=0,c_failed=0)
    affected={letter:set() for letter in 'abc'}
    exclusive={letter:0 for letter in 'abc'}
    for index,e in enumerate(components):
        result=classify(e)
        if result['isotrivial'] is not None:continue
        totals['unknown_components']+=1
        has_a=False;has_ab=False
        for b in result['isotriviality']['failed_block_systems']:
            totals['block_systems']+=1
            if not b['a_passing_entries']:
                totals['a_failed']+=1;affected['a'].add(index)
                continue
            has_a=True
            for c in b['candidate_failures']:
                totals['a_passing_candidates']+=1
                if 'b' in c['failed_conditions']:
                    totals['b_failed']+=1;affected['b'].add(index)
                    continue
                totals['ab_passing_candidates']+=1
                has_ab=True
                if 'c' in c['failed_conditions']:
                    totals['c_failed']+=1;affected['c'].add(index)
                else:raise ValueError('unknown component has a successful witness')
        exclusive['c' if has_ab else 'b' if has_a else 'a']+=1
    totals['affected_components']={k:len(v) for k,v in affected.items()}
    totals['exclusive_failures']=exclusive
    if sum(exclusive.values())!=totals['unknown_components']:
        raise ValueError('exclusive categories do not partition the unknown components')
    return totals


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('catalog',type=Path,help='path to all-smooth.json')
    p.add_argument('--details',action='store_true',help='also show overlapping partition/candidate counts')
    args=p.parse_args();t=count_failures(json.loads(args.catalog.read_text())['components'])
    print(f"Unknown components: {t['unknown_components']}")
    print('Each component is counted once, at the furthest stage it reaches:')
    print(f"(a): {t['exclusive_failures']['a']} — no choice of blocks/special entry passes (a)")
    print(f"(b): {t['exclusive_failures']['b']} — some choice passes (a), but none passes both (a),(b)")
    print(f"(c): {t['exclusive_failures']['c']} — some choice passes (a),(b), but none passes all three")
    print(f"Total: {sum(t['exclusive_failures'].values())}")
    if args.details:
        print('\nAdditional checks (different units; do not add these):')
        print(f"(a) fails: {t['a_failed']} / {t['block_systems']} block systems have no entry fixing all blocks")
        print(f"(b) fails: {t['b_failed']} / {t['a_passing_candidates']} special-entry candidates passing (a)")
        print(f"(c) fails: {t['c_failed']} / {t['ab_passing_candidates']} candidates passing both (a) and (b)")
        print('Component presence (can overlap): '+str(t['affected_components']))
