"""Report exact failures of Proposition 1 for an annotated local catalog."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from classify_smooth_isotrivial import classify
from count_isotrivial_failures import count_failures


def summarize(data):
    unknowns=[e for e in data['components'] if e['isotrivial'] is None]
    block_patterns=Counter();candidate_patterns=Counter();component_patterns=Counter()
    component_presence=Counter();block_count=0;candidate_count=0
    ab_components=set();elliptic_components=set();elliptic_candidates=0
    elliptic_profiles={};equations=Counter()
    for e in unknowns:
        expected=classify(e)['isotriviality']['failed_block_systems']
        if e['isotriviality'].get('failed_block_systems')!=expected:
            raise ValueError('missing or stale exact failure diagnostics')
        flags=set()
        for b in expected:
            block_count+=1
            if not b['a_passing_entries']:
                block_patterns['a']=block_patterns['a']+1;flags.add('a')
            else:
                patterns=set()
                for c in b['candidate_failures']:
                    if not c['failed_conditions']:raise ValueError('unknown has a successful witness')
                    candidate_count+=1
                    key='+'.join(c['failed_conditions']);candidate_patterns[key]+=1
                    patterns.add(key);flags.update(c['failed_conditions'])
                    if c['failed_conditions']==['c']:
                        ab_components.add(e['id'])
                        q=e['block_quotients'][b['block_quotient_index']]
                        for f in c['ordinary_cycle_count_failures']:
                            equations[(b['block_size'],f['quotient_cycles'],
                                f['required_cycles'],f['actual_cycles'])]+=1
                        if q['genus']==1:
                            elliptic_components.add(e['id']);elliptic_candidates+=1
                            profile=(b['block_count'],b['block_size'],tuple(sorted(
                                f['ramification_excess'] for f in c['ordinary_cycle_count_failures'])))
                            elliptic_profiles.setdefault(profile,set()).add(e['id'])
                block_patterns[' / '.join(sorted(patterns))]+=1
        component_patterns['+'.join(sorted(flags)) or 'no recorded block system']+=1
        component_presence.update(flags)
    stats=dict(dataset_version=data['dataset_version'],unknown_components=len(unknowns),
        exclusive_component_failure_counts=count_failures(data['components'])['exclusive_failures'],
        unknown_states=sum(e['states'] for e in unknowns),recorded_block_systems=block_count,
        candidate_special_entries_passing_a=candidate_count,
        block_system_failure_patterns=dict(sorted(block_patterns.items())),
        candidate_failure_patterns_after_a=dict(sorted(candidate_patterns.items())),
        component_failure_patterns=dict(sorted(component_patterns.items())),
        component_reason_presence=dict(sorted(component_presence.items())),
        components_with_a_b_pass_c_fail=len(ab_components),
        components_with_elliptic_a_b_pass_c_fail=len(elliptic_components),
        elliptic_a_b_passing_candidates=elliptic_candidates,
        elliptic_ordinary_ramification_profiles=[dict(block_count=r,block_size=m,
            positive_ordinary_excesses=list(ex),components=len(ids))
            for (r,m,ex),ids in sorted(elliptic_profiles.items())],
        ordinary_cycle_failure_equations_after_a_b=[dict(block_size=m,
            quotient_cycles=q,required_cycles=required,actual_cycles=actual,
            ramification_excess=required-actual,occurrences=count)
            for (m,q,required,actual),count in sorted(equations.items())],
        counting='Exclusive component failure counts assign each component once to its furthest reachable stage and sum to the unknown count. '
            'Component patterns are disjoint and sum to the number of unknown components. '
            'Reason presence can overlap. Block-system and candidate counts are not component counts. '
            'These are representative-specific failures of the sufficient criterion, not proofs of non-isotriviality.')
    lines=['# Unknown entries: exact criterion failures','',
        f"Dataset version {data['dataset_version']}; {len(unknowns)} unknown components.",'',
        'The diagnostics below concern Proposition 1(a)-(c) for the displayed representative. '
        'Failure does not itself determine isotriviality.','',
        '## Statistics','', '```json',json.dumps(stats,indent=2),'```','']
    for i,e in enumerate(unknowns,1):
        lines.extend([f"## {i}. Degree {e['degree']}, genus {e['genus']}, {e['states']:,} states",''])
        lines.extend(f"- g{j+1} = {p}" for j,p in enumerate(e['tuple']))
        lines.append('')
        for b in e['isotriviality']['failed_block_systems']:
            original=e['block_quotients'][b['block_quotient_index']]
            lines.extend([f"### {b['block_count']} blocks of size {b['block_size']}",'',
                'Blocks: '+'; '.join('{'+', '.join(map(str,B))+'}' for B in original['blocks'])+'.',''])
            if not b['a_passing_entries']:
                lines.append('(a) fails for every entry: quotient cycle counts '
                    +str(b['quotient_cycle_counts'])+f"; an identity on blocks needs {b['block_count']} cycles.")
            for c in b['candidate_failures']:
                lines.append(f"Special entry g{c['special_entry']}: (a) passes.")
                lines.append(f"(b) {'fails' if 'b' in c['failed_conditions'] else 'passes'}: "
                    +f"{len(c['active_blocks'])} active blocks {c['active_blocks']}; exactly 1 required.")
                failures=c['ordinary_cycle_count_failures']
                if not failures:lines.append('(c) passes for every ordinary entry.')
                for f in failures:
                    lines.append(f"(c) fails for g{f['entry']}: {f['actual_cycles']} cycles != "
                        +f"{b['block_size']} * {f['quotient_cycles']} = {f['required_cycles']}; "
                        +f"ramification excess {f['required_cycles']} - {f['actual_cycles']} = {f['ramification_excess']}.")
            lines.append('')
    return stats,'\n'.join(lines).rstrip()+'\n'


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',type=Path,required=True);p.add_argument('--out-dir',type=Path,required=True)
    args=p.parse_args();raw=args.catalog.read_bytes();data=json.loads(raw)
    stats,report=summarize(data)
    stats['source_sha256']=hashlib.sha256(raw).hexdigest()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    for name,body in [('STATISTICS.json',json.dumps(stats,indent=2)+'\n'),('CALCULATIONS.md',report)]:
        with (args.out_dir/name).open('x') as stream:stream.write(body)
    print(json.dumps(stats,indent=2))
