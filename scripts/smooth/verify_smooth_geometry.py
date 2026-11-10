#!/usr/bin/env python3
"""Independent exact hypotheses for the fixed-quotient geometric proof.

No producer construction, block-recovery, or permutation helpers are imported.
The numerical audit is not a formal verification of the geometric theorem.
Default is plan-only. --run freezes inputs and checks every closed state.
"""
import argparse
import collections
import csv
import datetime as dt
import fcntl
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

PROJECT = Path(__file__).resolve().parents[2]
CATALOG = PROJECT / 'evidence/smooth_special_20261012/CATALOG.json'
PRIOR_CHECKS = PROJECT / 'evidence/smooth_special_20261012/CHECKS.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def parse_tuple(text, n):
    require(text.startswith('[') and text.endswith(']'), 'tuple brackets')
    entries = text[1:-1].split(',')
    require(len(entries) == 4, 'four ordered entries required')
    result = []
    for entry in entries:
        entry = entry.strip()
        require(re.sub(r'\([^()]*\)', '', entry).strip() == '', 'invalid cycle syntax')
        p = list(range(n))
        used = set()
        for body in re.findall(r'\(([^()]*)\)', entry):
            cycle = [int(x)-1 for x in body.split()]
            require(all(0 <= x < n for x in cycle), 'point outside degree')
            require(len(cycle) == len(set(cycle)), 'repeated point in cycle')
            require(not used.intersection(cycle), 'overlapping cycles')
            used.update(cycle)
            for x, y in zip(cycle, cycle[1:] + cycle[:1]):
                p[x] = y
        result.append(tuple(p))
    return tuple(result)


def compose(a, b):
    return tuple(a[b[x]] for x in range(len(a)))


def inverse(p):
    q = [0]*len(p)
    for x, y in enumerate(p):
        q[y] = x
    return tuple(q)


def cycles(p):
    seen, out = set(), []
    for x in range(len(p)):
        if x in seen:
            continue
        c, y = [], x
        while y not in seen:
            seen.add(y)
            c.append(y)
            y = p[y]
        out.append(c)
    return out


def orbits(permutations):
    n = len(permutations[0])
    remaining, out = set(range(n)), []
    while remaining:
        x = min(remaining)
        found, queue = {x}, [x]
        for y in queue:
            for p in permutations:
                if p[y] not in found:
                    found.add(p[y])
                    queue.append(p[y])
        remaining.difference_update(found)
        out.append(sorted(found))
    return out


def hurwitz(t, m):
    i = m % 3
    a, b = t[i:i+2]
    out = list(t)
    if m < 3:
        out[i:i+2] = [compose(compose(a, b), inverse(a)), a]
    else:
        out[i:i+2] = [b, compose(compose(inverse(b), a), b)]
    return tuple(out)


def forced_partition(t, initial):
    """Least G-invariant equivalence relation containing the required joins."""
    n = len(t[0])
    parent = list(range(n))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def merge(x, y):
        x, y = root(x), root(y)
        if x == y:
            return False
        parent[max(x, y)] = min(x, y)
        return True

    for part in initial:
        for x in part:
            merge(part[0], x)
    changed = True
    while changed:
        changed = False
        for p in t:
            for x in range(n):
                changed = merge(p[x], p[root(x)]) or changed
    parts = collections.defaultdict(list)
    for x in range(n):
        parts[root(x)].append(x)
    return tuple(tuple(v) for _, v in sorted(parts.items()))


def recover_blocks(t, support, r, k):
    # For all ten higher-genus cases support already IS the whole block.
    # For smaller local supports enumerate invariant coarsenings, not point
    # subsets: exact equivalence closure and size pruning are sufficient.
    base = forced_partition(t, [sorted(support)])
    stack, seen = [base], set()
    anchor = min(support)
    while stack:
        partition = stack.pop()
        if partition in seen:
            continue
        seen.add(partition)
        require(len(seen) <= 10000, 'block coarsening cap: not verified')
        part = next(x for x in partition if anchor in x)
        if len(part) > k:
            continue
        if len(part) == k:
            require(len(partition) == r and all(len(x) == k for x in partition),
                    'transitive invariant partition has wrong sizes')
            b = next(p for p in t if all(p[x] != x for x in range(r*k)))
            blocks, current = [], set(part)
            for _ in range(r):
                blocks.append(sorted(current))
                current = {b[x] for x in current}
            if len({tuple(x) for x in blocks}) == r:
                require(current == set(blocks[0]), 'block return failed')
                return blocks
            continue
        for other in partition:
            if other != part and len(part) + len(other) <= k:
                stack.append(forced_partition(t, [list(part + other)]))
    raise ValueError('no required invariant cyclic block partition')


def analyze(t, r, k):
    n = r*k
    identity = tuple(range(n))
    require(len(t) == 4 and all(sorted(p) == list(identity) for p in t), 'invalid permutations')
    require(compose(compose(compose(t[0], t[1]), t[2]), t[3]) == identity, 'product is not one')
    require(len(orbits(t)) == 1, 'intransitive total group')
    supports = [{x for x in range(n) if p[x] != x} for p in t]
    small = [i for i, s in enumerate(supports) if 0 < len(s) <= k]
    require(len(small) == 1, 'unique nonidentity local entry required')
    s = small[0]
    require(all(len(supports[i]) == n for i in range(4) if i != s), 'large support is not full')
    blocks = recover_blocks(t, supports[s], r, k)
    membership = {x: j for j, block in enumerate(blocks) for x in block}
    require(len(membership) == n, 'blocks do not partition points')
    quotients, shifts, relative = [], [], []
    for i, p in enumerate(t):
        q = []
        for block in blocks:
            image = {membership[p[x]] for x in block}
            require(len(image) == 1, 'entry does not preserve block system')
            q.append(next(iter(image)))
        shift = q[0]
        require(all(q[j] == (j+shift) % r for j in range(r)), 'quotient is not a translation')
        if i == s:
            require(shift == 0 and supports[i].issubset(blocks[0]), 'local entry not localized')
        else:
            require(math.gcd(shift, r) == 1, 'large quotient translation not a unit')
            require(sorted(map(len, cycles(p))) == [r]*k, 'large cycle type not r^k')
        quotients.append(q)
        shifts.append(shift)
        # Each quotient cycle corresponds to ONE point of E above that branch
        # value. Its inverse image points have relative indices len(cycle)/eE.
        branch_fibers = []
        point_cycles = cycles(p)
        for qc in cycles(q):
            qcset = set(qc)
            relevant = [c for c in point_cycles if membership[c[0]] in qcset]
            require(all(len(c) % len(qc) == 0 for c in relevant), 'nonintegral relative ramification')
            e = sorted(len(c)//len(qc) for c in relevant)
            require(sum(e) == k, 'relative fiber degree mismatch')
            branch_fibers.append(dict(quotient_blocks=qc, relative_indices=e,
                                      ramified=any(x > 1 for x in e)))
        relative.append(branch_fibers)
    require(sum(shifts) % r == 0, 'quotient product not one')
    branch_counts = [sum(x['ramified'] for x in fibers) for fibers in relative]
    require(branch_counts == [int(i == s) for i in range(4)], 'not exactly one relative branch value')
    indices = [n-len(cycles(p)) for p in t]
    require(sum(indices) % 2 == 0, 'source genus nonintegral')
    genus = 1-n+sum(indices)//2
    quotient_indices = [r-len(cycles(q)) for q in quotients]
    require(sum(quotient_indices) % 2 == 0, 'quotient genus nonintegral')
    quotient_genus = 1-r+sum(quotient_indices)//2
    relative_index = sum(e-1 for fibers in relative for fiber in fibers for e in fiber['relative_indices'])
    require(relative_index == indices[s], 'relative branching index mismatch')
    require(2*genus-2 == k*(2*quotient_genus-2)+relative_index, 'relative Riemann-Hurwitz failed')
    require(quotient_genus == (r-1)//2 and r % 2 == 1, 'wrong cyclic three-point genus')
    left, right = len(orbits(t[:2])), len(orbits(t[2:]))
    product_index = n-len(cycles(compose(t[0], t[1])))
    smooth = product_index in [indices[0]+indices[1], indices[2]+indices[3]]
    defect = n-product_index-left-right+1
    require(smooth and defect == 0, 'state not smooth/defect zero')
    large = [i for i in range(4) if i != s]
    u, v = shifts[large[0]], shifts[large[1]]
    require(math.gcd(u+v, r) == 1, 'Kummer infinity not fully branched')
    return dict(degree=n, blocks=r, block_size=k, small_entry=s+1,
                block_partition=[[x+1 for x in b] for b in blocks],
                quotient_permutations=quotients, quotient_shifts=shifts,
                quotient_genus=quotient_genus, source_genus=genus,
                indices=indices, local_cycle_lengths=sorted(len(c) for c in cycles(t[s]) if len(c)>1),
                relative_ramification=relative, relative_branch_value_count=sum(branch_counts),
                relative_branching_index=relative_index, product_index=product_index,
                left_orbits=left, right_orbits=right, smooth=smooth, defect=defect,
                kummer=dict(r=r, u=u, v=v, equation=f'y^{r}=x^{u}*(x-1)^{v}',
                            interpretation='smooth projective normalization; fixed quotient, NOT an equation of C'),
                tangent_pullback_degree=k*(2-2*quotient_genus),
                geometry='nonisotrivial_higher_genus' if quotient_genus >= 2 else 'isotrivial_elliptic')


def audit(root):
    start = time.monotonic()
    plan = json.loads((root/'PLAN.json').read_text())
    for relative, digest in plan['snapshot_sha256'].items():
        require(sha(root/relative) == digest, f'snapshot hash: {relative}')
    catalog = json.loads((root/'input/CATALOG.json').read_text())
    prior = json.loads((root/'input/PRIOR_CHECKS.json').read_text())
    require(prior['independent_six_move_closure'] and prior['independent_group_orders'], 'prior closure/order checks absent')
    require(catalog['component_count'] == len(catalog['components']) == 52, 'scope must include all52')
    raw = Path(catalog['raw_root'])
    output, totals, ids = [], collections.Counter(), set()
    for component in catalog['components']:
        folder = raw/component['folder']
        table = folder/'STATES.tsv'
        digest = sha(table)
        require(digest == component['table_sha256'] == plan['state_table_sha256'][component['folder']], 'state table hash mismatch')
        require(component['component_id'] not in ids, 'duplicate component identity')
        ids.add(component['component_id'])
        n, r = int(component['summary']['degree']), component['blocks']
        require(n % r == 0, 'degree/block mismatch')
        k = n//r
        seed = analyze(parse_tuple(component['summary']['seed'], n), r, k)
        with table.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        require(len(rows) == int(component['summary']['states']), 'state count mismatch')
        require([int(row['state']) for row in rows] == list(range(len(rows))), 'state ids not complete/ordered')
        edges, histograms = [], {name: collections.Counter() for name in ['small_entry','quotient_shifts','geometry']}
        for row in rows:
            result = analyze(parse_tuple(row['tuple'], n), r, k)
            require(result['geometry'] == seed['geometry'] and result['source_genus'] == seed['source_genus'], 'geometry/genus changed within component')
            for field in ['genus', 'defect', 'product_index', 'left_orbits', 'right_orbits']:
                key = 'source_genus' if field == 'genus' else field
                require(int(row[field]) == result[key], f'cached {field} mismatch')
            require(row['smooth'] == 'true', 'cached smoothness mismatch')
            for name in histograms:
                value = result[name]
                histograms[name][str(value)] += 1
            e = [int(x) for x in row['edges'].split(',')]
            require(len(e) == 6 and all(0 <= x < len(rows) for x in e), 'invalid six edges')
            edges.append(e)
        for state, e in enumerate(edges):
            require(all(edges[target][(m+3) % 6] == state for m, target in enumerate(e)), 'edges not inverse paired')
        reached, queue = {0}, [0]
        for state in queue:
            for target in edges[state]:
                if target not in reached:
                    reached.add(target)
                    queue.append(target)
        require(len(reached) == len(rows), 'state graph not connected')
        require(seed['source_genus'] == int(component['summary']['genus']), 'summary genus mismatch')
        output.append(dict(folder=component['folder'], component_id=component['component_id'],
                           states_checked=len(rows), table_sha256=digest, seed=seed,
                           all_states_passed=True, state_histograms={k:dict(v) for k,v in histograms.items()},
                           group_order=component['independent_group_order'],
                           group_order_scope='reuses hash-bound C138 independent SymPy audit'))
        totals['components'] += 1
        totals['states'] += len(rows)
        totals[seed['geometry']+'_components'] += 1
        totals[seed['geometry']+'_states'] += len(rows)
        print(f"PASS {component['folder']}: {len(rows)} states, gC={seed['source_genus']}, gE={seed['quotient_genus']}, {seed['geometry']}", flush=True)
    require(totals['states'] == 22380 and totals['nonisotrivial_higher_genus_components'] == 10 and totals['nonisotrivial_higher_genus_states'] == 21648, 'unexpected declared totals')
    report = dict(status='all_geometric_hypotheses_passed', raw_root=str(raw),
                  quotient=catalog['quotient'], counts=dict(totals), components=output,
                  snapshot_sha256=plan['snapshot_sha256'], prior_checks_sha256=sha(root/'input/PRIOR_CHECKS.json'),
                  closure_scope='hash-bound C138 literal six-move closure reused; graph connectivity and inverse edges independently rechecked here',
                  classification_basis='PROOF.md: fixed quotient and exactly one moving branch value; not numerical moduli computation',
                  human_mathematical_review='pending; no approval inferred', formal_proof_assistant=False,
                  elapsed_seconds=time.monotonic()-start)
    write_json(root/'CHECKS.json', report)
    write_json(root/'FINISHED.json', dict(status='complete', finished=dt.datetime.now(dt.timezone.utc).isoformat(), checks_sha256=sha(root/'CHECKS.json')))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--output-parent', type=Path, default=Path('data'))
    parser.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        audit(args.worker)
        return
    catalog = json.loads(CATALOG.read_text())
    if not args.run:
        print(json.dumps(dict(mode='plan_only', components=catalog['component_count'],
                              states=sum(int(c['summary']['states']) for c in catalog['components']),
                              question='audit all52 geometric classifications, particularly all10 non-isotrivial claims',
                              cpu_workers=1, gpu=False), indent=2))
        return
    require(Path('data').is_dir() and args.output_parent.is_dir(), 'T7/output parent not mounted')
    lock = (PROJECT/'.run-control/.campaign.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = args.output_parent/('smooth_geometry_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    root.mkdir()
    snapshot = {}
    files = [(Path(__file__), 'source/scripts/verify_smooth_geometry.py'),
             (PROJECT/'scripts/smooth/test_smooth_geometry.py', 'source/scripts/test_smooth_geometry.py'),
             (CATALOG, 'input/CATALOG.json'), (PRIOR_CHECKS, 'input/PRIOR_CHECKS.json')]
    files += [(PROJECT/'evidence/smooth_geometry_20261012'/name, name) for name in ['PLAN.md','PROOF.md','SOURCES.json']]
    for source, destination in files:
        path = root/destination
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, path)
        snapshot[destination] = sha(path)
    raw = Path(catalog['raw_root'])
    write_json(root/'PLAN.json', dict(created=dt.datetime.now(dt.timezone.utc).isoformat(),
              cpu_workers=1, gpu=False, snapshot_sha256=snapshot,
              state_table_sha256={c['folder']:sha(raw/c['folder']/'STATES.tsv') for c in catalog['components']}))
    print(f'Audit root: {root}', flush=True)
    with (root/'RUN.log').open('w') as log:
        try:
            subprocess.run([sys.executable, str(root/'source/scripts/test_smooth_geometry.py')], stdout=log, stderr=subprocess.STDOUT, check=True)
            subprocess.run([sys.executable, str(root/'source/scripts/verify_smooth_geometry.py'), '--worker', str(root)], stdout=log, stderr=subprocess.STDOUT, check=True)
        except Exception as error:
            write_json(root/'FAILED.json', dict(status='failed', error=str(error)))
            raise
    print((root/'FINISHED.json').read_text(), end='')
    print(json.dumps(json.loads((root/'CHECKS.json').read_text())['counts'], indent=2))


if __name__ == '__main__':
    main()
