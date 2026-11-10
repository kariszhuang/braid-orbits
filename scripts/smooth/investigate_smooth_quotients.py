"""Audit quotient genera and test two small assumption-breaking lift strata.

The reported genus belongs to the explicitly selected block system. This is
not an enumeration of all block systems or a whole-degree existence census.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter, deque
import csv
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import tarfile
import time

from compact_smooth_outputs import sha
from verify_smooth_geometry import parse_tuple, compose, inverse, cycles, orbits, hurwitz, forced_partition
from run_s18_linked_roots import state_key


def index(p):
    return len(p) - len(cycles(p))


def index_regime(indices, n):
    a, b, c, d = sorted(indices)
    if a + c >= n or (a + d >= n and b + c >= n):
        raise ValueError('passport cannot be smooth throughout a braid component')
    return 'star-budget' if a + d < n else 'triangle-budget'


def quotient(t, blocks):
    n = len(t[0]); r = len(blocks)
    if len({len(b) for b in blocks}) != 1 or sorted(x for b in blocks for x in b) != list(range(n)):
        raise ValueError('invalid equal-size partition')
    k = len(blocks[0]); membership = {x: i for i, b in enumerate(blocks) for x in b}
    q = []
    for g in t:
        images = [{membership[g[x]] for x in b} for b in blocks]
        if any(len(s) != 1 for s in images): raise ValueError('partition is not invariant')
        q.append(tuple(next(iter(s)) for s in images))
    if len(orbits(t)) != 1 or len(orbits(q)) != 1: raise ValueError('not transitive')
    if compose(compose(compose(t[0], t[1]), t[2]), t[3]) != tuple(range(n)):
        raise ValueError('not product one')
    indices = [index(g) for g in t]; qindices = [index(g) for g in q]
    if sum(indices) % 2 or sum(qindices) % 2: raise ValueError('nonintegral genus')
    g = 1 - n + sum(indices) // 2; h = 1 - r + sum(qindices) // 2
    excess = [a - k * b for a, b in zip(indices, qindices)]
    if min(excess) < 0 or g != 1 + k * (h - 1) + sum(excess) // 2:
        raise ValueError('relative Riemann-Hurwitz identity')
    return dict(degree=n, blocks=r, block_size=k, genus=g, quotient_genus=h,
                point_indices=indices, quotient_indices=qindices,
                ramification_excess=excess, quotient_cycle_types=[sorted(map(len, cycles(p)), reverse=True) for p in q])


def profile_reports(root, archive=None):
    """Read bound original reports, including losslessly packed reports."""
    for path in sorted(root.glob('profile_*/PROFILES.json')):
        yield json.loads(path.read_text()), dict(path=str(path), sha256=sha(path))
    if archive is None:
        return
    recovery = json.loads((archive.parent / 'RECOVERY.json').read_text())
    if Path(recovery['root']).resolve() != root.resolve() or sha(archive) != recovery['archive_sha256']:
        raise ValueError('profile archive binding')
    records = {r['path']: r for r in recovery['members']}
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            parts = Path(member.name).parts
            if len(parts) != 2 or not parts[0].startswith('profile_') or parts[1] != 'PROFILES.json':
                continue
            if not member.isfile() or member.name not in records:
                raise ValueError('unexpected archived profile')
            raw = tar.extractfile(member).read()
            digest = hashlib.sha256(raw).hexdigest()
            if len(raw) != records[member.name]['size'] or digest != records[member.name]['sha256']:
                raise ValueError('archived profile checksum')
            yield json.loads(raw), dict(path=str(root / member.name), sha256=digest,
                                       archive=str(archive), archive_sha256=recovery['archive_sha256'])


def profile_summary(root, archive=None):
    queue_path = root / 'PROFILE_INPUT.json'; queue = json.loads(queue_path.read_text())
    binding = sha(queue_path); counts = Counter(); mechanisms = Counter(); examples = {}; sources = []
    seen = set(); towers = []; regimes = Counter()
    for report, source in profile_reports(root, archive):
        if report['queue_sha256'] != binding: raise ValueError('profile queue binding')
        for e in report['components']:
            i = e['queue_index']
            if i in seen or e['table_sha256'] != queue['components'][i]['table_sha256']:
                raise ValueError('duplicate or rebound profile')
            seen.add(i)
            representative = e['representative']
            if e['mechanism'] == 'ramified-double-lift':
                # Historical profiles expressed the partition in state0 labels,
                # but displayed the constructor seed. Read its bound state0.
                table = Path(e['raw_folder']) / 'STATES.tsv'
                if sha(table) != e['table_sha256']: raise ValueError('ramified state0 table changed')
                with table.open() as stream:
                    representative = next(csv.DictReader(stream, delimiter='\t'))['tuple']
            t = parse_tuple(representative, e['degree'])
            blocks = [[x - 1 for x in b] for b in e['block_partition']]
            value = quotient(t, blocks)
            if value['quotient_genus'] != e.get('quotient_genus', e.get('base_genus')):
                raise ValueError('independent quotient genus disagreement')
            counts[value['quotient_genus']] += 1
            regimes[index_regime(value['point_indices'], e['degree'])] += 1
            mechanisms[(e['mechanism'], value['quotient_genus'])] += 1
            examples.setdefault((e['mechanism'], value['quotient_genus']),
                                dict(value, representative=representative, block_partition=e['block_partition']))
            if e['mechanism'] == 'ramified-double-lift':
                membership = {x: j for j, b in enumerate(blocks) for x in b}
                q = tuple(tuple(membership[g[b[0]]] for b in blocks) for g in t)
                support = min((sorted(x for x in range(len(q[0])) if g[x] != x)
                               for g in q if g != tuple(range(len(q[0])))), key=len)
                coarse = forced_partition(q, [support])
                if len(coarse) != 3 or any(len(b) != 3 for b in coarse):
                    raise ValueError('degree18 example quotient tower changed')
                merged = [sorted(x for j in b for x in blocks[j]) for b in coarse]
                coarse_value = quotient(t, merged)
                towers.append(dict(queue_index=i, fine=value, coarse=coarse_value,
                                   coarse_block_partition=[[x+1 for x in b] for b in merged]))
        sources.append(source)
    if seen != set(range(len(queue['components']))): raise ValueError('incomplete profile population')
    return dict(passed=True, components=len(seen), by_quotient_genus=dict(counts),
                by_mechanism=[dict(mechanism=a, quotient_genus=b, components=c) for (a, b), c in sorted(mechanisms.items())],
                examples=list(examples.values()), input_queue_sha256=binding, source_reports=sources,
                example_quotient_towers=towers,
                index_budget_regimes=dict(regimes),
                scope='one certified chosen block system per queue component; not all possible factorization genera')


def smooth(t):
    p = index(compose(t[0], t[1])); inds = [index(g) for g in t]
    return p in [inds[0] + inds[1], inds[2] + inds[3]]


def text_tuple(t):
    return '[' + ', '.join(''.join('(' + ' '.join(str(x + 1) for x in c) + ')'
                                   for c in cycles(p) if len(c) > 1) or '()' for p in t) + ']'


def close_or_reject(t, cap, deadline):
    tuples = [t]; keys = {state_key(t): 0}; paths = [[]]; adjacency = []
    for i, current in enumerate(tuples):
        edges = []
        for m in range(6):
            moved = hurwitz(current, m)
            if not smooth(moved): return dict(status='nonsmooth', word=paths[i] + [m])
            key = state_key(moved)
            if key not in keys:
                if len(tuples) >= cap or time.monotonic() > deadline:
                    return dict(status='unresolved', states_seen=len(tuples))
                keys[key] = len(tuples); tuples.append(moved); paths.append(paths[i] + [m])
            edges.append(keys[key])
        adjacency.append(edges)
    return dict(status='closed_smooth', states=len(tuples), tuples=[text_tuple(t) for t in tuples], edges=adjacency)


def from_cycles(n, parts):
    p = list(range(n))
    for c in parts:
        for x, y in zip(c, c[1:] + c[:1]): p[x] = y
    return tuple(p)


def lift(voltages):
    r = len(voltages); k = len(voltages[0])
    return tuple(((j + 1) % r) * k + voltages[j][x] for j in range(r) for x in range(k))


def two_active_blocks():
    a = from_cycles(12, [[0, 1], [5, 6]]); b = lift([tuple(range(4))] * 3)
    for x, y in itertools.product(itertools.permutations(range(4)), repeat=2):
        z = inverse(compose(y, x)); c = lift([x, y, z])
        d = inverse(compose(compose(a, b), c))
        yield (a, b, c, d), [3] * 4


def genus_zero_quotient():
    a = from_cycles(12, [[0, 1, 2], [6, 7, 8]])
    b = from_cycles(12, [[2, 3, 4, 5], [8, 9, 10, 11]])
    for p in itertools.permutations(range(6)):
        for i, j in itertools.combinations(range(6), 2):
            u = from_cycles(6, [[i, j]]); q = compose(u, inverse(p)); c = lift([p, q])
            d = inverse(compose(compose(a, b), c))
            yield (a, b, c, d), [2, 2, 2, 2, 4]


def finite_pilot(generator, total, seconds=20):
    deadline = time.monotonic() + seconds; start = time.monotonic(); counts = Counter(); results = []; roots = {}
    for t, dtype in generator():
        counts['constructed'] += 1
        if sorted(map(len, cycles(t[3]))) != dtype:
            counts['fourth_type_reject'] += 1; continue
        if len(orbits(t)) != 1:
            counts['intransitive'] += 1; continue
        if not smooth(t): raise ValueError('constructed first pair is not smooth')
        counts['qualifying'] += 1; key = state_key(t)
        if key in roots:
            counts['seed_alias'] += 1; continue
        if time.monotonic() > deadline:
            result = dict(status='unresolved', reason='global pilot limit')
        else:
            result = close_or_reject(t, 5000, deadline)
        roots[key] = result; counts[result['status']] += 1
        results.append(dict(seed=text_tuple(t), **result))
    if counts['constructed'] != total: raise ValueError('finite construction coverage')
    return dict(counts=dict(counts), results=results, seconds=time.monotonic() - start,
                scope='fixed first pair and declared block-voltage population; not a whole-degree census')


def audit_pilot_negatives(result, out, auditor):
    """Replay literal negatives in the independent Rust kernel, not this code."""
    negatives = [e for e in result['results'] if e['status'] == 'nonsmooth']
    out.mkdir()
    leads = []; resolutions = ['candidate\tstatus\tword\tfolder\tseed\n']
    for i, e in enumerate(negatives, 1):
        leads.append(f"{i}\t{e['seed']}\n")
        resolutions.append(f"{i}\tnonsmooth\t{','.join(map(str, e['word']))}\t\t{e['seed']}\n")
    (out / 'LEADS.tsv').write_text(''.join(leads))
    (out / 'RESOLUTIONS.tsv').write_text(''.join(resolutions))
    count = len(negatives)
    (out / 'COUNTS.tsv').write_text(
        'degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'
        f'12\treplay\t0\t{count}\t{count}\t{count}\t{count}\t{count}\t0\t0\t0\t0\n')
    report = json.loads(subprocess.check_output([str(auditor), '--guided', str(out), '12'], text=True))
    if report['records'] != count or report['nonsmooth'] != count or report['unresolved'] or report['closed']:
        raise ValueError('independent literal negative audit')
    report.update(binary=str(auditor), binary_sha256=sha(auditor),
                  source_sha256=sha(out / 'RESOLUTIONS.tsv'),
                  scope='all saved literal negatives only; construction coverage separately declared')
    (out / 'CHECKS.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True); p.add_argument('--pilots', action='store_true')
    p.add_argument('--archive', type=Path, help='verified lossless archive containing finished profile reports')
    p.add_argument('--auditor', type=Path, help='independent smooth_negative_audit binary for pilot witnesses')
    args = p.parse_args()
    args.out.mkdir()
    summary = profile_summary(args.root.resolve(), args.archive.resolve() if args.archive else None)
    (args.out / 'QUOTIENT_GENERA.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('QUOTIENTS', summary['by_quotient_genus'], flush=True)
    if args.pilots:
        first = finite_pilot(two_active_blocks, 576)
        if args.auditor:
            first['independent_negative_audit'] = audit_pilot_negatives(first, args.out / 'two_active_negative_audit', args.auditor.resolve())
        (args.out / 'TWO_ACTIVE_BLOCKS.json').write_text(json.dumps(first, indent=2) + '\n')
        print('TWO_ACTIVE_BLOCKS', first['counts'], flush=True)
        second = finite_pilot(genus_zero_quotient, 10800)
        if args.auditor:
            second['independent_negative_audit'] = audit_pilot_negatives(second, args.out / 'genus_zero_negative_audit', args.auditor.resolve())
        (args.out / 'GENUS_ZERO_QUOTIENT.json').write_text(json.dumps(second, indent=2) + '\n')
        print('GENUS_ZERO_QUOTIENT', second['counts'], flush=True)


if __name__ == '__main__': main()
