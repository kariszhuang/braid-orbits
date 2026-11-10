#!/usr/bin/env python3
"""Bounded, local-only audit of a frozen smooth catalog; no discovery/deployment.

Recomputes tuple and quotient arithmetic from literal permutations, records
fixed/moving branch witnesses, and tests a proposed stronger geometric lemma.
The finite checks do not verify the analytic proof or rerun braid closure.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import time

PROJECT = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def permutation(text, n):
    if re.sub(r'\([^()]*\)', '', text).strip():
        raise ValueError('unexpected permutation syntax')
    p, used = list(range(n)), set()
    for body in re.findall(r'\(([^()]*)\)', text):
        c = [int(x)-1 for x in body.split()]
        if len(c) != len(set(c)) or any(x < 0 or x >= n or x in used for x in c):
            raise ValueError('invalid disjoint cycles')
        used.update(c)
        for x, y in zip(c, c[1:]+c[:1]):
            p[x] = y
    return tuple(p)


def compose(a, b):
    return tuple(a[x] for x in b)


def cycles(p):
    seen, result = set(), []
    for x in range(len(p)):
        if x in seen:
            continue
        c, y = [], x
        while y not in seen:
            seen.add(y)
            c.append(y)
            y = p[y]
        result.append(c)
    return result


def orbit_sizes(perms):
    unseen, result = set(range(len(perms[0]))), []
    while unseen:
        queue, seen = [min(unseen)], {min(unseen)}
        for x in queue:
            for p in perms:
                if p[x] not in seen:
                    seen.add(p[x])
                    queue.append(p[x])
        unseen.difference_update(seen)
        result.append(len(seen))
    return sorted(result, reverse=True)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def quotient_check(t, blocks):
    n, r = len(t[0]), len(blocks)
    m = len(blocks[0])
    require(r*m == n and all(len(b) == m for b in blocks), 'block sizes')
    require(sorted(x for b in blocks for x in b) == list(range(1, n+1)), 'partition')
    labels = {x-1: i for i, b in enumerate(blocks) for x in b}
    q = []
    for p in t:
        images = [{labels[p[x-1]] for x in b} for b in blocks]
        require(all(len(v) == 1 for v in images), 'invariant blocks')
        q.append(tuple(next(iter(v)) for v in images))
    qi = [r-len(cycles(p)) for p in q]
    indices = [n-len(cycles(p)) for p in t]
    excess = [a-m*b for a, b in zip(indices, qi)]
    require(sum(qi) % 2 == 0 and min(excess) >= 0, 'quotient parity/excess')
    h = 1-r+sum(qi)//2
    branch = []
    # A quotient cycle of length ell is one point Q of E. Upstairs lengths
    # ell*e give the complete local degree partition at Q, not merely excess.
    for p, qp in zip(t, q):
        pc, row = cycles(p), []
        for bc in cycles(qp):
            bs, ell = set(bc), len(bc)
            local = sorted([len(c)//ell for c in pc if labels[c[0]] in bs], reverse=True)
            require(sum(local) == m and all(len(c) % ell == 0 for c in pc if labels[c[0]] in bs), 'local degrees')
            e = sum(x-1 for x in local)
            if e:
                row.append(dict(quotient_cycle=[x+1 for x in bc], local_degrees=local, excess=e))
        branch.append(row)
    require([sum(a['excess'] for a in row) for row in branch] == excess, 'branch sum')
    return q, qi, excess, h, branch, labels


def inspect(entry):
    n = entry['degree']
    t = [permutation(s, n) for s in entry['tuple']]
    ident = tuple(range(n))
    require(all(p != ident for p in t), 'identity-free catalog')
    require(compose(compose(compose(*t[:2]), t[2]), t[3]) == ident, 'product one')
    require(orbit_sizes(t) == [n], 'transitivity')
    indices = [n-len(cycles(p)) for p in t]
    require(indices == entry['indices'] and sum(indices) % 2 == 0, 'indices')
    g = 1-n+sum(indices)//2
    require(g > 0 and g == entry['genus'], 'source genus')
    ab = compose(t[0], t[1])
    ip = n-len(cycles(ab))
    require(ip in (sum(indices[:2]), sum(indices[2:])), 'representative smoothness')
    pair_sizes = [orbit_sizes(t[i:i+2]) for i in range(3)]
    defect = len(cycles(ab))-len(pair_sizes[0])-len(pair_sizes[2])+1
    require(defect == 0, 'smooth implies defect zero control')
    require(entry['all_states_smooth'] is True, 'closure marker')
    systems, elliptic, hyperbolic, old = [], [], [], []
    for bi, b in enumerate(entry['block_quotients']):
        q, qi, excess, h, branch, labels = quotient_check(t, b['blocks'])
        m, r = b['block_size'], b['block_count']
        require(len(q[0]) == r and len(t[0]) == r*m, 'quotient shape')
        require(q == [permutation(s, r) for s in b['tuple']], 'quotient literal')
        require(qi == b['indices'] and excess == b['ramification_excess'] and h == b['genus'], 'quotient arithmetic')
        require(2*g-2 == m*(2*h-2)+sum(excess), 'relative Riemann-Hurwitz')
        zero = [k for k, p in enumerate(q) if p == tuple(range(r))]
        record = dict(index=bi, block_count=r, block_size=m, genus=h,
                      indices=qi, excess=excess, identity_entries=[k+1 for k in zero],
                      blocks=b['blocks'], quotient_tuple=b['tuple'])
        systems.append(record)
        for k in zero:
            active = sorted({labels[x] for x in range(n) if t[k][x] != x})
            fixed = sum(len(branch[j]) for j in range(4) if j != k)
            moving = len(branch[k])
            w = dict(**record, special_entry=k+1, active_blocks=[x+1 for x in active],
                     fixed_relative_branch_points=fixed, moving_relative_branch_points=moving,
                     relative_branch_local_degrees=branch)
            # All source entries are nonidentity, so an identity quotient entry
            # creates ramification over the moving base point.
            require(excess[k] > 0 and moving > 0, 'moving ramification')
            if len(active) == 1 and fixed == 0:
                require(h >= 1, 'original sufficient theorem genus bound')
                old.append(w)
            if h == 1 and fixed > 0:
                elliptic.append(w)
            if h >= 2:
                hyperbolic.append(w)
    require(bool(systems), 'known catalog is imprimitive')
    # A recorded quotient is coarsest when no other proper recorded partition
    # coarsens it. Completeness is inherited from the frozen export audit.
    coarsest = []
    for i, a in enumerate(entry['block_quotients']):
        if not any(b['block_count'] < a['block_count'] and
                   all(any(set(x) <= set(y) for y in b['blocks']) for x in a['blocks'])
                   for b in entry['block_quotients']):
            coarsest.append(systems[i])
    old_decisions = {w['genus'] == 1 for w in old}
    require(len(old_decisions) <= 1, 'conflicting original witnesses')
    old_value = next(iter(old_decisions)) if old_decisions else None
    require(old_value is entry['isotrivial'], 'original classification')
    require(not (old_value is True and (elliptic or hyperbolic)), 'new/old contradiction')
    a, b, c, d = sorted(indices)
    require(a+c <= n-1 and min(a+d, b+c) <= n-1, 'necessary passport budget')
    return dict(id=entry['id'], degree=n, genus=g, states=entry['states'],
                original_isotrivial=old_value, systems=systems, coarsest=coarsest,
                elliptic_fixed_and_moving_witnesses=elliptic,
                hyperbolic_moving_witnesses=hyperbolic,
                star_budget=a+d <= n-1, tuple=entry['tuple'], indices=indices,
                representative_pair_orbit_sizes=pair_sizes)


def audit(catalog):
    started = time.monotonic()
    data = json.loads(Path(catalog).read_text())
    entries = data['components']
    require(len({e['id'] for e in entries}) == len(entries), 'distinct IDs')
    rows = [inspect(e) for e in entries]
    unknown = [r for r in rows if r['original_isotrivial'] is None]
    elliptic_new = [r for r in unknown if r['elliptic_fixed_and_moving_witnesses']]
    hyperbolic_new = [r for r in unknown if r['hyperbolic_moving_witnesses']]
    counts = Counter('true' if r['original_isotrivial'] is True else
                     'false' if r['original_isotrivial'] is False else 'unknown' for r in rows)
    degree_rows = []
    for d in sorted({r['degree'] for r in rows}):
        rs = [r for r in rows if r['degree'] == d]
        degree_rows.append(dict(degree=d, components=len(rs), states=sum(r['states'] for r in rs),
            genus_range=[min(r['genus'] for r in rs), max(r['genus'] for r in rs)],
            true=sum(r['original_isotrivial'] is True for r in rs),
            false=sum(r['original_isotrivial'] is False for r in rs),
            unknown=sum(r['original_isotrivial'] is None for r in rs)))
    profile = Counter()
    for r in elliptic_new:
        w = min(r['elliptic_fixed_and_moving_witnesses'], key=lambda w: (w['block_count'], w['block_size'], w['special_entry']))
        k = w['special_entry']-1
        profile[(r['degree'], w['block_count'], w['block_size'], r['genus'],
                 tuple(sorted(e for j, e in enumerate(w['excess']) if j != k and e)),
                 w['moving_relative_branch_points'], w['fixed_relative_branch_points'])] += 1
    stats = dict(dataset_version=data['dataset_version'], snapshot_sha256=digest(catalog),
        code_sha256=digest(__file__), components=len(rows), states=sum(r['states'] for r in rows),
        original_classifications=dict(counts), by_degree=degree_rows,
        min_genus=min(r['genus'] for r in rows), max_genus=max(r['genus'] for r in rows),
        star_budget_components=sum(r['star_budget'] for r in rows),
        triangle_budget_components=sum(not r['star_budget'] for r in rows),
        coarsest_systems=sum(len(r['coarsest']) for r in rows),
        coarsest_without_identity=sum(not q['identity_entries'] for r in rows for q in r['coarsest']),
        proper_block_systems=sum(len(r['systems']) for r in rows),
        quotient_genus_incidence=dict(sorted(Counter(str(q['genus']) for r in rows for q in r['systems']).items())),
        components_with_rational_quotient=sum(any(q['genus'] == 0 for q in r['systems']) for r in rows),
        unknown_components=len(unknown), unknown_states=sum(r['states'] for r in unknown),
        unknown_with_elliptic_fixed_and_moving=len(elliptic_new),
        unknown_with_hyperbolic_moving=len(hyperbolic_new),
        unknown_with_single_active_elliptic_witness=sum(any(len(w['active_blocks']) == 1 for w in r['elliptic_fixed_and_moving_witnesses']) for r in unknown),
        unknown_without_either=sum(not r['elliptic_fixed_and_moving_witnesses'] and not r['hyperbolic_moving_witnesses'] for r in unknown),
        proposed_nonisotrivial_total=counts['false']+len({r['id'] for r in elliptic_new+hyperbolic_new}),
        original_nonisotrivial_min_genus=min(r['genus'] for r in rows if r['original_isotrivial'] is False),
        proposed_nonisotrivial_min_genus=min(r['genus'] for r in rows if r['original_isotrivial'] is False or r['elliptic_fixed_and_moving_witnesses'] or r['hyperbolic_moving_witnesses']),
        by_genus=[dict(genus=g, components=sum(r['genus'] == g for r in rows),
            true=sum(r['genus'] == g and r['original_isotrivial'] is True for r in rows),
            false=sum(r['genus'] == g and r['original_isotrivial'] is False for r in rows),
            unknown=sum(r['genus'] == g and r['original_isotrivial'] is None for r in rows))
            for g in sorted({r['genus'] for r in rows})],
        elliptic_unknown_profiles=[dict(degree=p[0], block_count=p[1], block_size=p[2], source_genus=p[3],
            ordinary_positive_excesses=p[4], moving_branch_points=p[5], fixed_branch_points=p[6], components=count)
            for p, count in sorted(profile.items())],
        elapsed_seconds=time.monotonic()-started,
        scope='Literal representatives and recorded complete block lattices. Existing full closures reused, not rerun. Geometric proof requires human review; no public labels changed.')
    return stats, unknown


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--catalog', type=Path, default=PROJECT/'exports/latest-discoveries-20261014/release-v5/all-smooth.json')
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    stats, unknown = audit(args.catalog)
    for name, value in [('SUMMARY.json', stats), ('UNKNOWN_WITNESSES.json', unknown)]:
        (args.out/name).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
