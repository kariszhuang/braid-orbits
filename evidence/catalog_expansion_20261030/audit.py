#!/usr/bin/env python3
"""Independent complete-certificate audit and exact Hodge classification.

The native auditor imports no producer modules and uses a different canonical
ordering. This Python layer reads the literal state table independently.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def atomic(p, value):
    p = Path(p)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(p)


def parse(s, n):
    pieces = s.strip()[1:-1].split(',')
    assert len(pieces) == 4
    out = []
    for part in pieces:
        p = list(range(n))
        used = set()
        for text in re.findall(r'\(([^)]*)\)', part):
            cyc = [int(x) - 1 for x in text.split()]
            assert len(set(cyc)) == len(cyc) and all(0 <= x < n for x in cyc)
            assert not used.intersection(cyc)
            used.update(cyc)
            for x, y in zip(cyc, cyc[1:] + cyc[:1]):
                p[x] = y
        out.append(tuple(p))
    return tuple(out)


def mul(a, b):
    return tuple(a[x] for x in b)


def lengths(p):
    seen = set()
    out = []
    for x in range(len(p)):
        if x in seen:
            continue
        k = 0
        while x not in seen:
            seen.add(x)
            k += 1
            x = p[x]
        out.append(k)
    return tuple(sorted(out, reverse=True))


def notation(p):
    seen = set()
    out = []
    for x in range(len(p)):
        if x in seen:
            continue
        c = []
        while x not in seen:
            seen.add(x)
            c.append(x + 1)
            x = p[x]
        if len(c) > 1:
            out.append('(' + ' '.join(map(str, c)) + ')')
    return ''.join(out) or '()'


def block_witness(t):
    """A proper minimal invariant equivalence proves imprimitive.

    Trying every pair (0,x) proves primitive if no proper one exists. This
    provides one witness, not an enumeration of the full block lattice.
    """
    n = len(t[0])
    for x in range(1, n):
        par = list(range(n))
        def root(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        def join(a, b):
            a, b = root(a), root(b)
            if a == b:
                return False
            par[b] = a
            return True
        join(0, x)
        changed = True
        while changed:
            changed = False
            cl = {}
            for u in range(n):
                cl.setdefault(root(u), []).append(u)
            for p in t:
                for row in cl.values():
                    for u in row[1:]:
                        changed = join(p[row[0]], p[u]) or changed
        cl = {}
        for u in range(n):
            cl.setdefault(root(u), []).append(u)
        blocks = sorted(cl.values())
        if 1 < len(blocks) < n:
            assert len({len(b) for b in blocks}) == 1
            where = {u: i for i, b in enumerate(blocks) for u in b}
            q = tuple(tuple(where[p[b[0]]] for b in blocks) for p in t)
            assert all(all(where[p[u]] == q[j][i] for u in b)
                       for j, p in enumerate(t) for i, b in enumerate(blocks))
            iq = [len(q[0]) - len(lengths(p)) for p in q]
            return False, dict(block_count=len(blocks), block_size=len(blocks[0]),
                blocks=[[u + 1 for u in b] for b in blocks], tuple=[notation(p) for p in q],
                indices=iq, genus=1-len(blocks)+sum(iq)//2,
                scope='one invariant block-system witness; full lattice not yet profiled')
    return True, None


def audit(folder, auditor, compress=True):
    folder = Path(folder)
    with (folder / 'SUMMARY.tsv').open() as f:
        summary = next(csv.DictReader(f, delimiter='\t'))
    assert summary['status'] == 'closed_smooth'
    n, count = int(summary['degree']), int(summary['states'])
    records = folder / 'INDEPENDENT_KEYS.tsv'
    if records.exists():
        raise ValueError('audit output already exists; inspect previous verification first')
    p = subprocess.run([str(auditor), '--component', str(folder), str(n), str(records)],
                       capture_output=True, text=True, check=True)
    native = json.loads(p.stdout)
    assert native['closed'] and native['all_six_edges_checked'] and native['states'] == count
    min_key = None
    for line in records.read_text().splitlines():
        key = json.loads(line.split('\t', 1)[0])
        if min_key is None or key < min_key:
            min_key = key
    component_id = hashlib.sha256(json.dumps([n, min_key], separators=(',', ':')).encode()).hexdigest()
    t = parse(summary['seed'], n)
    types = [lengths(p) for p in t]
    ix = [n - len(c) for c in types]
    genus = 1 - n + sum(ix) // 2
    assert genus == int(summary['genus']) and genus >= 2 and all(i > 0 for i in ix)
    hist = Counter()
    with (folder / 'STATES.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            u = parse(row['tuple'], n)
            assert row['smooth'] == 'true' and row['defect'] == '0'
            assert int(row['genus']) == genus
            hist[lengths(mul(u[0], u[1]))] += 1
    assert sum(hist.values()) == count
    rho = lambda c: sum((Fraction(1, k) for k in c), Fraction(0))
    # Full orbit (all six moves) makes the three cusp-sheet means equal.
    l12 = n - sum((rho(c) for c in types), Fraction(0))
    l12 += 3 * sum((num * rho(c) for c, num in hist.items()), Fraction(0)) / count
    assert l12 >= 0, 'negative Hodge degree: quarantine, do not classify'
    primitive, witness = block_witness(t)
    input_record = json.loads((folder / 'INPUT.json').read_text())
    result = dict(id=component_id, degree=n, genus=genus, states=count,
        tuple=[notation(p) for p in t], indices=ix, cycle_types=[list(c) for c in types],
        primitive=primitive, group_order=None, group_order_status='deferred',
        contained_in=('A' if all(i % 2 == 0 for i in ix) else 'S') + str(n),
        all_states_smooth=True, isotrivial=l12 == 0,
        isotriviality=dict(status='isotrivial' if l12 == 0 else 'non-isotrivial',
            method='exact full-orbit normalized Hodge degree', twelve_L=str(l12), L=str(l12 / 12)),
        product_cycle_distribution=[dict(cycle_lengths=list(c), states=num) for c, num in sorted(hist.items())],
        block_witness=witness, source=input_record['source'],
        last_updated=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        certificate=dict(folder=str(folder), states_sha256=sha(folder / 'STATES.tsv'),
            independent_keys_sha256=sha(records), independent_native=native,
            auditor_sha256=sha(auditor), audit_script_sha256=sha(Path(__file__))))
    if compress:
        for src in [folder / 'STATES.tsv', records]:
            dst = src.with_suffix(src.suffix + '.gz')
            with src.open('rb') as a, gzip.open(dst, 'wb', compresslevel=3) as b:
                shutil.copyfileobj(a, b)
            # Read back the entire compressed stream before dropping the plain copy.
            h = hashlib.sha256()
            with gzip.open(dst, 'rb') as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b''):
                    h.update(chunk)
            assert h.hexdigest() == sha(src)
            src.unlink()
        result['certificate']['states_file'] = 'STATES.tsv.gz'
        result['certificate']['independent_keys_file'] = 'INDEPENDENT_KEYS.tsv.gz'
    atomic(folder / 'VERIFIED.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder', required=True, type=Path)
    p.add_argument('--auditor', required=True, type=Path)
    p.add_argument('--no-compress', action='store_true')
    a = p.parse_args()
    r = audit(a.folder, a.auditor, not a.no_compress)
    print(json.dumps({k: r[k] for k in ['id', 'degree', 'genus', 'states', 'isotrivial', 'source']}))


if __name__ == '__main__':
    main()
