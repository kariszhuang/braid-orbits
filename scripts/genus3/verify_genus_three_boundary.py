#!/usr/bin/env python3
"""Independent coverage of a degree-16 two-path construction interval.

This verifies the finite constructor's coverage, not the global necessity of
the two chosen frames. Literal negative/positive certificates are audited by
the separate native auditor. Caps and timeouts remain unresolved.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import itertools
import json
from pathlib import Path

from run_genus_three import sha
from verify_genus_three_finite import mul, profile
from verify_smooth_geometry import parse_tuple

TOTAL = 504000


def candidates(start, stop):
    if not 0 <= start <= stop <= TOTAL:
        raise ValueError('boundary interval outside0..504000')
    offset = 0
    for split in (8, 7):
        pair = [list(range(16)), list(range(16))]
        for path, (base, length) in enumerate(((0, split), (split, 16-split))):
            for j in range(length-1):
                side = (j+int(split == 8 and path == 1)) % 2
                pair[side][base+j], pair[side][base+j+1] = base+j+1, base+j
        p = mul(pair[0], pair[1])
        omissions = itertools.product(range(8), range(8)) if split == 8 else itertools.combinations(range(7, 16), 2)
        for omitted in omissions:
            if offset+5040 <= start:
                offset += 5040
                continue
            if offset >= stop:
                return
            if split == 8:
                left = [x for x in range(8) if x != omitted[0]]
                right = [x for x in range(8, 16) if x != 8+omitted[1]]
            else:
                left = list(range(7))
                right = [x for x in range(7, 16) if x not in omitted]
            lo, hi = max(0, start-offset), min(5040, stop-offset)
            for rank, images in enumerate(itertools.islice(itertools.permutations(right), lo, hi), lo):
                c = list(range(16))
                for x, y in zip(left, images):
                    c[x], c[y] = y, x
                pc = mul(p, c)
                if len(profile(pc)) == 1:
                    fourth = [0]*16
                    for x, y in enumerate(pc):
                        fourth[y] = x
                    yield offset+rank, tuple(tuple(q) for q in (*pair, c, fourth))
            offset += 5040


def audit_boundary(folder, expected_start=None):
    folder = Path(folder)
    with (folder/'BOUNDARY.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if len(rows) != 1:
        raise ValueError('one boundary interval required')
    interval = {k: int(v) for k, v in rows[0].items()}
    start, stop = interval['start'], interval['next']
    if interval['total'] != TOTAL or expected_start is not None and start != expected_start:
        raise ValueError('boundary cursor binding')
    with (folder/'COUNTS.tsv').open() as stream:
        counts = next(csv.DictReader(stream, delimiter='\t'))
    if counts['lane'] != 'involution-boundary' or int(counts['degree']) != 16 or int(counts['constructed']) != stop-start:
        raise ValueError('boundary construction mass')
    qualified = 0
    with (folder/'LEADS.tsv').open() as stream:
        for cursor, t in candidates(start, stop):
            row = stream.readline().rstrip('\n').split('\t', 1)
            if len(row) != 2 or int(row[0]) != cursor-start+1:
                raise ValueError('missing/duplicate boundary candidate')
            if tuple(tuple(p) for p in parse_tuple(row[1], 16)) != t:
                raise ValueError('boundary tuple differs from independent constructor')
            qualified += 1
        if stream.readline():
            raise ValueError('extra boundary candidate')
    if qualified != interval['qualified'] or qualified != int(counts['scope']):
        raise ValueError('boundary qualified mass')
    if interval['unresolved'] != int(counts['unresolved']):
        raise ValueError('boundary unresolved mass')
    return dict(passed=True, **interval, input_sha256={name: sha(folder/name)
                for name in ('BOUNDARY.tsv', 'COUNTS.tsv', 'LEADS.tsv', 'RESOLUTIONS.tsv')},
                construction_complete=stop == TOTAL,
                qualification='selected two-path matching interval; caps retained, not whole degree16 census')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    print(json.dumps(audit_boundary(parser.parse_args().folder), indent=2))
