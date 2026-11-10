#!/usr/bin/env python3
"""Dependency-free descriptive reconciliation; no pretend independent walk samples."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import statistics

PROJECT = Path(__file__).resolve().parents[2]


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def spectral_summary(tables):
    blocks = defaultdict(list)
    seen = set()
    for table in tables:
        for row in rows(table):
            if row['case_id'] in seen:
                raise ValueError('Duplicate spectral case: ' + row['case_id'])
            seen.add(row['case_id'])
            n = int(row['quotient_states'])
            if row['group'] != 'S' or row['degree'] != '8' or row['kind'] != 'exact_component':
                raise ValueError('Wrong S8 spectral population')
            l2, minimum = float(row['lambda2']), float(row['lambda_min'])
            relaxation = 1-l2/6
            absolute = 1-max(abs(l2), abs(minimum))/6
            if not all(math.isfinite(x) for x in (relaxation, absolute)) or min(relaxation, absolute) < -1e-8:
                raise ValueError('Invalid eigenvalue/gap row')
            # These are floating-point values, not new interval certificates.
            row = dict(row, relaxation=max(0, relaxation), absolute=max(0, absolute))
            bucket = '1-24' if n <= 24 else '25-128' if n <= 128 else '129-512' if n <= 512 else '513+'
            for size in ('all', bucket):
                blocks[(row['genus'], row['method'], size)].append(row)
    result = []
    for (genus, method, size), block in sorted(blocks.items()):
        passports = defaultdict(list)
        for row in block:
            passport = '/'.join(sorted(row['passport'].split('/')))
            passports[passport].append(row)
        item = dict(genus=genus, method=method, size_bin=size, components=len(block), passports=len(passports))
        for key in ('relaxation', 'absolute'):
            item[key+'_mean_component_weighted'] = statistics.mean(r[key] for r in block)
            item[key+'_median_component_weighted'] = statistics.median(r[key] for r in block)
            item[key+'_mean_passport_balanced'] = statistics.mean(statistics.mean(r[key] for r in p) for p in passports.values())
        result.append(item)
    return result, len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plan = json.loads((args.root / 'PLAN.json').read_text())
    tables = [Path(p) for p in plan['input_hashes'] if 'S8_exact_shards/' in p and p.endswith('/SPECTRUM.tsv')]
    if len(tables) != 16:
        raise ValueError('Need all sixteen completed S8 spectral tables')
    summary, count = spectral_summary(tables)
    if count != 8676:
        raise ValueError('S8 census size changed')
    with (args.out / 'S8_GENUS_SPECTRA.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(summary)
    statuses = Counter()
    for marker in args.root.glob('jobs/*/output/draw_*/DONE'):
        if marker.read_text().startswith('no_eligible_seed'):
            statuses['no_eligible_seed'] += 1
    components, seen = [], set()
    for table in sorted(args.root.glob('jobs/*/output/*/SUMMARY.tsv')):
        if not (table.parent / 'DONE').is_file():
            continue
        row = rows(table)[0]
        statuses[row['status']] += 1
        # Closed components deduplicate by the minimum exact state key, not seed.
        if row['status'] == 'closed':
            key = (row['group'], row['degree'], row['representative'])
            if key in seen:
                continue
            seen.add(key)
            hist = rows(table.parent / 'DEFECT.tsv')
            n = int(row['S_quotient_states'])
            if sum(int(r['count']) for r in hist) != n:
                raise ValueError('Closed histogram count mismatch')
            row.update(defect_min=min(int(r['defect']) for r in hist), defect_max=max(int(r['defect']) for r in hist),
                       defect_mean=sum(int(r['defect'])*int(r['count']) for r in hist)/n,
                       zero_fraction=sum(int(r['count']) for r in hist if r['defect']=='0')/n)
            components.append(row)
    if components:
        with (args.out / 'EXACT_COMPONENTS.tsv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(components[0]), delimiter='\t')
            writer.writeheader(); writer.writerows(components)
    (args.out / 'STATUS_COUNTS.json').write_text(json.dumps(dict(statuses), indent=2)+'\n')
    (args.out / 'METHODS.md').write_text(
        '# Scope\n\nAll 8,676 stored S8 spectral rows, split by genus, graph size, and numeric method.\n'
        'Component- and passport-balanced means are descriptive finite-catalog summaries, not causal genus effects.\n'
        'Exact defect components are deduplicated by minimum full canonical-state key. Inner classes may split;\n'
        'S-quotient state percentages lift unchanged to the reachable G-quotient component.\n'
        'Partial BFS histograms and early-stopped lazy walks are not uniform component distributions.\n'
        'Status counts include repeated inputs; they are not distinct orbit counts. No-eligible-seed draws\n'
        'are recorded separately in each draw/DONE and imply no mathematical exclusion.\n'
        'The null stage identifies exact move-labeled classes and applies existing class-balanced analysis.\n'
        'Null selection is stratified and limited to graphs of at most 512 vertices; it is not all S8 classes.\n'
        'Lanczos rows are estimates; Jacobi rows are full floating-point diagonalizations, not exact algebra.\n')


if __name__ == '__main__':
    main()
