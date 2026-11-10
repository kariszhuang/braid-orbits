#!/usr/bin/env python3
"""Small checks for mission preservation, hash pinning and spectral populations."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from pathlib import Path
import tempfile

from analyze_defect_mission import spectral_summary
from run_defect_mission import completed_hashes, preserve_incomplete, sha, verify_hashes


def check():
    with tempfile.TemporaryDirectory(prefix='hurwitz-defect-check-') as temporary:
        root = Path(temporary)
        output = root / 'output'
        complete = output / 'case_000000'
        partial = output / 'case_000001'
        complete.mkdir(parents=True); partial.mkdir()
        (complete / 'DONE').write_text('closed')
        (complete / 'SUMMARY.tsv').write_text('unchanged')
        (partial / 'SUMMARY.tsv').write_text('partial retained')
        hashes = completed_hashes(output)
        preserve_incomplete(output)
        assert complete.is_dir() and not partial.exists()
        assert next((root / 'interrupted').glob('*/case_000001/SUMMARY.tsv')).read_text() == 'partial retained'
        verify_hashes(hashes)
        (complete / 'SUMMARY.tsv').write_text('modified')
        try:
            verify_hashes(hashes)
        except RuntimeError:
            pass
        else:
            raise AssertionError('changed checkpoint accepted')
        table = root / 'SPECTRUM.tsv'
        names = ['case_id','group','degree','kind','quotient_states','lambda2','lambda_min','genus','method','passport']
        with table.open('w', newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=names,delimiter='\t');writer.writeheader()
            base = dict(group='S',degree='8',kind='exact_component',quotient_states='12',genus='1',method='jacobi-exact',passport='a/b/c/d')
            writer.writerow(dict(base,case_id='1',lambda2='3',lambda_min='-6'))
            writer.writerow(dict(base,case_id='2',lambda2='0',lambda_min='-3'))
        summaries,n = spectral_summary([table])
        assert n == 2
        whole = next(row for row in summaries if row['size_bin']=='all')
        assert whole['relaxation_mean_component_weighted'] == .75
        assert whole['absolute_mean_component_weighted'] == .25
        try:
            spectral_summary([table, table])
        except ValueError:
            pass
        else:
            raise AssertionError('duplicated census accepted')
    print('defect mission checks passed')


if __name__ == '__main__':
    check()
