try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from audit_guided_native import audit_guided_native,known_keys
from run_smooth_degree9 import sha
from verified_smooth_bank import VerifiedSmoothBank
from test_guided_native import fixture,seed_text

PROJECT=Path(__file__).resolve().parents[2]
AUDITOR=Path(os.environ.get('GUIDED_AUDIT_BINARY',PROJECT/'target/theory-frontier/release/smooth_negative_audit'))


def positive_fixture(root,status='closed_smooth'):
    c=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
    fixture(root,status,'')
    for name in ['LEADS.tsv','RESOLUTIONS.tsv']:
        p=root/name;p.write_text(p.read_text().replace(seed_text(),c['representative']))
    shutil.copytree(c['raw_folder'],root/'component_1')
    if status!='closed_smooth':
        p=root/'component_1/SUMMARY.tsv';p.write_text(p.read_text().replace('closed_smooth',status))
    check=audit_guided_native(root,AUDITOR,collect_records=True)
    records=check.pop('_closed_records')
    (root/'CHECKS.json').write_text(json.dumps(check,indent=2)+'\n')
    return check,records


class TheoryFrontierTests(unittest.TestCase):
    def test_verified_publication_keeps_prior_snapshots_immutable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);job=root/'job';job.mkdir();check,records=positive_fixture(job)
            initial=root/'initial.tsv';initial.write_text('')
            bank=VerifiedSmoothBank(root/'bank',initial,sha(initial));old,digest=bank.snapshot(9)
            with self.assertRaisesRegex(ValueError,'unverified'):bank.publish(job,check,records)
            (job/'VERIFIED').write_text('checked\n')
            forged=records[:];forged[0]=(9,records[0][1],'false')
            with self.assertRaisesRegex(ValueError,'records binding'):bank.publish(job,check,forged)
            result=bank.publish(job,check,records);new,new_digest=bank.snapshot(9)
            self.assertEqual(result['added'],12);self.assertNotEqual(old,new)
            self.assertEqual(sha(old),digest);self.assertEqual(sha(new),new_digest)
            self.assertEqual(len(known_keys(new)),12)
            self.assertEqual(bank.publish(job,check,records)['added'],0)
            (job/'component_1/STATES.tsv').write_text('corrupted')
            with self.assertRaisesRegex(ValueError,'table changed'):bank.publish(job,check,records)

    def test_partial_not_published_and_ram_bound_never_discards_lead(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);job=root/'job';job.mkdir();check,records=positive_fixture(job,'capped')
            self.assertEqual(records,[])
            initial=root/'initial.tsv';initial.write_text('')
            bank=VerifiedSmoothBank(root/'bank',initial,sha(initial));(job/'VERIFIED').write_text('checked\n')
            with self.assertRaisesRegex(ValueError,'records binding'):bank.publish(job,check,[(9,(), 'false')])
            self.assertEqual(bank.publish(job,check,records)['added'],0)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);job=root/'job';job.mkdir();check,records=positive_fixture(job)
            initial=root/'initial.tsv';initial.write_text('');(job/'VERIFIED').write_text('checked\n')
            bank=VerifiedSmoothBank(root/'bank',initial,sha(initial),limit=1)
            self.assertTrue(bank.publish(job,check,records)['at_bound'])
            self.assertEqual(known_keys(bank.snapshot(9)[0]),set())
            self.assertTrue((job/'LEADS.tsv').exists())

    def test_general_profile_checks_genus_and_invariant_blocks(self):
        from profile_smooth_theory import profile
        c=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
        p=profile(c['representative'],9,3,3)
        self.assertEqual(p['quotient_genus'],1);self.assertEqual(p['group_order'],324)
        self.assertEqual(p['local_index'],2);self.assertTrue(p['relative_genus_checked'])

if __name__=='__main__':unittest.main()
