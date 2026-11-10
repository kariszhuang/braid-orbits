try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from audit_guided_native import audit_guided_native,prime_known_cache,known_keys
from run_forest_fast import mechanism_job
from run_s18_linked_roots import state_key
from odd_order_spin import parse_seed
from profile_mechanism_breaker import primitive,projection_profile,install_bases,ramified_double_profile
from audit_primitive_family import check

PROJECT=Path(__file__).resolve().parents[2]
BINARY=Path(os.environ.get('DESIGN_BINARY',PROJECT/'target/mechanism-frontier/release/smooth_forest_fast'))
AUDITOR=Path(os.environ.get('GUIDED_AUDIT_BINARY',PROJECT/'target/mechanism-frontier/release/smooth_negative_audit'))


class MechanismTests(unittest.TestCase):
    def test_schedule_keeps_six_slots_and_retires_only_finite_constructor(self):
        self.assertEqual(mechanism_job(0,9,{})[:2],(11,'primitive'))
        self.assertEqual(mechanism_job(2,0,{})[:2],(12,'mutate'))
        self.assertEqual(mechanism_job(3,0,{})[:2],(11,'forest'))
        self.assertEqual(mechanism_job(4,0,{})[:2],(24,'double_lift'))
        self.assertEqual(mechanism_job(4,1,{})[1],'quotient')
        self.assertEqual(mechanism_job(4,0,{}, {'24':dict(producer_complete=True)})[1],'quotient')
        with self.assertRaises(ValueError):mechanism_job(5,0,{})
        cursor={'psl11_11':dict(producer_complete=True)}
        self.assertEqual(mechanism_job(0,0,cursor)[3],'psl2_16')

    def test_exact_small_group_family_and_slice_restart(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for start,name in [(0,'full'),(1000,'slice')]:
                folder=root/name
                subprocess.run([str(BINARY),'--d','11','--lane','primitive','--primitive-group','psl11',
                    '--frame-reuse','1','--primitive-cursor',str(start),'--seconds','20','--out',str(folder)],check=True,capture_output=True)
                report=audit_guided_native(folder,AUDITOR,native_components=True)
                exact=check(folder,root/f'{name}.json')
                self.assertTrue(exact['passed']);self.assertEqual(exact['start'],start)
                self.assertEqual(int(report['construction']['constructed']),exact['next']-start)
                self.assertEqual(exact['frames'],8);self.assertEqual(exact['total'],5280)

    def test_mutation_first_pair_product_and_independent_primitive_labels(self):
        e=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
        # Supplied degree-nine control is used ONLY as a bounded algorithm test.
        blocks=[2,0,0,0,1,1,1,2,2]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seed=root/'seed.tsv';seed.write_text('9\t'+e['representative']+'\t'+','.join(map(str,blocks))+'\n')
            job=root/'job'
            subprocess.run([str(BINARY),'--d','9','--lane','mutate','--mutation-seeds',str(seed),'--frame-reuse','1',
                '--seed','20261013','--seconds','10','--attempts','3000','--out',str(job)],check=True,capture_output=True)
            audit_guided_native(job,AUDITOR,native_components=True)
            from run_smooth_degree9 import rows
            leads=dict(line.split('\t',1) for line in (job/'LEADS.tsv').read_text().splitlines())
            for r in rows(job/'MECHANISM.tsv'):
                t=parse_seed(leads[r['candidate']],9)
                self.assertEqual(primitive(t),r['primitive']=='true')

    def test_native_positive_key_cache_matches_python_and_rejects_corruption(self):
        e=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'bank.tsv'
            from run_smooth_degree9 import rows
            seeds=[r['tuple'] for r in rows(Path(e['raw_folder'])/'STATES.tsv')]
            p.write_text(''.join('9\t'+s+'\n' for s in seeds))
            self.assertEqual(prime_known_cache(p,AUDITOR)['records'],12)
            self.assertEqual(known_keys(p),{(9,state_key(parse_seed(s,9))) for s in seeds})
            p.write_text('9\t[(), (), (), ()]\n')
            with self.assertRaises(subprocess.CalledProcessError):prime_known_cache(p,AUDITOR)

    def test_marked_quotient_projection_requires_extra_lift_data(self):
        e=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
        q=projection_profile(e)
        self.assertEqual((q['pure_fiber'],q['marked_quotient_classes']),(3,1))
        self.assertFalse(q['faithful']);self.assertTrue(q['deterministic_factor'])

    def test_ramified_partition_uses_the_displayed_state_labels(self):
        base=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
        item=json.loads((PROJECT/'evidence/smooth_mechanism_breaker_20261013/CATALOG.json').read_text())['components'][0]
        install_bases([base])
        result=ramified_double_profile(item)
        t=parse_seed(result['representative'],18)
        blocks=[[x-1 for x in b] for b in result['block_partition']]
        membership={x:j for j,b in enumerate(blocks) for x in b}
        for p in t:
            self.assertTrue(all(len({membership[p[x]] for x in b})==1 for b in blocks))
        self.assertEqual(result['partition_for'],'STATES.tsv state0')
        self.assertEqual(result['constructor_seed'],item['representative'])
        self.assertEqual((result['genus'],result['base_genus']),(4,2))

if __name__=='__main__':unittest.main()
