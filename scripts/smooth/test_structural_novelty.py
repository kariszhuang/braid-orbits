try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from run_forest_fast import novelty_job,triangle_target_job,affine_interval_cover
from audit_guided_native import audit_guided_native
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]
BINARY=Path(os.environ.get('DESIGN_BINARY',PROJECT/'target/release/smooth_forest_fast'))
AUDITOR=Path(os.environ.get('GUIDED_AUDIT_BINARY',PROJECT/'target/release/smooth_negative_audit'))


class StructuralNoveltyTests(unittest.TestCase):
    def test_triangle_schedule_skips_reported_exclusions_and_covers_higher_r1_bands(self):
        spec=dict(minimum_degree=16,maximum_degree=63,passport_minimum_shards=[
            dict(passport_min=7,passport_max=7,minimum_degree=16,maximum_degree=16,
                 pair_index_min=7,pair_index_max=7,strong_divisor=8,exact_genus=3,
                 exact_passport=[7,7,7,15]),
            dict(passport_min=6,passport_max=12,minimum_degree=17,maximum_degree=63,
                 pair_index_min=6,pair_index_max=63,strong_divisor=8,exact_genus=3),
            dict(passport_min=13,passport_max=31,minimum_degree=27,maximum_degree=63,
                 pair_index_min=13,pair_index_max=63,strong_divisor=8,exact_genus=3)])
        shards=spec['passport_minimum_shards']
        for d in range(16,64):
            expected=([7] if d==16 else range(6,(d-1)//2+1))
            for r1 in expected:
                owners=[i for i,b in enumerate(shards) if b['passport_min']<=r1<=b['passport_max']
                        and b['minimum_degree']<=d<=b['maximum_degree']]
                self.assertEqual(len(owners),1,(d,r1,owners))
                n,lane,reuse,group,band=triangle_target_job(owners[0],0,spec)
                self.assertEqual(lane,'triangle');self.assertEqual(reuse,1);self.assertIsNone(group)
                self.assertIn(band['third_method'],['random','joint','mixed','strong','centralizer','pairwise_low','middle_repair'])
                self.assertEqual(band['exact_genus'],3)
                self.assertEqual(band['pair_index_min'],band['passport_min'])
                self.assertGreaterEqual(n,max(2*band['passport_min']+1,band['minimum_degree']))
                self.assertLessEqual(n,63)
                if owners[0]==0:self.assertEqual((n,band['exact_passport']),(16,[7,7,7,15]))
        self.assertEqual(triangle_target_job(0,0,spec)[0],16)
        self.assertEqual(triangle_target_job(0,688,spec)[4]['third_method'],'middle_repair')
        self.assertEqual(triangle_target_job(0,695,spec)[4]['third_method'],'middle_repair')

    def test_affine_coverage_requires_every_independent_interval(self):
        def proof(a,b,unresolved=False):
            return dict(passed=True,start=a,next=b,total=20,no_unresolved=not unresolved,
                        counts=dict(qualifying=b-a,literal_negative=b-a,accounted_survivor=0))
        self.assertFalse(affine_interval_cover([proof(0,10)],20)['negative_census'])
        self.assertTrue(affine_interval_cover([proof(0,10),proof(10,20)],20)['negative_census'])
        self.assertFalse(affine_interval_cover([proof(0,20,True)],20)['negative_census'])
        for parts in [[proof(0,10),proof(11,20)],[proof(0,10),proof(9,20)]]:
            with self.assertRaises(ValueError):affine_interval_cover(parts,20)

    def test_schedule_has_new_families_and_no_completed_degrees(self):
        for worker in range(6):
            for epoch in range(32):
                n,lane,reuse,group=novelty_job(worker,epoch,{})
                self.assertGreaterEqual(n,11)
                if worker<2:self.assertEqual(lane,'triangle')
                if worker==4:self.assertEqual((n,lane,group),(16,'primitive','asp16'))
        self.assertEqual(novelty_job(4,0,{'asp16_16':{'producer_complete':True}})[3],'agl16')

    def test_cached_primitive_frames_preserve_complete_journals_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name,flag in [('cached',[]),('uncached',['--uncached-primitive-frames'])]:
                subprocess.run([str(BINARY),'--d','16','--lane','primitive','--primitive-group','asp16','--frame-reuse','1',
                    '--primitive-cursor','150000','--seconds','30','--attempts','10000','--out',str(root/name),*flag],check=True,capture_output=True)
                audit_guided_native(root/name,AUDITOR,native_components=True)
            for name in ['LEADS.tsv','RESOLUTIONS.tsv','KNOWN.tsv','PAIR_FRAMES.tsv','GROUP.tsv']:
                self.assertEqual(sha(root/'cached'/name),sha(root/'uncached'/name))
            def count(path):
                with path.open() as f:r=next(csv.DictReader(f,delimiter='\t'));r.pop('seconds');return r
            self.assertEqual(count(root/'cached/COUNTS.tsv'),count(root/'uncached/COUNTS.tsv'))

    def test_native_omission_audit_detects_deleted_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'primitive'
            subprocess.run([str(BINARY),'--d','16','--lane','primitive','--primitive-group','asp16','--frame-reuse','1',
                '--seconds','20','--attempts','100','--out',str(folder)],check=True,capture_output=True)
            subprocess.run([str(AUDITOR),'--primitive-population',str(folder)],check=True,capture_output=True)
            p=folder/'PAIR_FRAMES.tsv';p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n')
            self.assertNotEqual(subprocess.run([str(AUDITOR),'--primitive-population',str(folder)],capture_output=True).returncode,0)

    def test_new_constructor_jobs_and_independent_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for lane in ['triangle','returns']:
                folder=root/lane
                subprocess.run([str(BINARY),'--d','12','--lane',lane,'--frame-reuse','1','--seconds','20',
                    '--attempts','2000','--full-closure-budget','--out',str(folder)],check=True,capture_output=True)
                r=audit_guided_native(folder,AUDITOR,native_components=True)
                self.assertEqual(int(r['construction']['constructed']),2000)


if __name__=='__main__':unittest.main()
