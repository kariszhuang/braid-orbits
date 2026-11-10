try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

from audit_guided_native import audit_component_native,audit_guided_native
from run_smooth_mission import audit_component
from test_theory_frontier import positive_fixture
from run_s18_linked_roots import state_key
from analyze_fixed_braid_constraints import product,point_orbits
from odd_order_spin import inverse

PROJECT=Path(__file__).resolve().parents[2]
AUDITOR=Path(os.environ.get('GUIDED_AUDIT_BINARY',PROJECT/'target/theory-native/release/smooth_negative_audit'))


class NativeComponentTests(unittest.TestCase):
    def test_closed_and_partial_match_full_python(self):
        for status in ['closed_smooth','capped']:
            with self.subTest(status=status),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);positive_fixture(root,status)
                self.assertEqual(audit_component_native(root/'component_1',AUDITOR,True),
                    audit_component(root/'component_1',True))
                fast=audit_guided_native(root,AUDITOR,native_components=True);old=audit_guided_native(root,AUDITOR)
                fast['independent_native'].pop('seconds');old['independent_native'].pop('seconds')
                self.assertEqual(fast,old)

    def test_corrupt_rows_edges_and_summary_fail_without_key_export(self):
        edits=[('STATES.tsv',lambda s:s.replace('\ttrue\t','\tfalse\t',1)),
               ('STATES.tsv',lambda s:s.replace('\t0\t2\t','\t1\t2\t',1)),
               ('STATES.tsv',lambda s:s.replace('\t0\t2\t','\t0\t7\t',1)),
               ('STATES.tsv',lambda s:s.replace('0\ttrue','1\ttrue',1)),
               ('STATES.tsv',lambda s:'\n'.join(s.splitlines()[:-1])+'\n'),
               ('SUMMARY.tsv',lambda s:s.replace('\t12\t','\t11\t',1))]
        for name,edit in edits:
            with self.subTest(name=name),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);positive_fixture(root);folder=root/'component_1';p=folder/name;p.write_text(edit(p.read_text()))
                export=root/'keys.tsv';r=subprocess.run([str(AUDITOR),'--component',str(folder),'9',str(export)],capture_output=True,text=True)
                self.assertNotEqual(r.returncode,0);self.assertFalse(export.exists())
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);positive_fixture(root);p=root/'component_1/STATES.tsv'
            lines=p.read_text().splitlines();row=lines[1].split('\t');edge=row[7].split(',');edge[0]=str((int(edge[0])+1)%12)
            row[7]=','.join(edge);lines[1]='\t'.join(row);p.write_text('\n'.join(lines)+'\n')
            with self.assertRaises(subprocess.CalledProcessError):audit_component_native(root/'component_1',AUDITOR)

    def test_larger_literal_keys_equal_python(self):
        rng=random.Random(2026101333)
        for n in [9,14,25]:
            p=subprocess.Popen([str(AUDITOR),'--key-server',str(n)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
            self.assertEqual(p.stdout.readline().strip(),'READY')
            try:
                for _ in range(100):
                    a,b,c=[tuple(rng.sample(range(n),n)) for _ in range(3)];t=(a,b,c,inverse(product(a,b,c)))
                    if len(point_orbits(t))!=1:continue
                    p.stdin.write('\t'.join([*(','.join(map(str,g)) for g in t),''])+'\n');p.stdin.flush()
                    key=tuple(tuple(g) for g in json.loads(p.stdout.readline()))
                    self.assertEqual(key,state_key(t))
            finally:
                p.stdin.close();self.assertEqual(p.wait(timeout=10),0);p.stdout.close()

    def test_graph_action_positive_control(self):
        from profile_smooth_actions import action_profile
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);check,_=positive_fixture(root);sha=check['components'][0]['table_sha256']
            p=action_profile(root/'component_1',sha)
            self.assertEqual(p['pure_braid_orbits'],4);self.assertEqual(p['pure_braid_fiber_size'],3)

if __name__=='__main__':unittest.main()
