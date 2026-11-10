try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest

from odd_order_spin import parse_seed, inverse, hurwitz, cycle_lengths
from analyze_fixed_braid_constraints import product, point_orbits
from run_smooth_frontier import audit_guided, values
from audit_guided_native import audit_guided_native

PROJECT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('GUIDED_AUDIT_BINARY', PROJECT / 'target/forest-second/release/smooth_negative_audit'))


def seed_text():
    # Independent conversion avoids requiring a hand-computed fourth entry.
    a,b,c,_ = parse_seed('[(1 2), (2 3)(4 5)(6 7), (1 2 3 4 5 6 7 8 9), ()]', 9)
    def cycle_string(p):
        seen=set();out=[]
        for x in range(len(p)):
            if x in seen:continue
            c=[]
            while x not in seen:seen.add(x);c.append(x+1);x=p[x]
            if len(c)>1:out.append('('+' '.join(map(str,c))+')')
        return ''.join(out) or '()'
    return '[' + ', '.join(cycle_string(p) for p in (a,b,c,inverse(product(a,b,c)))) + ']'


def fixture(folder, status='nonsmooth', word='4'):
    seed = seed_text()
    name = '' if status == 'nonsmooth' else 'component_1'
    (folder/'LEADS.tsv').write_text(f'1\t{seed}\n')
    (folder/'RESOLUTIONS.tsv').write_text(f'candidate\tstatus\tword\tfolder\tseed\n1\t{status}\t{word}\t{name}\t{seed}\n')
    closed = int(status == 'closed_smooth')
    partial = int(status in ['capped', 'timeout'])
    (folder/'COUNTS.tsv').write_text('degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'
                                   f'9\treplay\t1\t1\t1\t1\t1\t0\t0\t{closed}\t{partial}\t0.1\n')


class GuidedNativeTests(unittest.TestCase):
    def native(self, root):
        return subprocess.run([str(BINARY), '--guided', str(root), '9'], capture_output=True, text=True)

    def test_matches_python_and_checks_every_record(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fixture(root)
            old=audit_guided(root);new=audit_guided_native(root,BINARY)
            self.assertEqual(old['construction'],new['construction'])
            self.assertEqual(old['components'],new['components'])
            self.assertEqual(new['independent_native']['records'],1)
            self.assertEqual(new['independent_native']['nonsmooth'],1)

    def test_corruption_and_false_word_fail_closed(self):
        edits=[('RESOLUTIONS.tsv',lambda s:s.replace('\tnonsmooth\t4\t','\tnonsmooth\t0,3\t')),
               ('RESOLUTIONS.tsv',lambda s:s.replace('\tnonsmooth\t4\t','\tnonsmooth\t6\t')),
               ('RESOLUTIONS.tsv',lambda s:s.replace('\tnonsmooth\t4\t','\tnonsmooth\t\t')),
               ('LEADS.tsv',lambda s:s+'2\t'+seed_text()+'\n'),
               ('LEADS.tsv',lambda s:s.replace('1\t','2\t',1)),
               ('LEADS.tsv',lambda s:s[:-15]),
               ('COUNTS.tsv',lambda s:s.replace('9\treplay','8\treplay')),
               ('RESOLUTIONS.tsv',lambda s:s.replace('(1 2)', '(1 26)'))]
        for filename, edit in edits:
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as temp:
                root=Path(temp);fixture(root);p=root/filename;p.write_text(edit(p.read_text()))
                self.assertNotEqual(self.native(root).returncode,0)

    def test_partial_is_not_negative_and_still_requires_component_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fixture(root,'capped','')
            r=self.native(root);self.assertEqual(r.returncode,0,r.stderr)
            data=json.loads(r.stdout)
            self.assertEqual(data['nonsmooth'],0);self.assertEqual(data['unresolved'],1)
            self.assertEqual(data['exceptions'][0]['status'],'capped')
            with self.assertRaises(FileNotFoundError):audit_guided_native(root,BINARY)
            p=root/'RESOLUTIONS.tsv';p.write_text(p.read_text().replace('component_1','../unsafe'))
            self.assertNotEqual(self.native(root).returncode,0)

    def test_real_positive_and_partial_components_use_full_python_checks(self):
        catalog=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components'][0]
        for status in ['closed_smooth','capped']:
            with self.subTest(status=status), tempfile.TemporaryDirectory() as temp:
                root=Path(temp);fixture(root,status,'')
                for name in ['LEADS.tsv','RESOLUTIONS.tsv']:
                    p=root/name;p.write_text(p.read_text().replace(seed_text(),catalog['representative']))
                component=root/'component_1';component.mkdir()
                for name in ['SUMMARY.tsv','STATES.tsv']:
                    shutil.copy2(Path(catalog['raw_folder'])/name,component/name)
                if status=='capped':
                    p=component/'SUMMARY.tsv';p.write_text(p.read_text().replace('closed_smooth','capped'))
                check=audit_guided_native(root,BINARY)
                self.assertEqual(len(check['components']),1)
                self.assertEqual(check['components'][0]['independent_six_move_closure'],status=='closed_smooth')
                self.assertEqual(check['independent_native']['nonsmooth'],0)
                p=component/'STATES.tsv';p.write_text(p.read_text().replace('\ttrue\t','\tfalse\t',1))
                with self.assertRaisesRegex(ValueError,'smoothness mismatch'):audit_guided_native(root,BINARY)

    def test_literal_kernel_larger_degrees_and_positive_catalog_parser(self):
        rng=random.Random(20261012)
        for n in [11,14,25]:
            p=subprocess.Popen([str(BINARY),'--literal-server',str(n)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
            self.assertEqual(p.stdout.readline().strip(),'READY')
            try:
                for _ in range(100):
                    a,b,c=[tuple(rng.sample(range(n),n)) for _ in range(3)]
                    t=(a,b,c,inverse(product(a,b,c)));word=[rng.randrange(6) for _ in range(7)];moved=t
                    for m in word:moved=hurwitz(moved,m)
                    p.stdin.write('\t'.join([*(','.join(map(str,g)) for g in t),','.join(map(str,word))])+'\n');p.stdin.flush()
                    data=json.loads(p.stdout.readline())
                    self.assertEqual(data['tuple'],[list(g) for g in moved])
                    self.assertEqual(data['indices'],[n-len(cycle_lengths(g)) for g in moved])
                    self.assertEqual(data['pair_orbits'],[len(point_orbits(moved[i:i+2])) for i in range(3)])
                    self.assertEqual(data['smooth'],values(moved)['smooth'])
            finally:
                p.stdin.close();self.assertEqual(p.wait(timeout=10),0);p.stdout.close()
        catalog=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components']
        for c in catalog:
            with tempfile.TemporaryDirectory() as temp:
                root=Path(temp);fixture(root,'closed_smooth','')
                for name in ['LEADS.tsv','RESOLUTIONS.tsv']:
                    p=root/name;p.write_text(p.read_text().replace(seed_text(),c['representative']))
                p=root/'COUNTS.tsv';p.write_text(p.read_text().replace('9\treplay',f"{c['degree']}\treplay"))
                r=subprocess.run([str(BINARY),'--guided',str(root),str(c['degree'])],capture_output=True,text=True)
                self.assertEqual(r.returncode,0,r.stderr)
                self.assertEqual(json.loads(r.stdout)['closed'],1)


if __name__=='__main__':unittest.main()
