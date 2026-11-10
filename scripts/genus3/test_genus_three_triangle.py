"""Coverage/negative corruption controls; real binary paths supplied explicitly."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from verify_genus_three_triangle import audit
from verify_genus_three_finite import connection

class TriangleTests(unittest.TestCase):
    def test_character_counts_are_exact(self):
        for n,expected in [(5,5),(9,189),(13,19305)]:
            invol=(1,)+(2,)*((n-1)//2)
            self.assertEqual(connection(invol,invol),n)
            self.assertEqual(connection(invol,(n,)),expected)

    @unittest.skipUnless(os.environ.get('G3_TRIANGLE_BINARY')and os.environ.get('G3_TRIANGLE_AUDITOR'),'explicit native fixtures required')
    def test_complete_small_control_and_corruptions(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);base=root/'base';native=Path(os.environ['G3_TRIANGLE_AUDITOR'])
            subprocess.run([os.environ['G3_TRIANGLE_BINARY'],'5',str(base)],check=True,capture_output=True,text=True)
            self.assertEqual(audit(base,native)['coverage']['negative'],5)
            for kind in ['missing','false_word','duplicate']:
                folder=root/kind;shutil.copytree(base,folder)
                records=(folder/'RESOLUTIONS.tsv').read_text().splitlines();leads=(folder/'LEADS.tsv').read_text().splitlines()
                if kind=='missing':records.pop();leads.pop()
                elif kind=='false_word':
                    fields=records[1].split('\t');fields[2]='0';records[1]='\t'.join(fields)
                else:
                    first=records[1].split('\t');last=records[-1].split('\t');last[2]=first[2];last[4]=first[4];records[-1]='\t'.join(last)
                    leads[-1]=last[0]+'\t'+last[4]
                (folder/'RESOLUTIONS.tsv').write_text('\n'.join(records)+'\n');(folder/'LEADS.tsv').write_text('\n'.join(leads)+'\n')
                with self.subTest(kind=kind),self.assertRaises((AssertionError,subprocess.CalledProcessError,ValueError)):
                    audit(folder,native)

if __name__=='__main__':unittest.main()
