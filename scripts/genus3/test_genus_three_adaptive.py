"""Literal prefix-template provenance and deliberate corruption controls."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from pathlib import Path
import tempfile
import unittest

from verify_genus_three_adaptive import verify_adaptive
from verify_smooth_geometry import parse_tuple,hurwitz,cycles


class AdaptiveTests(unittest.TestCase):
    def fixture(self,root):
        text='[(8 9), (1 12 7 3 2)(4 5 6 11 10), (1 8 11 3 10)(2 4 9 12 5)(6 7), (1 5 2 11 7 9)(3 6 12 10 8 4)]'
        word=[0,2,4,0,0,0,0,0,0,4,4]
        t=parse_tuple(text,12)
        for m in word[:10]:t=hurwitz(t,m)
        def notation(p):return ''.join('('+' '.join(str(x+1)for x in c)+')'for c in cycles(p)if len(c)>1)or'()'
        seed='['+', '.join(notation(p)for p in t)+']'
        (root/'COUNTS.tsv').write_text('degree\tlane\n12\trefactor\n')
        (root/'ADAPTIVE_COUNTS.tsv').write_text('capacity\tadmitted\tadaptive_frame_draws\n4\t1\t5\n')
        (root/'RESOLUTIONS.tsv').write_text('candidate\tstatus\tword\tseed\n17\tnonsmooth\t'+','.join(map(str,word))+'\t'+text+'\n')
        (root/'ADAPTIVE.tsv').write_text('candidate\tprefix\tseed\n17\t10\t'+seed+'\n')

    def test_literal_prefix_and_capacity(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.fixture(root)
            check=verify_adaptive(root)
            self.assertTrue(check['passed']);self.assertEqual(check['templates'],1)
            with self.assertRaisesRegex(TimeoutError,'deadline'):verify_adaptive(root,deadline_epoch=0)
            (root/'ADAPTIVE_COUNTS.tsv').write_text('capacity\tadmitted\tadaptive_frame_draws\n2048\t1\t5\n')
            with self.assertRaisesRegex(ValueError,'capacity'):verify_adaptive(root)

    def test_corrupt_tuple_prefix_and_root_are_rejected(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp)
            for before,after in (('\t10\t','\t9\t'),('(1','(12'),('17\t10','18\t10')):
                self.fixture(root);p=root/'ADAPTIVE.tsv';p.write_text(p.read_text().replace(before,after,1))
                with self.assertRaises((ValueError,AssertionError)):verify_adaptive(root)
            self.fixture(root);p=root/'RESOLUTIONS.tsv';p.write_text(p.read_text().replace('nonsmooth','closed_smooth'))
            with self.assertRaisesRegex(ValueError,'negative witness'):verify_adaptive(root)


if __name__=='__main__':unittest.main()
