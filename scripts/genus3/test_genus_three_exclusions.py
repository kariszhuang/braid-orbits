try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from pathlib import Path
import tempfile
import unittest

from verify_genus_three_exclusions import audit_exclusions

SEED='[(1 3 4 2 5), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 12)(2 8 15)(3 10 14)(4 7 11)(5 6 13), (1 8 14)(2 10 13)(3 7 15)(4 6 12)(5 9 11)]'


class ExclusionTests(unittest.TestCase):
    def fixture(self,path,special=1,labels=None):
        if labels is None:labels=[i//5 for i in range(15)]
        path.write_text('degree\tcandidate\tspecial\tlabels\tseed\n15\t1\t'+str(special)+'\t'+','.join(map(str,labels))+'\t'+SEED+'\n')

    def test_valid_static_criterion_is_not_claimed_as_a_closure(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'proof.tsv';self.fixture(p)
            v=audit_exclusions(p,1)
            self.assertEqual(v['records'],1)
            self.assertIn('not a full-orbit',v['scope'])

    def test_wrong_special_wrong_partition_and_false_counts_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'proof.tsv';self.fixture(p,special=2)
            with self.assertRaises(ValueError):audit_exclusions(p)
            self.fixture(p,labels=[i%3 for i in range(15)])
            with self.assertRaises(ValueError):audit_exclusions(p)
            self.fixture(p)
            with self.assertRaisesRegex(ValueError,'counter'):audit_exclusions(p,2)


if __name__=='__main__':unittest.main()
