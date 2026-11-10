"""Arithmetic theorem-target skip controls, no synthetic braid certificates."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from pathlib import Path
import tempfile
import unittest
from verify_genus_three_transpositions import verify_theorem_exclusions


class TranspositionExclusionTests(unittest.TestCase):
    def fixture(self,root):
        t='[(8 9), (1 12 7 3 2)(4 5 6 11 10), (1 8 11 3 10)(2 4 9 12 5)(6 7), (1 5 2 11 7 9)(3 6 12 10 8 4)]'
        (root/'COUNTS.tsv').write_text('degree\tconstructed\n12\t99\n')
        (root/'TARGET_THEOREM.tsv').write_text('theorem\texcluded_occurrences\tsampled_occurrences\tavoided_sorted_profiles\nno-positive-genus-all-smooth-transposition\t1\t1\t5\n')
        (root/'TARGET_EXCLUSIONS.tsv').write_text('candidate\tseed\n3\t'+t+'\n')

    def test_exact_index_one_scope_and_zero_sample(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.fixture(root)
            self.assertTrue(verify_theorem_exclusions(root)['passed'])
            (root/'TARGET_EXCLUSIONS.tsv').write_text('candidate\tseed\n')
            (root/'TARGET_THEOREM.tsv').write_text('theorem\texcluded_occurrences\tsampled_occurrences\tavoided_sorted_profiles\nno-positive-genus-all-smooth-transposition\t0\t0\t5\n')
            self.assertEqual(verify_theorem_exclusions(root)['sampled_occurrences'],0)

    def test_wrong_tuple_scope_and_sample_mass_fail(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);self.fixture(root);p=root/'TARGET_EXCLUSIONS.tsv'
            p.write_text(p.read_text().replace('(8 9)','()'))
            with self.assertRaisesRegex(ValueError,'genus-three'):verify_theorem_exclusions(root)
            self.fixture(root);p=root/'TARGET_THEOREM.tsv';p.write_text(p.read_text().replace('\t1\t1\t5','\t0\t1\t5'))
            with self.assertRaisesRegex(ValueError,'metadata'):verify_theorem_exclusions(root)


if __name__=='__main__':unittest.main()
