try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
from smooth_block_lattice import block_lattice
from verify_smooth_geometry import parse_tuple


class BlockLatticeTests(unittest.TestCase):
    def test_regular_elementary_abelian_action_has_all_fourteen_systems(self):
        t=parse_tuple('[(1 2)(3 4)(5 6)(7 8), (1 3)(2 4)(5 7)(6 8), (1 5)(2 6)(3 7)(4 8), (1 8)(2 7)(3 6)(4 5)]',8)
        r=block_lattice(t)
        self.assertTrue(r['complete']);self.assertEqual(r['partitions_seen'],16)
        self.assertEqual(len(r['proper_systems']),14);self.assertEqual(len(r['primitive_quotient_indices']),7)
        self.assertFalse(r['primitive'])
    def test_known_degree_nine_quotient(self):
        t=parse_tuple('[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]',9)
        r=block_lattice(t)
        self.assertTrue(r['complete']);self.assertEqual([(s['blocks'],s['genus']) for s in r['proper_systems']],[(3,1)])
    def test_cap_is_not_primitivity(self):
        t=parse_tuple('[(1 2 3 4), (), (1 4 3 2), ()]',4)
        r=block_lattice(t,cap=1)
        self.assertFalse(r['complete']);self.assertIsNone(r['primitive'])
    def test_prime_transitive_action_is_primitive(self):
        t=parse_tuple('[(1 2 3 4 5), (), (1 5 4 3 2), ()]',5)
        r=block_lattice(t)
        self.assertTrue(r['complete']);self.assertTrue(r['primitive'])


if __name__=='__main__':unittest.main()
