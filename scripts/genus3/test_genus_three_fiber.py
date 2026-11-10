try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest

from genus_three_fiber_geometry import full_fiber_witness,partial_fiber_witness,translation_decks
from verify_smooth_geometry import parse_tuple


class FiberGeometryTests(unittest.TestCase):
    def test_elliptic_quotient_translation_decks_use_inertia_fixed_points(self):
        ident=(0,1,2);cycle=(1,2,0)
        self.assertEqual(translation_decks([ident,cycle,cycle,cycle]),[])
        # An elliptic degree-four V4 cover: a,b have fixed points upstairs,
        # but their product is a nonidentity translation.
        a=(1,0,3,2);b=(2,3,0,1)
        self.assertEqual(translation_decks([a,a,b,b]),[(3,2,1,0)])

    def test_partial_fibre_argument_is_separate_from_closure(self):
        t=parse_tuple('[(1 2 3)(8 9)(14 15), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 13)(2 8 14)(3 10 15)(4 6 11)(5 7 12), (1 9 15)(2 6 13)(3 7 14)(4 8 11)(5 10 12)]',15)
        blocks=[list(range(1,6)),list(range(6,11)),list(range(11,16))]
        result=partial_fiber_witness(t,blocks,0)
        self.assertEqual(result['active_blocks'],3);self.assertNotIn('all_smooth',result)
        self.assertIsNone(partial_fiber_witness(t,blocks,1))
    def test_full_fibre_has_literal_hypotheses_without_claiming_smoothness(self):
        t = parse_tuple('[(1 2 3)(8 9)(14 15), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 13)(2 8 14)(3 10 15)(4 6 11)(5 7 12), (1 9 15)(2 6 13)(3 7 14)(4 8 11)(5 10 12)]', 15)
        blocks = [list(range(1,6)), list(range(6,11)), list(range(11,16))]
        result = full_fiber_witness(t, blocks, 0)
        self.assertEqual(result['relative_excess'], [4,0,0,0])
        self.assertEqual(result['active_blocks'], 3)
        self.assertNotIn('all_smooth', result)
        self.assertIsNone(full_fiber_witness(t, blocks, 1))

    def test_old_single_local_isotrivial_family_is_not_full_fibre(self):
        t = parse_tuple('[(1 3 4 2 5), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 12)(2 8 15)(3 10 14)(4 7 11)(5 6 13), (1 8 14)(2 10 13)(3 7 15)(4 6 12)(5 9 11)]', 15)
        blocks = [list(range(1,6)), list(range(6,11)), list(range(11,16))]
        self.assertIsNone(full_fiber_witness(t, blocks, 0))


if __name__ == '__main__':
    unittest.main()
