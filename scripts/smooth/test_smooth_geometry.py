try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import random
import unittest
from verify_smooth_geometry import analyze, compose, forced_partition, hurwitz, inverse, parse_tuple

DEGREE15 = '''[(1 2 3),
(1 4 7 10 13)(2 5 8 11 14)(3 6 9 12 15),
(1 5 9 12 15)(2 6 8 11 14)(3 4 7 10 13),
(1 12 5 14 8)(2 10 4 15 9)(3 11 6 13 7)]'''
ELLIPTIC12 = '''[(1 3)(2 4),
(1 5 9)(2 6 10)(3 7 11)(4 8 12),
(1 6 12)(2 7 11)(3 8 10)(4 5 9),
(1 7 10)(2 6 11)(3 5 12)(4 8 9)]'''
ELLIPTIC15 = '''[(2 4)(3 5),
(1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15),
(1 7 12)(2 8 15)(3 9 14)(4 10 13)(5 6 11),
(1 6 12)(2 9 13)(3 8 14)(4 7 15)(5 10 11)]'''


class GeometryChecks(unittest.TestCase):
    def test_fixed_higher_genus_quotient(self):
        result = analyze(parse_tuple(DEGREE15, 15), 5, 3)
        self.assertEqual(result['source_genus'], 5)
        self.assertEqual(result['quotient_genus'], 2)
        self.assertEqual(result['quotient_shifts'], [0, 1, 1, 3])
        self.assertEqual(result['relative_branch_value_count'], 1)
        self.assertEqual(result['relative_branching_index'], 2)
        self.assertEqual(result['tangent_pullback_degree'], -6)
        self.assertEqual(result['geometry'], 'nonisotrivial_higher_genus')

    def test_elliptic_not_nonisotrivial(self):
        result = analyze(parse_tuple(ELLIPTIC12, 12), 3, 4)
        self.assertEqual(result['source_genus'], 2)
        self.assertEqual(result['quotient_genus'], 1)
        self.assertEqual(result['local_cycle_lengths'], [2, 2])
        self.assertEqual(result['relative_branch_value_count'], 1)
        self.assertEqual(result['geometry'], 'isotrivial_elliptic')
        self.assertEqual(result['tangent_pullback_degree'], 0)

    def test_partial_block_support(self):
        result = analyze(parse_tuple(ELLIPTIC15, 15), 3, 5)
        self.assertEqual(result['block_partition'][0], [1, 2, 3, 4, 5])
        self.assertEqual(result['local_cycle_lengths'], [2, 2])
        self.assertEqual(result['source_genus'], 2)

    def test_all_six_moves_and_relabeling(self):
        t = parse_tuple(DEGREE15, 15)
        rng = random.Random(139)
        h = list(range(15))
        for _ in range(30):
            t = hurwitz(t, rng.randrange(6))
            rng.shuffle(h)
            relabeled = tuple(compose(compose(h, p), inverse(h)) for p in t)
            for m in range(6):
                moved = hurwitz(relabeled, m)
                self.assertEqual(hurwitz(moved, (m+3) % 6), relabeled)
                result = analyze(moved, 5, 3)
                self.assertEqual(result['geometry'], 'nonisotrivial_higher_genus')
                self.assertEqual(result['relative_branch_value_count'], 1)

    def test_bad_cycle_notation_rejected(self):
        for text in ['[(1 2)(2 3), (), (), ()]', '[(16), (), (), ()]',
                     '[(1 1), (), (), ()]', '[(1 2)x, (), (), ()]']:
            with self.assertRaises(ValueError):
                parse_tuple(text, 15)

    def test_bad_product_rejected(self):
        t = list(parse_tuple(DEGREE15, 15))
        t[0] = inverse(t[0])
        with self.assertRaisesRegex(ValueError, 'product is not one'):
            analyze(tuple(t), 5, 3)

    def test_wrong_block_parameters_rejected(self):
        with self.assertRaises(ValueError):
            analyze(parse_tuple(DEGREE15, 15), 3, 5)

    def test_extra_relative_branching_rejected(self):
        # An r*k full cycle has unit block shift but nontrivial return map.
        # Recompute the fourth entry, preserving product one; still must fail
        # the transversal ramification test rather than be called certified.
        t = list(parse_tuple(DEGREE15, 15))
        twist = parse_tuple('[(1 2 3), (), (), ()]', 15)[0]
        t[1] = compose(twist, t[1])
        t[3] = inverse(compose(compose(t[0], t[1]), t[2]))
        with self.assertRaisesRegex(ValueError, r'large cycle type not r\^k'):
            analyze(tuple(t), 5, 3)

    def test_forced_partition_equivariance(self):
        t = parse_tuple(ELLIPTIC15, 15)
        parts = forced_partition(t, [[1, 2, 3, 4]])
        membership = {x:i for i,part in enumerate(parts) for x in part}
        self.assertEqual(sorted(map(len, parts)), [5, 5, 5])
        for p in t:
            self.assertTrue(all(len({membership[p[x]] for x in part}) == 1 for part in parts))


if __name__ == '__main__':
    unittest.main()
