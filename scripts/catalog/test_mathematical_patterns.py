"""Local controls for the report audit, not a verification of its geometry."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import copy
import json
from itertools import permutations
import unittest

from audit_mathematical_patterns import PROJECT, compose, cycles, inspect, orbit_sizes, permutation, quotient_check


class PatternAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((PROJECT/'exports/latest-discoveries-20261014/release-v5/all-smooth.json').read_text())['components']
        cls.old = next(e for e in cls.data if e['degree'] == 9)
        cls.ramified = next(e for e in cls.data if e['degree'] == 18 and e['group_order'] == 41472 and e['states'] == 1296)

    def test_original_elliptic_control(self):
        row = inspect(self.old)
        self.assertIs(row['original_isotrivial'], True)
        self.assertEqual(row['elliptic_fixed_and_moving_witnesses'], [])

    def test_ramified_elliptic_witness(self):
        row = inspect(self.ramified)
        self.assertIsNone(row['original_isotrivial'])
        w, = row['elliptic_fixed_and_moving_witnesses']
        self.assertEqual(w['excess'], [4, 0, 1, 1])
        self.assertEqual(w['fixed_relative_branch_points'], 2)
        self.assertEqual(w['moving_relative_branch_points'], 1)
        self.assertEqual(w['relative_branch_local_degrees'][2][0]['local_degrees'], [2, 1, 1, 1, 1])
        self.assertEqual(len(w['blocks']), 3)
        self.assertEqual(w['quotient_tuple'], ['()', '(1 2 3)', '(1 2 3)', '(1 2 3)'])

    def test_identity_corruption_rejected(self):
        e = copy.deepcopy(self.old)
        e['tuple'][0] = '()'
        with self.assertRaises(ValueError):
            inspect(e)

    def test_excess_corruption_rejected(self):
        e = copy.deepcopy(self.ramified)
        e['block_quotients'][0]['ramification_excess'][0] += 1
        with self.assertRaises(ValueError):
            inspect(e)

    def test_noninvariant_partition_rejected(self):
        e = self.old
        b = copy.deepcopy(e['block_quotients'][0]['blocks'])
        b[0][0], b[1][0] = b[1][0], b[0][0]
        with self.assertRaises(ValueError):
            quotient_check([permutation(s, e['degree']) for s in e['tuple']], b)

    def test_overlapping_cycles_rejected(self):
        with self.assertRaises(ValueError):
            permutation('(1 2)(2 3)', 3)

    def test_multi_active_hyperbolic_control(self):
        e = next(e for e in self.data if e['id'] == '670786f56741a481d24d5250d90af2426f85aad2fe62d29bcd87e7c68e60321c')
        row = inspect(e)
        self.assertIsNone(row['original_isotrivial'])
        w, = row['hyperbolic_moving_witnesses']
        self.assertEqual((w['genus'], w['block_count'], w['block_size']), (2, 12, 2))
        self.assertEqual(len(w['active_blocks']), 2)
        self.assertEqual(w['excess'], [2, 0, 0, 0])

    def test_additive_pair_forests_all_s4_pairs(self):
        ps = list(permutations(range(4)))
        for a in ps:
            for b in ps:
                cc = len(cycles(a))+len(cycles(b))
                components = len(orbit_sizes([a, b]))
                is_forest = cc-4 == components
                is_additive = 4-len(cycles(compose(a, b))) == 8-cc
                self.assertEqual(is_forest, is_additive)

    def test_smooth_complement_lemma_all_s4_triples(self):
        ps = list(permutations(range(4)))
        for a in ps:
            for b in ps:
                ab = compose(a, b)
                iab = 4-len(cycles(ab))
                first = iab == 8-len(cycles(a))-len(cycles(b))
                for c in ps:
                    abc = compose(ab, c)
                    d = tuple(abc.index(x) for x in range(4))
                    second = iab == 8-len(cycles(c))-len(cycles(d))
                    if not (first or second) or orbit_sizes([a, b, c]) != [4]:
                        continue
                    left, right = len(orbit_sizes([a, b])), len(orbit_sizes([c, d]))
                    self.assertEqual(len(cycles(ab))-left-right+1, 0)
                    if first:
                        self.assertEqual(right, 1)
                    if second:
                        self.assertEqual(left, 1)

    def test_end_moves_preserve_all_s4_split_statistics(self):
        ps = list(permutations(range(4)))
        def inverse(p):
            return tuple(p.index(x) for x in range(4))
        def stats(t):
            ab = compose(t[0], t[1])
            left, right = len(orbit_sizes(t[:2])), len(orbit_sizes(t[2:]))
            indices = [4-len(cycles(p)) for p in t]
            smooth = 4-len(cycles(ab)) in (sum(indices[:2]), sum(indices[2:]))
            return ab, left, right, sum(indices[:2]), sum(indices[2:]), smooth
        for a in ps:
            for b in ps:
                for c in ps:
                    d = inverse(compose(compose(a, b), c))
                    t = [a, b, c, d]
                    expected = stats(t)
                    for i in (0, 2):
                        for positive in (True, False):
                            u = list(t)
                            x, y = u[i:i+2]
                            u[i:i+2] = [compose(compose(x, y), inverse(x)), x] if positive else [y, compose(compose(inverse(y), x), y)]
                            self.assertEqual(stats(u), expected)


if __name__ == '__main__':
    unittest.main()
