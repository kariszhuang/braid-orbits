try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
from run_smooth_noncyclic import (finite_population, direct_population, analyze,
                                 compose, inverse, hurwitz, BETA, GAMMA, orbit)


class NoncyclicChecks(unittest.TestCase):
    def test_literal_exhaustive_same_family_comparison(self):
        population = finite_population()
        examined, direct = direct_population()
        self.assertEqual(examined, 46656)
        self.assertEqual(len(population), 1296)
        self.assertEqual(len(set(population)), 1296)
        self.assertEqual(set(population), direct)

    def test_every_qualifying_occurrence_and_six_neighbors(self):
        accepted = 0
        for t in finite_population():
            if t[0] == tuple(range(15)):
                continue
            profile = analyze(t)
            accepted += 1
            self.assertEqual(profile['source_genus'], 5)
            self.assertEqual(profile['quotient_genus'], 2)
            for m in range(6):
                self.assertTrue(analyze(hurwitz(t, m))['smooth'])
        self.assertGreater(accepted, 0)

    def test_relabeling_and_negative_control(self):
        t = next(t for t in finite_population() if t[0] != tuple(range(15)))
        h = tuple(reversed(range(15)))
        relabeled = tuple(compose(compose(h, p), inverse(h)) for p in t)
        self.assertEqual(analyze(relabeled)['source_genus'], 5)
        with self.assertRaises(ValueError):
            analyze((tuple(range(15)), *t[1:]))

    def test_noncyclic_quotient(self):
        self.assertNotEqual(compose(BETA, GAMMA), compose(GAMMA, BETA))
        self.assertEqual(len(orbit(compose(BETA, GAMMA))), 5)
        from sympy.combinatorics import Permutation, PermutationGroup
        self.assertEqual(PermutationGroup([Permutation(list(BETA)), Permutation(list(GAMMA))]).order(), 60)


if __name__ == '__main__':
    unittest.main()
