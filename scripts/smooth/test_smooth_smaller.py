try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
from run_smooth_smaller import (finite_population, direct_population, analyze,
                               compose, inverse, hurwitz, orbits)


class SmallerChecks(unittest.TestCase):
    def test_exhaustive_direct_comparison(self):
        population = finite_population()
        examined, direct = direct_population()
        self.assertEqual(examined, 331776)
        self.assertEqual(len(population), 576)
        self.assertEqual(len(set(population)), 576)
        self.assertEqual(set(population), direct)

    def test_every_accepted_frame_and_six_neighbors(self):
        accepted = 0
        for t in finite_population():
            if t[0] == tuple(range(12)) or len(orbits(t)) != 1:
                continue
            accepted += 1
            for h in [t, *(hurwitz(t, m) for m in range(6))]:
                profile = analyze(h)
                self.assertEqual(profile['source_genus'], 2)
                self.assertEqual(profile['quotient_genus'], 1)
                self.assertTrue(profile['smooth'])
                self.assertEqual(profile['defect'], 0)
        self.assertGreater(accepted, 0)

    def test_relabeling_and_rejection(self):
        t = next(t for t in finite_population() if t[0] != tuple(range(12)) and len(orbits(t)) == 1)
        h = tuple(reversed(range(12)))
        relabeled = tuple(compose(compose(h, p), inverse(h)) for p in t)
        self.assertEqual(analyze(relabeled)['geometry'], 'isotrivial_elliptic')
        with self.assertRaises(ValueError):
            analyze((tuple(range(12)), *t[1:]))


if __name__ == '__main__':
    unittest.main()
