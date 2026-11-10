try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import unittest
from run_smooth_frontier import split_geodesics
from run_smooth_small_census import partition
from analyze_fixed_braid_constraints import product
from odd_order_spin import inverse, cycle_lengths


class FrontierTests(unittest.TestCase):
    def test_split_enumerator_equals_literal_all_permutations(self):
        for n in range(2,7):
            perms = list(itertools.permutations(range(n)))
            for lens in partition(n):
                p = list(range(n)); off = 0
                for length in lens:
                    for j in range(length): p[off+j] = off+(j+1)%length
                    off += length
                p = tuple(p); ip = n-len(lens)
                literal = {a for a in perms if 2*n-len(cycle_lengths(a))-len(cycle_lengths(product(inverse(a),p))) == ip}
                split = list(split_geodesics(p))
                self.assertEqual(len(split), len(set(split)))
                self.assertEqual(set(split),literal)


if __name__ == '__main__': unittest.main()
