"""Independent state comparison and square-root frame action regressions."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import unittest

from analyze_fixed_braid_constraints import conjugate, power, product
from odd_order_spin import inverse, hurwitz
from run_s18_linked_roots import frame_group, state_key


class LinkedRootChecks(unittest.TestCase):
    def test_state_keys_ignore_relabeling_but_not_entry_order(self):
        a=(1,2,3,4,0)
        b=(1,0,2,3,4)
        c=(0,2,1,3,4)
        t=(a,b,c,inverse(product(a,b,c)))
        expected=state_key(t)
        for g in itertools.permutations(range(5)):
            self.assertEqual(state_key(tuple(conjugate(p,g) for p in t)),expected)
        self.assertNotEqual(state_key((b,a,c,inverse(product(b,a,c)))),expected)

    def test_regular_frame_centralizer(self):
        h=tuple(b*9+((i+1)%3)*3+j for b in range(2) for i in range(3) for j in range(3))
        k=tuple(b*9+i*3+(j+1)%3 for b in range(2) for i in range(3) for j in range(3))
        group=set(frame_group(h,k))
        self.assertEqual(len(group),162)
        self.assertEqual(power(h,3),tuple(range(18)))
        self.assertTrue(all(conjugate(h,g)==h and conjugate(k,g)==k for g in group))
        self.assertTrue(all(inverse(g) in group for g in group))
        self.assertTrue(all(product(a,b) in group for a in group for b in group))

    def test_hurwitz_word_and_inverse(self):
        a=(1,2,3,4,0); b=(1,0,2,3,4); c=(0,2,1,3,4)
        t=(a,b,c,inverse(product(a,b,c)))
        for mv in range(3):
            self.assertEqual(hurwitz(hurwitz(t,mv),mv+3),t)
            self.assertEqual(product(*hurwitz(t,mv)),tuple(range(5)))


if __name__=='__main__':
    unittest.main()
