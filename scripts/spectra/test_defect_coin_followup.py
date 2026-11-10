"""Independent tiny controls for coin recurrences and literal seed handling."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import math
import unittest
from fractions import Fraction

from defect_coin_followup import compare, decode, cycle_type, passport
from test_defect_coin_model import coins


class CoinControls(unittest.TestCase):
    def test_binary_outcomes(self):
        for g in range(9):
            hist=[Fraction(0)]*(g+1)
            for bits in itertools.product((0,1),repeat=g):
                p=Fraction(1)
                for j,b in enumerate(bits,1):p*=Fraction(1 if b else j,j+1)
                hist[sum(bits)]+=p
            self.assertEqual(hist,coins(g))
            if g:self.assertEqual(sum(hist[::2]),Fraction(1,2))

    def test_random_permutation_identity(self):
        for n in range(1,8):
            hist=[0]*n
            for p in itertools.permutations(range(n)):hist[len(cycle_type(p))-1]+=1
            self.assertEqual([Fraction(x,math.factorial(n)) for x in hist],coins(n-1))

    def test_full_cycles_parity(self):
        for n in range(3,6):
            cycles=[p for p in itertools.permutations(range(n)) if cycle_type(p)==(n,)]
            for a in cycles:
                for b in cycles:
                    product=tuple(a[b[i]] for i in range(n))
                    self.assertEqual((len(cycle_type(product))-1)%2,(n-1)%2)

    def test_fixed_points_and_padding(self):
        t=decode('[(1 2), (2 3), (1 3 2), ()]',4)
        self.assertEqual(cycle_type(t[0]),(2,1,1))
        self.assertEqual(passport(t),passport(t[::-1]))
        # Reference beyond the observed genus must NOT disappear in zip().
        self.assertAlmostEqual(compare({0:1},1,4)['tv_degree_cycles'],.75)
        self.assertAlmostEqual(compare({0:1},1,5,'A')['tv_coin'],.5)


if __name__=='__main__':unittest.main()
