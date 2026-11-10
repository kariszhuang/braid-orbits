try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import combinations
import random
import unittest

import genus_three_map_core as maps
import genus_three_map_markings as marks


class MarkingTests(unittest.TestCase):
    def test_exact_clique_population_against_literal_all_supports(self):
        n=9;b=tuple((x+3)%n for x in range(n));e=tuple((x+6)%n for x in range(n));row=tuple(range(n))
        for k in(3,4):
            actual=set(marks.cliques(b,e,row,k));expected=set()
            for tail in combinations(range(1,n),k-1):
                support=(0,)+tail
                if not all(all(len(set(c)&{(x+t)%n for x in support})<=1 for c in maps.cycles(p))for p in(b,e)for t in range(n)):continue
                canonical=min(tuple(sorted((x-s)%n for x in support))for s in support)
                if support==canonical:expected.add(support)
            self.assertEqual(actual,expected)

    def test_middle_colour_identity_against_literal_both_square_moves(self):
        rng=random.Random(202610250247);count=0
        for n in(14,17,19,23,31):
            for _ in range(25):
                b,e=maps.grow(n,rng)
                for k in(3,4):
                    t=maps.quadruple(b,e,k,rng,rainbow=False)
                    if t is None:continue
                    positive,negative=marks.middle_square_failures(t)
                    self.assertEqual(positive,not maps.smooth(maps.braid(maps.braid(t,2),2)))
                    self.assertEqual(negative,not maps.smooth(maps.braid(maps.braid(t,3),3)))
                    count+=1
        self.assertGreater(count,100)


if __name__=='__main__':unittest.main()
