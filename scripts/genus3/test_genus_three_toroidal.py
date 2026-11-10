"""Genus-preserving operations and exact small-map coverage controls."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import itertools
from pathlib import Path
import random
import unittest

from genus_three_intrinsic import cycles,compose
from genus_three_toroidal import leaf,split_face,validate,grow,marked_key


class ToroidalTests(unittest.TestCase):
    def test_every_small_insertion_is_elliptic(self):
        beta=gamma=(1,2,0)
        for colour in('black','white'):
            for j in range(3):validate(*leaf(beta,gamma,colour,j))
        for u,k in itertools.product(range(3),repeat=2):validate(*split_face(beta,gamma,u,k))

    def test_seeded_growth_is_diverse_and_elliptic(self):
        rng=random.Random(20261023);profiles=set()
        for q in range(5,22):
            for _ in range(5):
                beta,gamma=grow(q,rng,need_two_unramified=True);validate(beta,gamma)
                profiles.add(tuple(tuple(sorted(map(len,cycles(p))))for p in(beta,gamma,compose(beta,gamma))))
        self.assertGreater(len(profiles),60)

    def test_marked_key_is_relabeling_invariant(self):
        beta,gamma=grow(10,random.Random(27));mask=3;fixed=(('b',0),('d',0));h=(9,2,5,0,8,1,7,3,6,4)
        def change(p):
            out=[0]*10
            for x in range(10):out[h[x]]=h[p[x]]
            return tuple(out)
        bb,gg=change(beta),change(gamma);old={'b':cycles(beta),'g':cycles(gamma),'d':cycles(compose(beta,gamma))};new={'b':cycles(bb),'g':cycles(gg),'d':cycles(compose(bb,gg))}
        ff=tuple((kind,next(j for j,row in enumerate(new[kind])if set(row)=={h[x]for x in old[kind][i]}))for kind,i in fixed)
        self.assertEqual(marked_key(beta,gamma,mask,fixed),marked_key(bb,gg,sum(1<<h[x]for x in(0,1)),ff))


if __name__=='__main__':unittest.main()
