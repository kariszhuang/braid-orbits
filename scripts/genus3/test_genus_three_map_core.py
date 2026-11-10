try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import random
import unittest

import genus_three_map_core as maps


class CoreTests(unittest.TestCase):
    def test_all_seed_pair_and_leaf_face_invariants(self):
        for n in(7,8,9,12):
            b,e=maps.seed(n);self.assertEqual(len(maps.cycles(maps.compose(b,e))),1)
            for corner in range(n):
                for colour in range(2):
                    bb,ee=maps.leaf(b,e,colour,corner);maps.validate(bb,ee)
                    self.assertEqual(len(maps.cycles(maps.compose(bb,ee))),1)
            for u in range(n):
                for k in range(n):
                    bb,ee=maps.split_face(b,e,u,k);maps.validate(bb,ee)
                    self.assertEqual(len(maps.cycles(maps.compose(bb,ee))),2)

    def test_genus_preserving_mixes_are_literal_and_can_break_commutation(self):
        b,e=maps.seed(8);accepted=0;noncommuting=0
        for x in range(8):
            for y in range(x+1,8):
                result=maps.mix(b,e,x,y)
                if result is None:continue
                bb,ee=result;maps.validate(bb,ee);accepted+=1
                self.assertEqual(maps.compose(bb,ee),maps.compose(b,e))
                noncommuting+=maps.compose(bb,ee)!=maps.compose(ee,bb)
        self.assertGreater(accepted,0);self.assertGreater(noncommuting,0)

    def test_grown_quadruple_product_genus_and_exact_state_key(self):
        rng=random.Random(202610240245);faces=set();starts=0
        for n in(14,17,19,23,31):
            for _ in range(6):
                b,e=maps.grow(n,rng);faces.add(len(maps.cycles(maps.compose(b,e))))
                for k in(3,4):
                    t=maps.quadruple(b,e,k,rng,rainbow=False)
                    if t is None:continue
                    starts+=1;self.assertTrue(maps.smooth(t));self.assertEqual(sum(map(maps.index,t)),2*n+4)
                    lab=list(range(n));rng.shuffle(lab);inv=maps.inverse(lab)
                    relabelled=tuple(maps.compose(maps.compose(lab,p),inv)for p in t)
                    self.assertEqual(maps.canonical(t),maps.canonical(relabelled))
        self.assertGreater(starts,20);self.assertGreater(len(faces),1)

    def test_cyclic_rainbow_support_and_no_failed_try_as_exclusion(self):
        rng=random.Random(202610240245);n=15
        b=tuple((x+3)%n for x in range(n));e=tuple((x+6)%n for x in range(n))
        row=tuple(range(n));support=maps.initial_rainbow_support(b,e,row,3,rng)
        self.assertIsNotNone(support)
        for t in range(n):
            shifted={(x+t)%n for x in support}
            for p in(b,e):
                self.assertTrue(all(len(shifted&set(c))<=1 for c in maps.cycles(p)))
        full=tuple((x+1)%n for x in range(n))
        self.assertIsNone(maps.initial_rainbow_support(full,full,row,3,rng))


if __name__=='__main__':unittest.main()
