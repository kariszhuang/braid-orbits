try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import unittest
from run_smooth_mission import REPORTED, values
from odd_order_spin import parse_seed, hurwitz, inverse
from analyze_fixed_braid_constraints import product, point_orbits
from run_s18_linked_roots import state_key


class SmoothChecks(unittest.TestCase):
    def test_reported_all_smooth(self):
        seed=parse_seed(REPORTED,9);states=[seed];seen={state_key(seed)}
        for t in states:
            for m in range(6):
                h=hurwitz(t,m);key=state_key(h)
                if key not in seen:seen.add(key);states.append(h)
        self.assertEqual(len(states),12)
        self.assertTrue(all(values(t)['smooth'] and values(t)['defect']==0 and values(t)['genus']==2 for t in states))

    def test_all_transitive_product_one_S3(self):
        perms=list(itertools.permutations(range(3)));tested=0
        for a,b,c in itertools.product(perms,repeat=3):
            t=(a,b,c,inverse(product(a,b,c)))
            if len(point_orbits(t))!=1:continue
            v=values(t);tested+=1
            if v['smooth']:
                self.assertEqual(v['defect'],0)
                self.assertTrue(v['left_orbits']==1 or v['right_orbits']==1)
        self.assertGreater(tested,0)

    def test_nontransitivity_matters(self):
        t=parse_seed('[(1 2), (1 2), (), ()]',3)
        self.assertTrue(values(t)['smooth']);self.assertEqual(values(t)['defect'],-1)

    def test_block_construction_algebra(self):
        for k in range(3,9):
            x=tuple((i+1)%k for i in range(k));y=(1,0,*range(2,k));z=inverse(product(y,x))
            tau=product(inverse(x),y,x,inverse(y));n=3*k
            a=tuple(tau)+tuple(range(k,n))
            b=tuple(((j+1)%3)*k+i for j in range(3) for i in range(k))
            c=tuple(((j+1)%3)*k+[x,y,z][j][i] for j in range(3) for i in range(k))
            t=(a,b,c,inverse(product(a,b,c)))
            self.assertEqual(product(*t),tuple(range(n)))
            self.assertEqual(len(point_orbits(t)),1)
            self.assertEqual(sorted(__import__('odd_order_spin').cycle_lengths(a)),[1]*(n-3)+[3])
            for p in t[1:]:self.assertEqual(__import__('odd_order_spin').cycle_lengths(p),[3]*k)
            for j in range(100):
                self.assertTrue(values(t)['smooth']);self.assertEqual(values(t)['defect'],0)
                self.assertEqual(values(t)['genus'],2);t=hurwitz(t,(j*7+j//3)%6)


if __name__=='__main__':unittest.main()
