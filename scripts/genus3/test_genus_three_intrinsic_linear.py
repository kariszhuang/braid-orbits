"""Controls for literal inverse charts and oriented F3 transition matrices."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import unittest

from genus_three_intrinsic import uniform3_quotient,cycles
from genus_three_intrinsic_linear import LinearFrame,PROJECTIVE,PUSH_WORDS,LOCAL_WORDS,apply_word,project,matvec,bad_forms,initial_smooth,greatest_good,triangular_quotient,euclidean_quotient,minimal_block_genus,double_isogeny_quotients
from verify_genus_three_intrinsic import support_orbit,equivalent,connected
from verify_smooth_geometry import orbits


class LinearTests(unittest.TestCase):
    def setUp(self):self.frame=LinearFrame(uniform3_quotient(6),(1,3,4,0,5,2));self.mask=sum(1<<j for j in(0,1,3,5))

    def test_all_160_connected_lifts_and_four_zero_frames(self):
        for signs in itertools.product(range(2),repeat=2):
            for v in PROJECTIVE:
                t=self.frame.lift(self.mask,signs,v)
                self.assertEqual(len(orbits(t)),1)
                self.assertEqual(self.frame.extract(t),(self.mask,signs,v))
            zero=self.frame.lift(self.mask,signs,(0,0,0,0))
            self.assertGreater(len(orbits(zero)),1)
            self.assertEqual(self.frame.extract(zero),(self.mask,signs,(0,0,0,0)))

    def test_every_projective_point_matches_literal_point_push_matrix(self):
        for signs in itertools.product(range(2),repeat=2):
            for word in PUSH_WORDS:
                mask,target,matrix=self.frame.matrix(self.mask,signs,word)
                for v in PROJECTIVE:
                    actual=self.frame.extract(apply_word(self.frame.lift(self.mask,signs,v),word))
                    self.assertEqual(actual[:2],(mask,target))
                    self.assertEqual(actual[2],matvec(matrix,v))
                    self.assertIsNotNone(project(actual[2]))

    def test_local_hyperplanes_match_every_literal_index_test(self):
        for signs in itertools.product(range(2),repeat=2):
            for word in LOCAL_WORDS:
                forms=bad_forms(self.frame,self.mask,signs,word)
                for v in PROJECTIVE:
                    predicted=all(sum(a*b for a,b in zip(form,v))%3!=0 for form in forms)
                    self.assertEqual(predicted,initial_smooth(apply_word(self.frame.lift(self.mask,signs,v),word)))

    def test_reverse_bad_paths_are_literal_and_old_positive_family_survives(self):
        mask=next(m for m in range(64)if m.bit_count()==4 and support_orbit(self.frame.frame.beta,self.frame.frame.gamma,m)[1])
        orbit,eligible=support_orbit(self.frame.frame.beta,self.frame.frame.gamma,mask);self.assertTrue(eligible)
        result=greatest_good(self.frame,orbit);self.assertFalse(result['survivors'])
        for i,state in enumerate(result['states']):
            t=self.frame.lift(*state)
            if initial_smooth(t):self.assertFalse(initial_smooth(apply_word(t,result['witness'](i))))
        beta,gamma=triangular_quotient(2,0,2);f=LinearFrame(beta,gamma);mask=sum(1<<(3*j)for j in range(4))
        orbit,eligible=support_orbit(beta,gamma,mask);self.assertTrue(eligible);self.assertEqual(len(orbit),3)
        positive=greatest_good(f,orbit);self.assertEqual(len(positive['survivors']),480)
        self.assertEqual(minimal_block_genus(beta,gamma,mask)[0],1)
        self.assertEqual(minimal_block_genus(self.frame.frame.beta,self.frame.frame.gamma,self.mask)[0],0)

    def test_euclidean_lattice_frames_match_small_exact_classification(self):
        for order in(3,4,6):
            for m in range(1,8):
                if order*m>21:continue
                distinct=[]
                for a in range(1,m+1):
                    if m%a:continue
                    for b in range(a):
                        beta,gamma=euclidean_quotient(order,a,b,m//a)
                        self.assertTrue(connected(beta,gamma))
                        f=LinearFrame(beta,gamma)
                        expected=[3,3,3]if order==3 else([4,4,2]if order==4 else[6,3,2])
                        for p,length in zip((beta,gamma,f.frame.delta),expected):self.assertEqual(set(map(len,cycles(p))),{length})
                        if not any(equivalent((beta,gamma),old)for old in distinct):distinct.append((beta,gamma))
                if order==3 and m in(2,3):self.assertEqual(len(distinct),1 if m==2 else 2)

    def test_all_three_double_isogenies_recover_small_exact_families(self):
        for order,expected in((3,1),(4,2)):
            beta,gamma=euclidean_quotient(order,1,0,1);lifts=double_isogeny_quotients(beta,gamma);self.assertEqual(len(lifts),3);distinct=[]
            for b,c,h in lifts:
                self.assertTrue(connected(b,c));LinearFrame(b,c)
                if not any(equivalent((b,c),old)for old in distinct):distinct.append((b,c))
            self.assertEqual(len(distinct),expected)


if __name__=='__main__':unittest.main()
