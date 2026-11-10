"""Literal small-state tests for the index-family fixed-point obstruction."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import permutations
import unittest
from verify_genus_three_finite import mul,inv,profile,braid
from verify_genus_three_twists import fixed_point_word,smooth


class TwistTests(unittest.TestCase):
    def test_triangle_boundary_and_degree14_overlap_word(self):
        n=14;budgets=[]
        for a in range(1,n):
            for b in range(a,n):
                for c in range(b,n):
                    e=2*n+4-a-b-c
                    if c<=e<n and a+c<n and b+c<n:budgets.append((a,b,c,e))
        self.assertEqual(budgets,[(6,6,7,13)])
        c=tuple((-i-1)%n for i in range(n));a=tuple((-i)%n for i in range(n))
        b=list((i+7)%n for i in range(n));b[0]=0;b[7]=7;b=tuple(b)
        p=mul(c,a);d=inv(mul(p,b));t=(c,a,b,d)
        self.assertEqual([n-len(profile(q))for q in t],[7,6,6,13]);self.assertTrue(smooth(t))
        u=braid(braid(t,3),4);self.assertEqual(u[0],a);self.assertEqual(u[1],b);self.assertFalse(smooth(u))

    def test_degree15_support_bound_profile_lists(self):
        def partitions(total,k,minimum=1):
            if k==0:
                if total==0:yield ()
                return
            for x in range(minimum,total//k+1):
                for tail in partitions(total-x,k-1,x):yield (x,)+tail
        profiles=list(partitions(15,5));self.assertEqual(len(profiles),30)
        self.assertEqual([q for q in profiles if min(q)>=3],[(3,)*5])
        self.assertEqual([q for q in profiles if min(q)>=2 and max(q)<=7],
            [(2,2,2,2,7),(2,2,2,3,6),(2,2,2,4,5),(2,2,3,3,5),(2,2,3,4,4),(2,3,3,3,4),(3,3,3,3,3)])

    def test_all_small_half_rank_frames(self):
        n=5;p=tuple((x+1)%n for x in range(n));candidates=[a for a in permutations(range(n))if len(profile(a))==3]
        count=0;noninvolutions=0
        for a in candidates:
            b=mul(inv(a),p)
            if len(profile(b))!=3:continue
            for c in candidates:
                d=inv(mul(p,c))
                if profile(d)!=(n,):continue
                t=(a,b,c,d);word=fixed_point_word(t)
                self.assertIsNotNone(word);self.assertLessEqual(len(word),n);u=t
                for m in word:u=braid(u,m)
                self.assertFalse(smooth(u));count+=1
                noninvolutions+=any(x==3 for x in profile(a))
        self.assertGreater(count,0);self.assertGreater(noninvolutions,0)

    def test_unsupported_is_undecided(self):
        identity=tuple(range(5))
        self.assertIsNone(fixed_point_word((identity,)*4))


if __name__=='__main__':unittest.main()
