"""Native chart cross-checks against independently parameterized small fibres."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import unittest

from genus_three_intrinsic_native import NativeFrame,ConjugacyFrame,MarkedCycleFrame
from genus_three_intrinsic_four import QuarticFrame,FOUR,MUL4,INV4
from genus_three_intrinsic_linear import SIGN_STATES,PROJECTIVE


class NativeTests(unittest.TestCase):
    def test_marked_cycle_B_classes_are_exact_centralizer_orbits(self):
        frame=MarkedCycleFrame((1,2,3,0),(1,2,3,0),m=5,lazy=True)
        actual={min(c[b]for c in frame.conjugate.values())for b in frame.second_handles()}
        expected={min(c[b]for c in frame.conjugate.values())for b in range(frame.size)}
        self.assertEqual(actual,expected);self.assertEqual(len(actual),len(frame.second_handles()))

    def test_marked_cycle_fibre_matches_every_direct_S5_class(self):
        beta=gamma=(1,2,3,0);direct=NativeFrame(beta,gamma,m=5);fast=MarkedCycleFrame(beta,gamma,m=5,lazy=True)
        expected={direct.rooted_key(meta['native_data'])for _,_,meta in direct.covers(15)}
        actual={fast.rooted_key(meta['native_data'])for _,_,meta in fast.covers(15)}
        self.assertEqual(len(expected),18304);self.assertEqual(actual,expected)

    def test_conjugacy_constructor_matches_entire_direct_S5_fibre(self):
        beta=gamma=(1,2,3,0);direct=NativeFrame(beta,gamma,m=5);fast=ConjugacyFrame(beta,gamma,m=5)
        def key(frame,data):return min(tuple(c[p]for p in data)for c in frame.conjugate.values())
        expected={key(direct,meta['native_data'])for _,_,meta in direct.covers(15)}
        actual=set()
        for _,_,meta in fast.covers(15):
            self.assertTrue(fast.full_group_reference(meta['native_data']))
            actual.add(key(fast,meta['native_data']))
        self.assertEqual(len(expected),18304);self.assertEqual(actual,expected)
        self.assertEqual(direct.population['canonical_all'],fast.population['canonical_all'])

    def test_cubic_fibre_has160_classes(self):
        frame=NativeFrame((1,2,0,4,5,3),(1,3,4,0,5,2),m=3)
        self.assertEqual(len(list(frame.covers(29))),160)
        self.assertEqual(frame.population['centralizer_size'],2)

    def test_quartic_native_matches_every_affine_class(self):
        beta=(1,2,0,4,5,3);gamma=(1,3,4,0,5,2);mask=29;active=[j for j in range(6)if mask&(1<<j)]
        native=NativeFrame(beta,gamma,m=4)
        def canonical(row):
            return min(tuple(tuple(h[p[h.index(x)]]for x in range(4))for p in row)for h in native.perms)
        expected=set()
        affine=QuarticFrame(beta,gamma)
        for signs in SIGN_STATES:
            for vector in PROJECTIVE:
                codes,values,_=affine.extension_space(mask,signs,vector)
                for code in codes:
                    local,handles=values(code)
                    expected.add(canonical(tuple(FOUR[p]for p in(*handles,*(local[j]for j in active)))))
        actual={canonical(tuple(native.perms[p]for p in meta['native_data']))for _,_,meta in native.covers(mask)}
        self.assertEqual(len(actual),2400);self.assertEqual(actual,expected)


if __name__=='__main__':unittest.main()
