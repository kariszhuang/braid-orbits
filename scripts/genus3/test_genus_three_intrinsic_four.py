"""Exact quartic binary-extension controls; never filter cubic-bad frames."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest

from genus_three_intrinsic import uniform3_quotient
from genus_three_intrinsic_four import QuarticFrame,FOUR,MUL4,INV4
from verify_smooth_geometry import compose,orbits


class QuarticTests(unittest.TestCase):
    def test_affine_table_is_literal_S4(self):
        self.assertEqual(len(set(FOUR)),24)
        for i in range(24):
            self.assertEqual(FOUR[MUL4[i][INV4[i]]],(0,1,2,3))
            for j in range(24):self.assertEqual(FOUR[MUL4[i][j]],compose(FOUR[i],FOUR[j]))

    def test_all160_root_maps_and2400_full_S4_connected_covers(self):
        frame=QuarticFrame(uniform3_quotient(6),(1,3,4,0,5,2));mask=sum(1<<j for j in(0,1,3,5));count=0
        for i,t,meta in frame.covers(mask,check_all=True):
            count+=1;self.assertEqual(len(orbits(t)),1);self.assertEqual(compose(compose(t[0],t[1]),compose(t[2],t[3])),tuple(range(24)))
        self.assertEqual(count,2400)


if __name__=='__main__':unittest.main()
