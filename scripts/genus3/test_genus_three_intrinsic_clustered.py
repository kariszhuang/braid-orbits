try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest

from genus_three_intrinsic import cycles
from genus_three_intrinsic_clustered import ClusteredFrame


class ClusteredTests(unittest.TestCase):
    def test_three_vertex_markings_same_complete_S5_fibre_count(self):
        beta=gamma=(1,2,0)
        for fixed in(('b',0),('g',0),('d',0)):
            frame=ClusteredFrame(beta,gamma,fixed,m=5);rows=list(frame.covers(3))
            self.assertEqual(len(rows),3600)
            for _,t,meta in rows:
                self.assertTrue(frame.full_group_reference(meta['native_data']))
                indices=[15-len(cycles(p))for p in t]
                self.assertEqual(sorted(indices),[2,10,10,12])

    def test_normal_transposition_closure_matches_explicit_group_for_all_small_data(self):
        from itertools import product
        frame=ClusteredFrame((1,2,0),(1,2,0),('d',0),m=3)
        for a,b,t,s in product(range(frame.size),range(frame.size),frame.trans,frame.threes):
            data=(a,b,frame.fixed,t,s)
            self.assertEqual(frame.full_group(data),frame.full_group_reference(data))


if __name__=='__main__':unittest.main()
