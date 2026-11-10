"""Mixed-frame literal local return and complete native fibre controls."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
from unittest.mock import patch

from genus_three_intrinsic_mixed import MixedFrame,fixed_vertex_pairs


class MixedTests(unittest.TestCase):
    def test_every_q4_fixed_vertex_pair_has160_cubic_covers(self):
        beta=gamma=(1,2,3,0)
        for fixed in fixed_vertex_pairs(beta,gamma):
            frame=MixedFrame(beta,gamma,fixed,m=3);rows=list(frame.covers(3));self.assertEqual(len(rows),160)
            self.assertEqual(frame.population['canonical_full'],160)

    def test_all_four_punctures_including_face_have2400_quartic_classes(self):
        frame=MixedFrame((1,2,3,0),(1,2,3,0),(('b',0),('d',1)),m=4)
        self.assertEqual(len(list(frame.covers(3))),2400)

    def test_three_moving_one_fixed_face_complete_quartic_fibre(self):
        from genus_three_intrinsic import cycles
        frame=MixedFrame((1,2,3,0),(1,2,3,0),(('d',0),),m=4)
        rows=list(frame.covers(7));self.assertEqual(len(rows),2400)
        self.assertEqual(frame.population['canonical_full'],2400)
        for _,t,_ in rows:
            self.assertEqual([16-len(cycles(p))for p in t],[3,12,12,9])

    def test_three_moving_marked_handle_and_old_stabilizer_enumerations_agree(self):
        from genus_three_intrinsic_native import ConjugacyFrame
        frame=MixedFrame((1,2,3,0),(1,2,3,0),(('d',0),),m=4)
        actual={frame.rooted_key(meta['native_data'])for _,_,meta in frame.covers(7)}
        reference={frame.rooted_key(meta['native_data'])for _,_,meta in ConjugacyFrame.covers(frame,7)}
        self.assertEqual(actual,reference);self.assertEqual(len(actual),2400)

    def test_explicit_deadline_is_still_enforced_and_default_has_no_expired_constant(self):
        frame=MixedFrame((1,2,3,0),(1,2,3,0),(('d',0),),m=3,deadline_epoch=100)
        with patch('genus_three_intrinsic_native.time.time',return_value=101):
            with self.assertRaises(TimeoutError):frame.check_deadline()
        frame.deadline_epoch=None
        with patch('genus_three_intrinsic_native.time.time',return_value=10**12):frame.check_deadline()


if __name__=='__main__':unittest.main()
