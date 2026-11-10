try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from pathlib import Path
import sys
import unittest
from export_research_catalog import tuple_details,block_details,mechanism,identifier
from verify_smooth_geometry import parse_tuple, hurwitz
from run_s18_linked_roots import state_key


class ExportTests(unittest.TestCase):
    def test_fixed_points_and_defect_are_explicit(self):
        t=parse_tuple('[(1 2), (2 3), (1 2 3), (1 2 3)]',3);v=tuple_details(t)
        self.assertEqual(v['genus'],1);self.assertEqual(v['defect'],0)
        self.assertTrue(v['smooth_at_representative']);self.assertEqual(v['indices'],[1,1,2,2])
        self.assertEqual(v['permutations'][0]['fixed_points_1_based'],[3])
        self.assertEqual(v['permutations'][0]['cycle_lengths'],[2,1])

    def test_invalid_product_is_not_exported(self):
        with self.assertRaises(ValueError):tuple_details(parse_tuple('[(1 2), (2 3), (1 3 2), (1 3 2)]',3))

    def test_known_nine_point_block_quotient(self):
        t=parse_tuple('[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]',9)
        b=block_details(t);self.assertTrue(b['complete']);self.assertFalse(b['primitive'])
        self.assertEqual(mechanism(b),'single-local/unramified')
        self.assertEqual(b['proper_systems'][0]['genus'],1)
        self.assertEqual(b['proper_systems'][0]['quotient_permutations'][0]['images_1_based'],[1,2,3])

    def test_identity_hash_is_degree_sensitive_and_deterministic(self):
        key=((0,1),(1,0),(0,1))
        self.assertEqual(identifier(2,key),identifier(2,key));self.assertNotEqual(identifier(2,key),identifier(3,key))

    def test_complete_component_id_survives_braid_move_and_relabeling(self):
        t=parse_tuple('[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]',9)
        def full_id(seed):
            queue=[seed];seen={state_key(seed)}
            for row in queue:
                for m in range(6):
                    next=hurwitz(row,m);key=state_key(next)
                    if key not in seen:seen.add(key);queue.append(next)
            self.assertEqual(len(seen),12)
            return identifier(9,min(seen))
        identity=full_id(t)
        labels=(1,0,2,3,4,5,6,7,8)
        renamed=tuple(tuple(labels[p[labels[x]]] for x in range(9)) for p in t)
        self.assertEqual(identity,full_id(renamed));self.assertEqual(identity,full_id(hurwitz(t,1)))
        self.assertEqual(identity,'70ca636db54773894178bd67b570d7336ce0b71bfc2aa49f0ea462582e17cd61')


if __name__=='__main__':unittest.main()
