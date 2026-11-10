try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
from clean_smooth_catalog import clean_component


class CleanCatalogTests(unittest.TestCase):
    def fixture(self):
        return dict(degree=9,genus=2,quotient_state_count=12,all_states_smooth=True,
            group=dict(order_decimal='324'),mechanism_at_representative='single-local/unramified',
            defect_distribution=[dict(defect=0,states=12)],
            representative=dict(permutations=[dict(cycle_notation=s) for s in ['(5 6 7)','(1 2 5)(3 7 9)(4 6 8)','(1 3 7)(2 6 9)(4 5 8)','(1 4 5)(2 7 8)(3 6 9)']],
                passport_ordered=[[3,1,1,1,1,1,1],[3,3,3],[3,3,3],[3,3,3]],indices=[2,6,6,6],
                adjacent_pairs=[dict(orbit_sizes=[9],product_index=8)]*3),
            block_lattice=dict(complete=True,primitive=False,primitive_quotient_indices=[0],proper_systems=[dict(
                blocks=3,block_size=3,partition=[[1,8,9],[2,3,4],[5,6,7]],genus=1,indices=[0,2,2,2],ramification_excess=[2,0,0,0],
                quotient_permutations=[dict(cycle_notation='()',images_1_based=[1,2,3])]+[dict(cycle_notation='(1 2 3)',images_1_based=[2,3,1])]*3)]))

    def test_compact_results_keep_primitive_and_numeric_filter_fields(self):
        e=clean_component(self.fixture());self.assertFalse(e['primitive']);self.assertEqual(e['group_order'],324)
        self.assertEqual(e['quotient_genera'],[1]);self.assertEqual(e['block_sizes'],[3])
        self.assertEqual(e['block_quotients'][0]['group_order'],3)
        self.assertEqual(e['block_quotients'][0]['kernel_order'],108)
        self.assertEqual(e['contained_in'],'A9')
        self.assertEqual(e['representative_g1_g2'],dict(permutation='(1 2 6 8 4 7 9 3 5)',index=8,cycle_count=1))
        for field in ['structure_type','generates_full_A9','pure_braid_orbits','evidence','geometry']:self.assertNotIn(field,e)

    def test_wrong_scope_and_incomplete_primitive_analysis_rejected(self):
        for key,value in [('genus',0),('all_states_smooth',False)]:
            e=self.fixture();e[key]=value
            with self.assertRaises(ValueError):clean_component(e)
        e=self.fixture();e['block_lattice']['complete']=False
        with self.assertRaises(ValueError):clean_component(e)

    def test_unsafe_numeric_order_is_not_silently_rounded(self):
        e=self.fixture();e['group']['order_decimal']=str(2**53)
        with self.assertRaises(ValueError):clean_component(e)

    def test_primitive_candidate_needs_no_invented_block_quotient(self):
        e=self.fixture();e['mechanism_at_representative']='unclassified'
        e['block_lattice'].update(primitive=True,proper_systems=[],primitive_quotient_indices=[])
        result=clean_component(e)
        self.assertTrue(result['primitive']);self.assertEqual(result['block_quotients'],[])
        self.assertNotIn('structure_type',result)


if __name__=='__main__':unittest.main()
