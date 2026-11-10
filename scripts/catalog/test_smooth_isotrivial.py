try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import copy
import json
from pathlib import Path
import unittest
from classify_smooth_isotrivial import classify
from clean_smooth_catalog import clean_component, public_component
from test_clean_smooth_catalog import CleanCatalogTests
from verify_smooth_geometry import parse_tuple, hurwitz, cycles


class ClassificationTests(unittest.TestCase):
    def fixture(self):
        return clean_component(CleanCatalogTests().fixture())

    def test_elliptic_control(self):
        result=classify(self.fixture())
        self.assertIs(result['isotrivial'],True)
        self.assertNotIn('theorem_review',result['isotriviality'])
        w=result['isotriviality']['witnesses'][0]
        self.assertEqual((w['claim'],w['special_entry'],w['quotient_genus']),(6,1,1))

    def test_missing_criterion_is_unknown_not_false(self):
        e=self.fixture();e['block_quotients']=[]
        result=classify(e)
        self.assertIsNone(result['isotrivial'])
        self.assertEqual(result['isotriviality']['witnesses'],[])
        self.assertEqual(result['isotriviality']['failed_block_systems'],[])

    def test_no_trust_in_supplied_boolean(self):
        e=self.fixture();e['isotrivial']=False
        self.assertIs(public_component(e)['isotrivial'],True)

    def test_full_published_population_and_braid_invariance(self):
        path=Path(__file__).resolve().parents[2]/'exports/versioned-20261014/release-v1/all-smooth.json'
        if not path.exists():self.skipTest('local frozen catalog not present')
        source=json.loads(path.read_text());counts={True:0,False:0,None:0}
        def notation(p):
            return ''.join('('+' '.join(str(x+1) for x in c)+')' for c in cycles(p) if len(c)>1) or '()'
        for e in source['components']:
            result=classify(e);counts[result['isotrivial']]+=1
            t=parse_tuple('['+','.join(e['tuple'])+']',e['degree'])
            if result['isotrivial'] is None:
                systems=result['isotriviality']['failed_block_systems']
                self.assertEqual(len(systems),len(e['block_quotients']))
                for b,check in zip(e['block_quotients'],systems):
                    q=parse_tuple('['+','.join(b['tuple'])+']',b['block_count'])
                    expected=[i+1 for i,p in enumerate(q) if p==tuple(range(b['block_count']))]
                    self.assertEqual(check['a_passing_entries'],expected)
                    self.assertEqual([c['special_entry'] for c in check['candidate_failures']],expected)
                    for candidate in check['candidate_failures']:
                        self.assertTrue(candidate['failed_conditions'])
                        for failure in candidate['ordinary_cycle_count_failures']:
                            j=failure['entry']-1
                            self.assertEqual(failure['actual_cycles'],len(cycles(t[j])))
                            self.assertEqual(failure['required_cycles'],b['block_size']*len(cycles(q[j])))
                            self.assertGreater(failure['ramification_excess'],0)
            else:
                self.assertNotIn('failed_block_systems',result['isotriviality'])
            if e['id']=='1ed71cb93a88a54f0d395b71c0a35538a1d66e848642e0b436835e9f7aea3a02':
                fine,coarse=result['isotriviality']['failed_block_systems']
                self.assertEqual(fine['a_passing_entries'],[])
                self.assertEqual(fine['quotient_cycle_counts'],[7,3,3,3])
                c=coarse['candidate_failures'][0]
                self.assertEqual(c['special_entry'],1)
                self.assertEqual(c['active_blocks'],[1])
                self.assertEqual(c['failed_conditions'],['c'])
                self.assertEqual([(f['entry'],f['actual_cycles'],f['required_cycles'])
                    for f in c['ordinary_cycle_count_failures']],[(3,5,6),(4,5,6)])
            for move in range(6):
                new=copy.deepcopy(e);new['tuple']=[notation(p) for p in hurwitz(t,move)]
                for b in new['block_quotients']:
                    q=parse_tuple('['+','.join(b['tuple'])+']',b['block_count'])
                    b['tuple']=[notation(p) for p in hurwitz(q,move)]
                self.assertIs(classify(new)['isotrivial'],result['isotrivial'])
        self.assertEqual(counts,{True:5126,False:2895,None:74})


if __name__=='__main__':unittest.main()
