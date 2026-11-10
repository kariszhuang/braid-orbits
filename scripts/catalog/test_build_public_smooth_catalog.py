try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
import json
from pathlib import Path
import tempfile
from build_public_smooth_catalog import build, delta, release
from clean_smooth_catalog import public_component, clean_component
from test_clean_smooth_catalog import CleanCatalogTests


class PublicationTests(unittest.TestCase):
    def fixture(self):
        e=clean_component(CleanCatalogTests().fixture());e['id']='a'*64
        return e

    def source(self):return dict(component_count=1,components=[self.fixture()])

    def test_public_envelope_and_scope(self):
        result=build(dict(component_count=1,components=[self.fixture()]))
        self.assertEqual(result['component_count'],1)
        self.assertEqual(result['total_state_count'],12)
        self.assertIn('+00:00',result['last_updated'])
        self.assertTrue(result['components'][0]['all_states_smooth'])
        self.assertNotIn('structure_type',result['components'][0])
        self.assertEqual(result['dataset_version'],1)
        self.assertEqual(result['components'][0]['revision'],1)
        self.assertEqual(result['components'][0]['last_updated'],result['last_updated'])

    def test_unchanged_publication_is_noop(self):
        first=build(self.source(),now='2026-10-14T17:00:00Z')
        again=build(self.source(),first,now='2026-10-14T18:00:00Z')
        self.assertEqual(first,again)
        self.assertEqual(delta(again,first)['upserted'],[])

    def test_changed_representative_preserves_id_but_increments_revision(self):
        first=build(self.source(),now='2026-10-14T17:00:00Z')
        source=self.source();source['components'][0]['tuple'][0]='(6 7 5)'
        next=build(source,first,now='2026-10-14T18:00:00Z')
        self.assertEqual(next['dataset_version'],2)
        self.assertEqual(next['components'][0]['id'],first['components'][0]['id'])
        self.assertEqual(next['components'][0]['revision'],2)
        self.assertEqual(len(delta(next,first)['upserted']),1)

    def test_delete_restore_and_delete_have_monotone_revisions(self):
        first=build(self.source(),now='2026-10-14T17:00:00Z')
        removed=build(dict(component_count=0,components=[]),first,now='2026-10-14T18:00:00Z')
        self.assertEqual(delta(removed,first)['deleted'][0]['revision'],2)
        restored=build(self.source(),removed,now='2026-10-14T19:00:00Z')
        self.assertEqual(restored['components'][0]['revision'],3)
        self.assertEqual(restored['deleted_components'],[])

    def test_invalid_ids_and_backdated_changes_rejected(self):
        e=self.source();e['components'][0]['id']='not-a-component-id'
        with self.assertRaises(ValueError):build(e)
        first=build(self.source(),now='2026-10-14T17:00:00Z')
        with self.assertRaises(ValueError):build(dict(component_count=0,components=[]),first,now='2026-10-14T16:00:00Z')

    def test_immutable_history_retained_and_delta_reconstructs_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a';b=Path(tmp)/'b'
            initial=release(self.source(),a,now='2026-10-14T17:00:00Z')
            source=self.source();source['components'][0]['tuple'][0]='(6 7 5)'
            # Explicit older baseline makes test independent of subsecond timing.
            old=json.loads((a/'all-smooth.json').read_text())
            next=build(source,old,now='2099-01-01T00:00:00Z')
            changes=delta(next,old);state={e['id']:e for e in old['components']}
            state.update({e['id']:e for e in changes['upserted']})
            for e in changes['deleted']:state.pop(e['id'])
            self.assertEqual(list(state.values()),next['components'])
            second=release(source,b,a,now='2099-01-01T00:00:00Z')
            self.assertEqual(second['latest_version'],2)
            self.assertEqual((a/'versions/1.json').read_bytes(),(b/'versions/1.json').read_bytes())
            self.assertEqual(json.loads((b/'changes/2.json').read_text()),changes)
            self.assertIsNone(release(source,Path(tmp)/'noop',b))
            with (a/'changes/1.json').open('ab') as stream:stream.write(b' ')
            with self.assertRaises(ValueError):release(source,Path(tmp)/'bad',a)

    def test_multiple_product_cycles(self):
        e=dict(degree=24,genus=4,states=576,primitive=False,group_order=98304,
            tuple=['(5 8)(21 23)',
                '(1 21 17)(2 22 18)(3 23 19)(4 24 20)(5 13 9)(6 14 10)(7 15 11)(8 16 12)',
                '(1 8 16 4 7 14)(2 6 15 3 5 13)(9 23 17)(10 22 20)(11 24 18)(12 21 19)',
                '(1 23 14 2 24 16)(3 21 15 4 22 13)(5 19 9)(6 20 11)(7 18 10)(8 17 12)'],
            cycle_types=[[2,2]+[1]*20,[3]*8,[6,6,3,3,3,3],[6,6,3,3,3,3]],indices=[2,16,18,18],
            defect_distribution=[dict(defect=0,states=576)],
            representative_pair_orbit_sizes=dict(g1_g2=[6,6,3,3,3,3],g2_g3=[24],g3_g4=[24]),
            quotient_genera=[],block_sizes=[],block_quotients=[])
        result=public_component(e)
        self.assertEqual(result['representative_g1_g2']['cycle_count'],6)
        self.assertEqual(result['representative_g1_g2']['index'],18)
        self.assertEqual(result['contained_in'],'A24')

    def test_corrupt_input_fails_closed(self):
        e=self.fixture();e['indices']=[0,6,6,6]
        with self.assertRaises(ValueError):public_component(e)
        with self.assertRaises(ValueError):build(dict(component_count=2,components=[self.fixture()]))
        e=self.fixture();e['block_quotients'][0]['tuple'][1]='(1 3 2)'
        with self.assertRaises(ValueError):public_component(e)
        e=self.fixture();e['genus']=0
        with self.assertRaises(ValueError):public_component(e)


if __name__=='__main__':unittest.main()
