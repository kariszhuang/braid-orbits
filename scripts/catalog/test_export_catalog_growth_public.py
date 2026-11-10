try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import copy
import unittest
from clean_smooth_catalog import public_component, clean_component
from test_clean_smooth_catalog import CleanCatalogTests

class GrowthPublicationTests(unittest.TestCase):
    def entry(self):
        r=clean_component(CleanCatalogTests().fixture())
        r.update(group_order=None,group_order_status='not_computed',block_quotients_complete=False)
        for b in r['block_quotients']:
            b.update(group_order=None,kernel_order=None,orders_status='not_computed')
        r['hodge_degree']=dict(normalization='12L=d-sum(rho(entries))+3*mean(rho(g1*g2))',
            twelve_L='0',L='0',pair_product_cycle_distribution=[dict(cycle_lengths=[9],states=12)])
        r['isotriviality']=dict(status='isotrivial',method='full-orbit-hodge-degree-v1',
            source=None,explanation='Exact Hodge degree zero.',witnesses=[])
        return r

    def test_deferred_orders_and_exact_hodge_are_preserved(self):
        r=public_component(self.entry())
        self.assertIsNone(r['group_order']);self.assertTrue(r['isotrivial'])
        self.assertEqual(r['isotriviality']['method'],'full-orbit-hodge-degree-v1')
        self.assertIs(r['block_quotients_complete'],False)

    def test_missing_order_requires_explicit_status(self):
        r=self.entry();r.pop('group_order_status')
        with self.assertRaises(ValueError):public_component(r)
        r=self.entry();r['block_quotients'][0].pop('orders_status')
        with self.assertRaises(ValueError):public_component(r)

    def test_false_hodge_value_or_sign_is_rejected(self):
        for field,value in [('twelve_L','1/2'),('L','1/24')]:
            r=self.entry();r['hodge_degree'][field]=value
            with self.assertRaises(ValueError):public_component(r)
        r=self.entry();r['isotrivial']=False;r['isotriviality']['status']='non-isotrivial'
        with self.assertRaises(ValueError):public_component(r)

    def test_distribution_mass_and_cycle_lengths_checked(self):
        r=self.entry();r['hodge_degree']['pair_product_cycle_distribution'][0]['states']=11
        with self.assertRaises(ValueError):public_component(r)
        r=self.entry();r['hodge_degree']['pair_product_cycle_distribution'][0]['cycle_lengths']=[8]
        with self.assertRaises(ValueError):public_component(r)

if __name__=='__main__':unittest.main()
