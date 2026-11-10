try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest

from verify_genus_three_saturated import arithmetic_rules,violations


def partitions(n,minimum=1,maximum=None):
    if n==0:yield();return
    for k in range(minimum,min(n,maximum if maximum is not None else n)+1):
        for tail in partitions(n-k,k,maximum):yield(k,)+tail


class SaturatedTests(unittest.TestCase):
    def test_arithmetic_scopes_and_norm_versus_root_run(self):
        self.assertIn(('four-tile-prime-periodicity',2,4),arithmetic_rules(17,4))
        self.assertIn(('four-tile-prime-periodicity',2,4),arithmetic_rules(23,4))
        self.assertEqual(arithmetic_rules(9,4),set())
        self.assertEqual(arithmetic_rules(31,4),{('four-tile-prime-periodicity',2,4)})
        self.assertEqual(arithmetic_rules(31,5),{('prime-tile-valuation',5,5)})

    def test_designated_partner_is_not_equalized(self):
        shapes=[(4,)+(1,)*19,(5,5,5,8),(5,5,5,5,1,1,1),(5,5,5)+(1,)*8]
        self.assertEqual(violations(23,shapes),set())
        shapes[2]=(2,2,2,2,5,5,5)
        self.assertTrue(violations(23,shapes))

    def test_prime_five_keeps_different_nonzero_residues(self):
        shape=(2,)*8+(3,)*5
        shapes=[(5,)+(1,)*26,(6,6,6,6,7),shape,shape]
        self.assertFalse(violations(31,shapes))
        shapes[2]=(1,)*8+(2,3,4,5,9)
        self.assertTrue(any(reason=='ordinary-cycle-divisibility'for *_,reason,v in violations(31,shapes)))

    def test_small_single_three_residue_slices_exhaust_integer_partitions(self):
        expected={17:{tuple(sorted(((2,5,5,5),(1,1,1,1,1,4,4,4)))),tuple(sorted(((1,4,4,4,4),(2,2,2,2,2,2,5))))},
                  19:{((1,1,1,4,4,4,4),(1,1,1,4,4,4,4))}}
        for n in(17,19):
            all_parts=list(partitions(n,maximum=n//3));kept=set()
            for b in all_parts:
                for e in all_parts:
                    if len(b)+len(e)!=n-5:continue
                    shapes=[(3,)+(1,)*(n-3),(1,1,n-2),b,e]
                    if not violations(n,shapes):kept.add(tuple(sorted((b,e))))
            self.assertEqual(kept,expected[n])


if __name__=='__main__':unittest.main()
