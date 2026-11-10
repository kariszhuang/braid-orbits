try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
from math import gcd
import unittest
from analyze_fixed_braid_constraints import product, point_orbits
from odd_order_spin import hurwitz, cycle_lengths, inverse
from run_smooth_mission import values
from run_smooth_special import construct, recover_blocks, JOBS, key_digest, known_folders
from run_s18_linked_roots import state_key


class SpecialChecks(unittest.TestCase):
    def test_all_small_voltages(self):
        tested = 0
        for r, k in [(3, 3), (3, 4), (5, 3), (5, 4), (7, 3)]:
            n = r*k
            for q in range(1, r):
                if gcd(q, r) != 1 or gcd(q+1, r) != 1:
                    continue
                for y in itertools.permutations(range(k)):
                    t = construct(r, k, q, y)
                    self.assertEqual(product(*t), tuple(range(n)))
                    self.assertEqual(len(point_orbits(t)), 1)
                    if t[0] == tuple(range(n)):
                        continue
                    tested += 1
                    self.assertIsNotNone(recover_blocks(t, r, k))
                    self.assertTrue(all(cycle_lengths(p) == [r]*k for p in t[1:]))
                    for m in range(6):
                        h = hurwitz(t, m); v = values(h)
                        self.assertIsNotNone(recover_blocks(h, r, k))
                        self.assertTrue(v['smooth']); self.assertEqual(v['defect'], 0)
        self.assertGreater(tested, 30)

    def test_higher_genus_and_new_small_type(self):
        # A disjoint-transposition commutator in S4, not a single3-cycle.
        t = construct(3, 4, 1, (2, 1, 0, 3))
        self.assertEqual(sorted(cycle_lengths(t[0])), [1]*8+[2, 2])
        self.assertEqual(values(t)['genus'], 2)
        t = construct(5, 3, 1, (1, 0, 2))
        self.assertEqual(values(t)['genus'], 5)
        self.assertTrue(values(t)['smooth'])

    def test_relabeling_and_digest(self):
        t = construct(5, 3, 1, (1, 0, 2)); h = tuple(reversed(range(15)))
        relabeled = tuple(product(h, p, inverse(h)) for p in t)
        self.assertIsNotNone(recover_blocks(relabeled, 5, 3))
        self.assertEqual(state_key(t), state_key(relabeled))
        self.assertEqual(key_digest([state_key(t)]), key_digest([state_key(relabeled)]))

    def test_plan_excludes_completed_template(self):
        self.assertEqual(len(JOBS), 9)
        self.assertNotIn(dict(blocks=3, size=3), JOBS)
        self.assertTrue(all(j['blocks']*j['size'] <= 25 for j in JOBS))

    def test_known_folder_filter_excludes_audit_json(self):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root/'block_family_12'; folder.mkdir()
            (folder/'SUMMARY.tsv').write_text('test fixture\n')
            (root/'block_family_12_CHECKS.json').write_text('{}\n')
            self.assertEqual(known_folders(root), [folder])


if __name__ == '__main__':
    unittest.main()
