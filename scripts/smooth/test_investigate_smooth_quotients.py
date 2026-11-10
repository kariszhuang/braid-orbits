try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from investigate_smooth_quotients import quotient, two_active_blocks, genus_zero_quotient, index, smooth, profile_reports, index_regime
from compact_smooth_outputs import compact
from verify_smooth_geometry import parse_tuple, cycles, compose, orbits


class QuotientTests(unittest.TestCase):
    def test_star_triangle_and_impossible_passport_budgets(self):
        self.assertEqual(index_regime([4, 12, 13, 13], 18), 'star-budget')
        self.assertEqual(index_regime([2, 2, 2, 4], 5), 'triangle-budget')
        with self.assertRaisesRegex(ValueError, 'cannot be smooth'):
            index_regime([3, 3, 3, 3], 5)

    def test_reported_degree_nine(self):
        t = parse_tuple('[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]', 9)
        q = quotient(t, [[4, 5, 6], [0, 7, 8], [1, 2, 3]])
        self.assertEqual(q['quotient_genus'], 1)
        self.assertEqual(q['quotient_indices'], [0, 2, 2, 2])
        self.assertEqual(q['genus'], 2)
        self.assertEqual(q['ramification_excess'], [2, 0, 0, 0])

    def test_multi_local_frame_and_population(self):
        count = 0
        for t, dtype in two_active_blocks():
            self.assertTrue(smooth(t)); self.assertEqual(index(compose(t[0], t[1])), 10)
            self.assertEqual(sorted(map(len, orbits(t[:2]))), [3, 9])
            count += 1
        self.assertEqual(count, 576)

    def test_rational_quotient_frame(self):
        for t, dtype in itertools.islice(genus_zero_quotient(), 20):
            self.assertTrue(smooth(t))
            self.assertEqual(index(t[0]), 4); self.assertEqual(index(t[1]), 6)
            self.assertEqual(sorted(map(len, cycles(t[2]))), [2, 2, 2, 2, 4])
            self.assertEqual(len(orbits(t)), 1)

    def test_compacted_profile_is_available_without_restoring(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp); root = base / 'run'; root.mkdir()
            folder = root / 'profile_000000'; folder.mkdir()
            path = folder / 'PROFILES.json'
            payload = dict(start=0, next=1, components=[dict(genus=2)])
            path.write_text(json.dumps(payload))
            out = base / 'packed'
            compact(root, out, [path], [folder])
            reports = list(profile_reports(root, out / 'REPEATED_OUTPUTS.tar.gz'))
            self.assertEqual(len(reports), 1)
            self.assertEqual(reports[0][0], payload)
            self.assertEqual(reports[0][1]['path'], str(path))
            self.assertFalse(path.exists())

    def test_invalid_partition_is_not_silently_reinterpreted(self):
        t = parse_tuple('[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]', 9)
        with self.assertRaisesRegex(ValueError, 'not invariant'):
            quotient(t, [[0, 1, 2], [3, 4, 5], [6, 7, 8]])

    def test_additivity_descends_in_three_two_point_blocks(self):
        group = []
        for q in itertools.permutations(range(3)):
            for offsets in itertools.product(range(2), repeat=3):
                group.append(tuple(2*q[x//2]+((x%2)^offsets[x//2]) for x in range(6)))
        for a, b in itertools.product(group, repeat=2):
            if index(compose(a, b)) != index(a) + index(b):
                continue
            qa = tuple(a[2*j]//2 for j in range(3))
            qb = tuple(b[2*j]//2 for j in range(3))
            self.assertEqual(index(compose(qa, qb)), index(qa) + index(qb))

    def test_multiple_return_vertices_forest_criterion(self):
        # Base pair (12),(23) has one incidence tree. Its 3-sheet lifts have
        # arbitrary returns alpha,beta at two vertices; all other returns id.
        for alpha, beta in itertools.product(itertools.permutations(range(3)), repeat=2):
            a = tuple((3+x if x < 3 else alpha[x-3] if x < 6 else x) for x in range(9))
            b = tuple((x if x < 3 else 3+x if x < 6 else 3+beta[x-6]) for x in range(9))
            excess = index(alpha) + index(beta)
            components = len(orbits((alpha, beta)))
            rank = excess - 3 + components
            self.assertGreaterEqual(rank, 0)
            self.assertEqual(len(orbits((a, b))), components)
            additive = index(compose(a, b)) == index(a) + index(b)
            self.assertEqual(additive, rank == 0)


if __name__ == '__main__': unittest.main()
