try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import unittest
from math import factorial

from analyze_fixed_braid_constraints import product, point_orbits
from odd_order_spin import inverse, cycle_lengths
from run_smooth_mission import values
from run_smooth_small_census import audit_frames, pair_key_and_centralizer, partition


class SmallCensusChecks(unittest.TestCase):
    def test_exported_pair_frames_parser_and_population(self):
        import csv
        from pathlib import Path
        import tempfile
        n = 3; perms = [(p, n-len(cycle_lengths(p))) for p in itertools.permutations(range(n))]
        frames = {}
        for (a, ia), (b, ib) in itertools.product(perms, repeat=2):
            if ia+ib == n-len(cycle_lengths(product(a, b))):
                key, z = pair_key_and_centralizer(a, b)
                frames[key] = (a, b, factorial(n)//z)
        def text(p):
            seen = set(); cycles = []
            for root in range(n):
                if root in seen:
                    continue
                points = []; x = root
                while x not in seen:
                    seen.add(x); points.append(x+1); x = p[x]
                if len(points)>1:
                    cycles.append('('+' '.join(map(str, points))+')')
            return ''.join(cycles) or '()'
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            with (folder/'FRAMES.tsv').open('w') as f:
                writer = csv.writer(f, delimiter='\t'); writer.writerow(['frame', 'pair_orbit', 'a', 'b'])
                for i, (a, b, weight) in enumerate(frames.values()):
                    writer.writerow([i, weight, text(a), text(b)])
            self.assertEqual(len(audit_frames(folder, n, perms)), len(frames))

    def test_partitions_are_exact_and_unique(self):
        for n, count in [(3, 3), (4, 5), (5, 7), (6, 11), (7, 15), (8, 22)]:
            parts = list(partition(n))
            self.assertEqual(len(set(parts)), count)
            self.assertTrue(all(sum(p) == n and tuple(sorted(p)) == p for p in parts))

    def test_literal_centralizers_all_s3_pairs(self):
        perms = list(itertools.permutations(range(3)))
        for a, b in itertools.product(perms, repeat=2):
            key, z = pair_key_and_centralizer(a, b)
            literal = sum(product(h, a, inverse(h)) == a and product(h, b, inverse(h)) == b for h in perms)
            self.assertEqual(z, literal)
            for h in perms:
                self.assertEqual(pair_key_and_centralizer(product(h, a, inverse(h)), product(h, b, inverse(h))), (key, z))

    def test_every_literal_s4_tuple_matches_weighted_reduction(self):
        n = 4; perms = list(itertools.permutations(range(n)))
        frames = {}
        for a, b in itertools.product(perms, repeat=2):
            ip = n-len(cycle_lengths(product(a, b)))
            if 2*n-len(cycle_lengths(a))-len(cycle_lengths(b)) == ip:
                key, z = pair_key_and_centralizer(a, b)
                frames[key] = (a, b, factorial(n)//z)
        reduced = [0, 0]
        for a, b, weight in frames.values():
            for c in perms:
                t = (a, b, c, inverse(product(a, b, c)))
                if len(point_orbits(t)) != 1:
                    continue
                reduced[0] += weight
                v = values(t)
                self.assertTrue(v['smooth']); self.assertEqual(v['defect'], 0)
                indices = [n-len(cycle_lengths(p)) for p in t]
                if sum(indices[2:]) == v['product_index']:
                    reduced[1] += weight
        literal = [0, 0, 0]
        for a, b, c in itertools.product(perms, repeat=3):
            t = (a, b, c, inverse(product(a, b, c)))
            if len(point_orbits(t)) != 1:
                continue
            v = values(t); indices = [n-len(cycle_lengths(p)) for p in t]
            l = sum(indices[:2]) == v['product_index']; r = sum(indices[2:]) == v['product_index']
            literal[0] += l; literal[1] += l and r; literal[2] += l or r
            if v['smooth']:
                self.assertEqual(v['defect'], 0)
        self.assertEqual(literal, [*reduced, 2*reduced[0]-reduced[1]])


if __name__ == '__main__':
    unittest.main()
