try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import struct
from pathlib import Path
import tempfile
import unittest

from odd_order_spin import parse_seed, hurwitz
from run_smooth_mission import values
from run_smooth_degree9 import binary_records


class DegreeNineChecks(unittest.TestCase):
    def test_compact_records_and_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'records'
            path.write_bytes(struct.pack('<IIII', 10, 4, 362879, 25))
            self.assertEqual(list(binary_records(path)), [(10, 4), (362879, 25)])
            path.write_bytes(path.read_bytes()+b'\x00')
            with self.assertRaisesRegex(ValueError, 'truncated'):
                list(binary_records(path))

    def test_known_positive_genus_control_and_all_neighbors(self):
        t = parse_seed('[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]', 9)
        from run_s18_linked_roots import state_key
        seen = {state_key(t)}; queue = [t]
        for u in queue:
            v = values(u)
            self.assertTrue(v['smooth']); self.assertEqual(v['defect'], 0); self.assertEqual(v['genus'], 2)
            for move in range(6):
                w = hurwitz(u, move); key = state_key(w)
                if key not in seen:
                    seen.add(key); queue.append(w)
        self.assertEqual(len(seen), 12)


if __name__ == '__main__':
    unittest.main()
