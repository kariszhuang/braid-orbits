"""Small checks for exact spectral certification and the regression rank guard."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass

import unittest

import numpy as np
import pandas as pd

from analyze_spectral_followups import design_matrix
from run_exact_spectral_certificates import adjacency_from_code, exact_certificate


class SpectralTheoremLogicTests(unittest.TestCase):
    def test_two_vertex_certificate_and_invalid_code(self) -> None:
        # Six labeled moves join the two vertices; spectrum is {6, -6}.
        code = ",".join(["1"] * 6 + ["0"] * 6)
        result = exact_certificate(code, 2)
        self.assertEqual(result["target_multiplicity"], 0)
        self.assertFalse(result["lambda2_certified"])
        with self.assertRaises(ValueError):
            adjacency_from_code("0,1", 2)

    def test_collinear_design_is_rejected(self) -> None:
        rng = np.random.default_rng(7)
        n = 40
        frame = pd.DataFrame({
            "group_degree": ["A5"] * 20 + ["S5"] * 20,
            "log_vertices": rng.normal(size=n),
            "loop_rate": rng.random(n),
            "parallel_rate": rng.random(n),
            "genus": np.arange(n) % 4,
        })
        frame["generator_cycle_density"] = frame["loop_rate"]
        with self.assertRaisesRegex(ValueError, "nonidentifiable design"):
            design_matrix(frame, True)


if __name__ == "__main__":
    unittest.main()
