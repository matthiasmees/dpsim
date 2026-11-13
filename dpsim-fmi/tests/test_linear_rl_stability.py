import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "examples" / "linear_rl_stability"))
from run_benchmark import cs_multiplier, stability_limit


class StabilityAnalysisTest(unittest.TestCase):
    def test_boundary_is_40_13_microseconds(self):
        self.assertAlmostEqual(stability_limit() * 1e6, 40.134143, places=4)

    def test_multiplier_crosses_minus_one_at_boundary(self):
        self.assertAlmostEqual(cs_multiplier(stability_limit()), -1.0, places=12)

    def test_steps_below_and_above_boundary(self):
        self.assertLess(abs(cs_multiplier(40e-6)), 1.0)
        self.assertGreater(abs(cs_multiplier(40.25e-6)), 1.0)
        self.assertGreater(abs(cs_multiplier(50e-6)), 1.0)


if __name__ == "__main__":
    unittest.main()
