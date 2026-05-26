from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
MRSF_GRAD = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"


class MrsfOvovGradientSignTest(unittest.TestCase):
    def test_ovov_spin_pair_gradient_uses_positive_sign(self):
        src = MRSF_GRAD.read_text()
        self.assertIn("df1 = df1 + sgnk*qfspcp2*db2", src)
        self.assertNotIn("df1 = df1 - sgnk*qfspcp2*db2", src)


if __name__ == "__main__":
    unittest.main()
