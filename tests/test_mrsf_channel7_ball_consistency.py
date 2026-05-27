import math
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Z_VECTOR = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
MRSF_LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"


class MrsfChannel7BallConsistencyTests(unittest.TestCase):
    def test_z_vector_preserves_mrsfcbc_channel7_ball_density(self):
        """MRSF Z-vector must not replace mrsfcbc's channel-7 ball density.

        mrsfcbc builds channel 7 (ball) with explicit open-open spin-adapted
        terms.  For CH2O S2/root 3, the dominant DRF has strong open-open
        [7,8] character, so overwriting ball with sfdmat(mrsfxvec(...)) loses
        the density actually used by the MRSF channel decomposition.
        """
        source = Z_VECTOR.read_text()

        self.assertRegex(
            source,
            r"call\s+mrsfcbc\s*\([^\n]*fmrst1\(1,\s*:\s*,\s*:\s*,\s*:\s*\)\)",
            "Z-vector should build all seven MRSF channel densities with mrsfcbc.",
        )
        self.assertNotRegex(
            source,
            r"fmrst1\s*\(\s*1\s*,\s*7\s*,\s*:\s*,\s*:\s*\)\s*=\s*td_abxc",
            "Do not overwrite mrsfcbc channel-7 ball with td_abxc; the two are "
            "not equivalent for deliberate open-open MRSF amplitudes.",
        )

    def test_mrsfcbc_documents_open_open_channel7_terms(self):
        """Guard the open-open terms that make mrsfcbc ball distinct."""
        source = MRSF_LIB.read_text()
        body = re.search(
            r"subroutine\s+mrsfcbc\b(?P<body>.*?)end\s+subroutine\s+mrsfcbc",
            source,
            re.IGNORECASE | re.DOTALL,
        )
        self.assertIsNotNone(body, "mrsfcbc subroutine must exist")
        assert body is not None
        text = body.group("body")

        self.assertRegex(text, r"ball\s*=>\s*fmrsf\(\s*7\s*,\s*:\s*,\s*:\s*\)")
        self.assertRegex(text, r"ball\s*=\s*ball\s*\+\s*bo2v\s*\+\s*bo1v\s*\+\s*bco1\s*\+\s*bco2")
        self.assertRegex(text, r"bvec\(lr1,\s*lr1\)\s*\*\s*isqrt2")
    def test_open_open_weight_makes_mrsfxvec_sfdmat_differ_from_mrsfcbc_ball(self):
        """Dependency-light numeric check for the channel-7 consistency hazard."""
        nbf = 4
        noca = 3
        nocb = 1
        lr1 = noca - 1  # zero-based equivalent of Fortran lr1=nocca-1
        lr2 = noca      # zero-based equivalent of Fortran lr2=nocca
        oo_weight = 0.5

        # Fortran mrsfxvec triplet mapping for the open-open diagonal:
        # X(lr1,lr1) -> sqrt(2)*X in both open-open diagonal slots.
        sfdmat_ball = [[0.0 for _ in range(nbf)] for _ in range(nbf)]
        sfdmat_ball[lr1][lr1] = math.sqrt(2.0) * oo_weight
        sfdmat_ball[lr2][lr2] = math.sqrt(2.0) * oo_weight

        # Fortran mrsfcbc triplet channel-7 ball term for the same original
        # open-open amplitude with identity alpha/beta MO coefficients:
        # (C_lr1 C_lr1 + C_lr2 C_lr2) * X(lr1,lr1) / sqrt(2).
        mrsfcbc_ball = [[0.0 for _ in range(nbf)] for _ in range(nbf)]
        mrsfcbc_ball[lr1][lr1] = oo_weight / math.sqrt(2.0)
        mrsfcbc_ball[lr2][lr2] = oo_weight / math.sqrt(2.0)

        self.assertNotEqual(sfdmat_ball, mrsfcbc_ball)
        self.assertAlmostEqual(2.0, sfdmat_ball[lr1][lr1] / mrsfcbc_ball[lr1][lr1])
        self.assertAlmostEqual(2.0, sfdmat_ball[lr2][lr2] / mrsfcbc_ball[lr2][lr2])


if __name__ == "__main__":
    unittest.main()
