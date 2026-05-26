import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"


class MrsfZVectorOperatorTests(unittest.TestCase):
    @unittest.expectedFailure
    def test_gmres_operator_uses_mrsf_specific_lhs_not_sf_lhs(self):
        """MRSF GMRES route must not apply the SF response operator as its LHS."""
        source = SOURCE.read_text()
        match = re.search(
            r"subroutine\s+apply_z_operator\b(?P<body>.*?)end\s+subroutine\s+apply_z_operator",
            source,
            flags=re.IGNORECASE | re.DOTALL,
        )
        self.assertIsNotNone(match, "Expected an explicit MRSF GMRES operator routine")
        body = match.group("body")

        self.assertNotRegex(
            body,
            r"\bsfrolhs\s*\(",
            "MRSF z-vector GMRES must not reuse the SF LHS operator; root-3 FD "
            "diagnostics show SF and MRSF diverge, so the MRSF operator/preconditioner "
            "needs its own validated route.",
        )
        self.assertRegex(
            body,
            r"\bmrsf(?:cbc|mntoia|esum|sp)\s*\(",
            "MRSF z-vector GMRES operator should build/apply MRSF spin-pair response "
            "terms rather than only TD/SF A+B pieces.",
        )
    def test_target_state_mrsf_density_uses_alpha_and_beta_orbitals(self):
        """MRSF density construction should match energy code's mo_a/mo_b usage."""
        source = SOURCE.read_text()
        self.assertIn(
            "call mrsfcbc(infos, mo_a, mo_b, wrk1, fmrst1(1,:,:,:))",
            source,
            "MRSF z-vector RHS/density setup should use beta MO coefficients for "
            "the beta side, consistent with tdhf_mrsf_energy.F90; using mo_a twice "
            "can make the gradient route inconsistent with the response eigenvectors.",
        )
        self.assertIn(
            "call mrsfsp(hxa, hxb, mo_a, mo_b, wrk3, fmrst2(1,:,:,:), nocca, noccb)",
            source,
            "MRSF spin-pair gradient terms should also use alpha/beta MO coefficients "
            "rather than mo_a for both sides.",
        )


if __name__ == "__main__":
    unittest.main()
