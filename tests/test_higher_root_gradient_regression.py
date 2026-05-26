import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "higher_root_gradient" / "h2o_bhhlyp_631gstar_nstate4_higher_root_gradient.json"


class HigherRootGradientRegressionTests(unittest.TestCase):
    @unittest.expectedFailure
    def test_live_diagnostic_fixture_marks_higher_root_mismatches_as_failures(self):
        """Stored FD diagnostic reproduces the remaining MRSF higher-root mismatch."""
        data = json.loads(FIXTURE.read_text())

        tolerance = 5.0e-4
        failing = []
        for case in data["cases"]:
            if not case["should_satisfy_gradient_tolerance"]:
                continue
            if case["max_abs_diff"] > tolerance:
                failing.append(
                    f"{case['method']} root {case['root']}: "
                    f"max_abs_diff={case['max_abs_diff']:.9f}, "
                    f"rms_diff={case['rms_diff']:.9f}, "
                    f"0_z ratio={case['ratio_fd_over_analytic']}"
                )

        self.assertEqual(
            [],
            failing,
            "Higher-root analytic gradients must agree with central finite differences "
            f"within {tolerance}; current mismatches: " + "; ".join(failing),
        )

    def test_tdhf_gradient_uses_target_state_transition_vectors(self):
        """TDDFT gradient must not flatten all roots and silently reuse root-1 X/Y."""
        source = (ROOT / "source" / "modules" / "tdhf_gradient.F90").read_text()

        self.assertRegex(
            source,
            r"xpy\s*\(:\s*,\s*:\s*\).*xmy\s*\(:\s*,\s*:\s*\)",
            "OQP_td_xpy/OQP_td_xmy are stored as [lexc, nstates]; the gradient "
            "reader must keep them rank-2 so target_state can select the requested root.",
        )
        self.assertRegex(
            source,
            r"xpy\s*\(:\s*,\s*infos%tddft%target_state\s*\)",
            "TDDFT gradient must build transition-density terms from the requested "
            "target_state, not from the first flat column of OQP_td_xpy.",
        )
        self.assertRegex(
            source,
            r"xmy\s*\(:\s*,\s*infos%tddft%target_state\s*\)",
            "TDDFT gradient must build transition-density terms from the requested "
            "target_state, not from the first flat column of OQP_td_xmy.",
        )

    @unittest.expectedFailure
    def test_mrsf_z_vector_operator_has_mrsf_specific_lhs(self):
        """MRSF Z-vector linear solve should not use only the SF LHS operator."""
        source = (ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90").read_text()
        sf_lhs_calls = len(re.findall(r"call\s+sfrolhs\s*\(", source, flags=re.IGNORECASE))
        mrsf_operator_calls = len(
            re.findall(r"call\s+mrsf(?:mntoia|esum|rolhs)\s*\(", source, flags=re.IGNORECASE)
        )

        self.assertGreater(
            mrsf_operator_calls,
            0,
            "MRSF Z-vector RHS and gradient include MRSF spin-pair terms, but the "
            f"current linear solve has {sf_lhs_calls} SF sfrolhs call(s) and no "
            "MRSF-specific LHS/operator application. Add/route through an MRSF "
            "operator before considering the higher-root MRSF gradient fixed.",
        )


if __name__ == "__main__":
    unittest.main()
