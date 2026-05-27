import importlib.util
import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "gradient_state_character.py"


def load_module():
    spec = importlib.util.spec_from_file_location("gradient_state_character", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MRSF_LOG = textwrap.dedent(
    """
    Spin-adapted spin-flip excitations

     State #   1  Energy =   -2.992748 eV
                   <S^2> =    0.0000
            DRF    Coeff        OCC       VIR
            ---  --------     ------    ------
              9 -0.978591         9  ->     8

     State #   3  Energy =    5.868203 eV
                   <S^2> =    0.0000
            DRF    Coeff        OCC       VIR
            ---  --------     ------    ------
              7 -0.924997         7  ->     8
              9  0.147693         9  ->     8
             63  0.063459         9  ->    14

         Summary table
    """
)

TDDFT_LOG = textwrap.dedent(
    """
        Summary of the TD-DFT calculation

      Transition          dE(eV)          DX          DY          DZ      Osc.str.
       0 -> 1            4.187711      0.0000     -0.0000      0.0000      0.0000
       0 -> 4           10.188995      0.9131      0.0000     -0.0000      0.2081
       0 -> 5           10.783632      0.0000      0.0000      0.0000      0.0000
       0 -> 6           12.134527     -0.1418     -0.0000      0.0000      0.0060
     ==============================================================================
    """
)


class GradientStateCharacterDiagnosticsTests(unittest.TestCase):
    def test_parses_mrsf_spin_adapted_state_configurations(self):
        diagnostics = load_module()

        states = diagnostics.parse_mrsf_states(MRSF_LOG)

        root3 = states[3]
        self.assertEqual(root3["root"], 3)
        self.assertAlmostEqual(root3["energy_ev"], 5.868203)
        self.assertAlmostEqual(root3["s2"], 0.0)
        self.assertEqual(root3["top_configuration"]["occ"], 7)
        self.assertEqual(root3["top_configuration"]["vir"], 8)
        self.assertAlmostEqual(root3["top_configuration"]["weight"], 0.85561945, places=7)

    def test_parses_tddft_summary_transition_signature(self):
        diagnostics = load_module()

        transitions = diagnostics.parse_tddft_transitions(TDDFT_LOG)

        root4 = transitions[4]
        self.assertAlmostEqual(root4["energy_ev"], 10.188995)
        self.assertEqual(root4["dominant_dipole_axis"], "x")
        self.assertAlmostEqual(root4["dominant_dipole_abs"], 0.9131)
        self.assertAlmostEqual(root4["oscillator_strength"], 0.2081)

    def test_cli_reports_stable_formaldehyde_mrsf_signature_for_matching_perturbations(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            base = tmpdir / "base.log"
            plus = tmpdir / "plus.log"
            minus = tmpdir / "minus.log"
            base.write_text(MRSF_LOG)
            plus.write_text(MRSF_LOG.replace("5.868203", "5.864345"))
            minus.write_text(MRSF_LOG.replace("5.868203", "5.872063"))

            out = subprocess.check_output(
                [
                    sys.executable,
                    str(SCRIPT),
                    "mrsf",
                    "--target-root",
                    "3",
                    str(base),
                    str(plus),
                    str(minus),
                ],
                text=True,
            )

        data = json.loads(out)
        self.assertTrue(data["stable"])
        self.assertEqual(data["target_root"], 3)
        self.assertEqual(data["signatures"][0]["top_occ_vir"], [7, 8])
        self.assertLess(data["max_top_weight_delta"], 1e-12)


if __name__ == "__main__":
    unittest.main()
