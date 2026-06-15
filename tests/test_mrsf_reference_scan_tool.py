import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_TOOL = ROOT / "tools" / "mrsf_reference_scan.py"


def load_scan_tool():
    spec = importlib.util.spec_from_file_location("mrsf_reference_scan_under_test", SCAN_TOOL)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TestMrsfReferenceScanTool(unittest.TestCase):
    def setUp(self):
        self.scan = load_scan_tool()

    def test_render_gap_softmax_input_contains_reference_controls(self):
        text = self.scan.render_input(
            self.scan.h2o_triplet_geometry(1.0),
            self.scan.VARIANTS["gap_softmax"],
        )

        self.assertIn("method=tdhf", text)
        self.assertIn("maxit=60", text)
        self.assertIn("conv=1.0e-6", text)
        self.assertIn("[mrsf_ref]", text)
        self.assertIn("mode=ensemble", text)
        self.assertIn("open_pairs=auto", text)
        self.assertIn("weights=gap_softmax", text)
        self.assertIn("weight_temperature=0.05", text)
        self.assertIn("max_refs=6", text)

    def test_render_input_accepts_manual_open_pairs(self):
        text = self.scan.render_input(
            self.scan.ethylene_torsion_geometry(90.0),
            self.scan.VARIANTS["equal"],
            open_pairs="8:9;7:10",
            max_refs=4,
        )

        self.assertIn("open_pairs=8:9;7:10", text)
        self.assertIn("max_refs=4", text)

    def test_render_rohf_input_avoids_mrsf_reference_section(self):
        text = self.scan.render_input(
            self.scan.h2o_triplet_geometry(0.98),
            self.scan.VARIANTS["rohf"],
        )

        self.assertIn("method=hf", text)
        self.assertNotIn("[mrsf_ref]", text)
        self.assertNotIn("[tdhf]", text)

    def test_render_mrsf_input_uses_single_reference_response(self):
        text = self.scan.render_input(
            self.scan.ethylene_torsion_geometry(90.0),
            self.scan.VARIANTS["mrsf"],
        )

        self.assertIn("method=tdhf", text)
        self.assertIn("[tdhf]", text)
        self.assertIn("type=mrsf", text)
        self.assertIn("maxit=60", text)
        self.assertNotIn("[mrsf_ref]", text)

    def test_ethylene_torsion_rotates_ch2_groups_symmetrically(self):
        planar = self.scan.ethylene_torsion_geometry(0.0)
        twisted = self.scan.ethylene_torsion_geometry(90.0)

        self.assertAlmostEqual(planar[2][1], 0.0)
        self.assertAlmostEqual(planar[4][1], 0.0)
        self.assertLess(twisted[2][1], -0.6)
        self.assertGreater(twisted[4][1], 0.6)
        self.assertEqual(twisted[0], planar[0])
        self.assertEqual(twisted[1], planar[1])

    def test_parse_log_extracts_final_scf_and_applied_weights(self):
        log = """
   PyOQP MRSF pair selection:          auto/frontier_window
   PyOQP MRSF weight model:            gap_softmax
   PyOQP MRSF weight temperature (Eh): 0.05
   PyOQP MRSF reference open pairs:    [[5, 6], [4, 7]]
   PyOQP MRSF reference weights:       [0.886, 0.114]
   PyOQP MRSF SCF applied pairs:       [[5, 6], [4, 7]]
   PyOQP MRSF SCF applied weights:     [0.987, 0.013]
   PyOQP MRSF min frontier gap (Eh):   0.050
   PyOQP MRSF response status:         implemented_energy_only
   PyOQP MRSF response model:          state_interaction_overlap
   PyOQP MRSF response coupled:        yes
   PyOQP MRSF full response kernel:    no
   PyOQP MRSF response energy only:    yes
   PyOQP MRSF selected states:         [{'energy': -0.0125, 'rank': 1, 'dominant_open_pair': [5, 6]}]
   PyOQP MRSF candidate states:        2
   PyOQP MRSF raw candidate states:    3
   PyOQP MRSF skipped blocks:          [2]
   PyOQP MRSF SI common dimension:     42
   PyOQP: SCF not converged; escalating to soscf
 Final ROHF energy is      -75.1000000000 after  7 iterations
          SCF convergence achieved ....
 Final ROHF energy is      -75.2000000000 after  9 iterations
   MRSF-TD-DFT energies converged in    4 iterations
   MRSF-TD-DFT energies converged in    5 iterations
   PyOQP state 0      -75.20000000
   PyOQP state 1      -75.21250000
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.log"
            path.write_text(log)

            parsed = self.scan.parse_log(path)

        self.assertEqual(parsed["pair_selection"], "auto/frontier_window")
        self.assertEqual(parsed["weight_model"], "gap_softmax")
        self.assertEqual(parsed["open_pairs"], [[5, 6], [4, 7]])
        self.assertEqual(parsed["reference_weights"], [0.886, 0.114])
        self.assertEqual(parsed["applied_weights"], [0.987, 0.013])
        self.assertEqual(parsed["scf_energy"], -75.2)
        self.assertEqual(parsed["scf_iterations"], 9)
        self.assertTrue(parsed["scf_converged"])
        self.assertTrue(parsed["scf_escalated"])
        self.assertEqual(parsed["min_frontier_gap_hartree"], 0.05)
        self.assertEqual(parsed["response_energy"], -0.0125)
        self.assertEqual(parsed["state_energy"], -75.2125)
        self.assertEqual(parsed["mrsf_converged_blocks"], 2)
        self.assertEqual(parsed["mrsf_block_iterations"], [4, 5])
        self.assertEqual(parsed["response_status"], "implemented_energy_only")
        self.assertEqual(parsed["response_model"], "state_interaction_overlap")
        self.assertTrue(parsed["response_coupled"])
        self.assertFalse(parsed["full_response_kernel"])
        self.assertTrue(parsed["response_energy_only"])
        self.assertEqual(parsed["response_candidate_count"], 2)
        self.assertEqual(parsed["response_raw_candidate_count"], 3)
        self.assertEqual(parsed["response_skipped_blocks"], [2])
        self.assertEqual(parsed["response_si_common_dimension"], 42)
        self.assertEqual(parsed["dominant_open_pair"], [5, 6])


if __name__ == "__main__":
    unittest.main()
