import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Z_VECTOR = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"


class MrsfSpcOperatorConsistencyTests(unittest.TestCase):
    def test_z_vector_ovov_only_control_maps_to_channel5(self):
        """Document the current OV-OV-only Z-vector control seam.

        Live CH2O/H2O/H2S FD checks can now fail even when the channel-7 ball
        density is preserved.  The next diagnostic is to compare all-SPC-on
        against OV-OV-only behavior, so keep the source-level mapping explicit:
        `spc_ovov` must scale only the channel-5 MRSF density in the Z-vector
        operator, while CO-CO and CO-OV controls map to their separate channel
        families.
        """
        source = Z_VECTOR.read_text()
        self.assertRegex(
            source,
            r"fmrst2\(:,\s*6,\s*:,\s*:\)\s*=\s*fmrst2\(:,\s*6,\s*:,\s*:\)\s*\*\s*infos%tddft%spc_coco\s*/\s*infos%tddft%hfscale",
            "CO-CO SPC should scale only Z-vector channel 6.",
        )
        self.assertRegex(
            source,
            r"fmrst2\(:,\s*5,\s*:,\s*:\)\s*=\s*fmrst2\(:,\s*5,\s*:,\s*:\)\s*\*\s*infos%tddft%spc_ovov\s*/\s*infos%tddft%hfscale",
            "OV-OV SPC should scale only Z-vector channel 5 for OV-OV-only diagnostics.",
        )
        self.assertRegex(
            source,
            r"fmrst2\(:,\s*1:4,\s*:,\s*:\)\s*=\s*fmrst2\(:,\s*1:4,\s*:,\s*:\)\s*\*\s*infos%tddft%spc_coov\s*/\s*infos%tddft%hfscale",
            "CO-OV SPC should scale only Z-vector channels 1:4.",
        )
        self.assertNotRegex(
            source,
            r"spc_ovov[^\n]*(?:fmrst2\(:,\s*[1467]|td_mrsf_den\(\s*[1467])",
            "OV-OV diagnostics should not silently rescale non-OV-OV channels.",
        )

    def test_gradient_spcscale_order_matches_reported_coco_ovov_coov_controls(self):
        """Keep gradient pair-class controls aligned with input/report order."""
        source = GRADIENT.read_text()
        self.assertRegex(source, r"\|\s+CO-CO\s+\|\s+OV-OV\s+\|\s+CO-OV\s+\|")
        self.assertRegex(
            source,
            r"spcscale\s*=\s*\[\s*infos%tddft%spc_coco,\s*&\s*\n\s*infos%tddft%spc_ovov,\s*&\s*\n\s*infos%tddft%spc_coov\]",
            "Gradient spcscale order must remain [CO-CO, OV-OV, CO-OV].",
        )
        self.assertRegex(source, r"qfspcp1\s*=\s*this%spcscale\(1\)")
        self.assertRegex(source, r"qfspcp2\s*=\s*this%spcscale\(2\)")
        self.assertRegex(source, r"qfspcp3\s*=\s*this%spcscale\(3\)")
        self.assertRegex(source, r"df1\s*=\s*df1\s*\+\s*sgnk\*qfspcp1\*db1")
        self.assertRegex(source, r"df1\s*=\s*df1\s*-\s*sgnk\*qfspcp2\*db2")
        self.assertRegex(source, r"df1\s*=\s*df1\s*\+\s*sgnk\*qfspcp3\*\(-dc1-dc2-dc3-dc4\s*&\s*\n\s*\+dd1\+dd2\+dd3\+dd4\)")


if __name__ == "__main__":
    unittest.main()
