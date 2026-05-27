import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TDDFT_GRADIENT = ROOT / "source" / "modules" / "tdhf_gradient.F90"
MRSF_GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
XC_GRADIENT = ROOT / "source" / "dftlib" / "dft_gridint_tdxc_grad.F90"


def _call_block(source: str, callee: str) -> str:
    start = re.search(rf"call\s+{callee}\s*\(", source, re.IGNORECASE)
    if not start:
        raise AssertionError(f"Could not find call to {callee}")
    tail = source[start.start():]
    end = re.search(r"\n\s*infos\s*=\s*infos\s*\)", tail, re.IGNORECASE)
    if not end:
        raise AssertionError(f"Could not find end of {callee} keyword call")
    return tail[: end.end()]


class MrsfXcGradientDensityDiagnostics(unittest.TestCase):
    def test_conventional_tddft_passes_root_selected_transition_density_to_xc_gradient(self):
        """Conventional TDDFT enables dF_xc/dR by passing root-selected X+Y."""
        source = TDDFT_GRADIENT.read_text()
        self.assertRegex(source, r"xpy\(:,\s*infos%tddft%target_state\)")
        self.assertRegex(source, r"xpy2\(:,:,1\)")
        block = _call_block(source, "tddft_xc_gradient")
        self.assertRegex(block, r"xa\s*=\s*xpy2\(:,:,1:1\)")

    def test_unrestricted_xc_gradient_uses_optional_xa_xb_to_enable_grad_x(self):
        """The unrestricted XC gradient computes grad_X only when xa/xb are present."""
        source = XC_GRADIENT.read_text()
        self.assertRegex(source, r"doFxc\s*=\s*present\(xa\)")
        self.assertRegex(source, r"if \(doFxc\) nxcder\s*=\s*3")
        self.assertRegex(source, r"dat%xa\s*=>\s*xa")
        self.assertRegex(source, r"dat%xb\s*=>\s*xb")
        self.assertRegex(source, r"grad_X")

    def test_mrsf_currently_omits_xc_transition_density_handoff(self):
        """Document the MRSF/TDDFT XC-gradient path difference.

        A RED version of this diagnostic required `xa=v(:,:,1:1)` and failed
        on the pre-fix source.  The direct trial handoff (`v(:,:,2)=td_abxc`,
        `xa=v(:,:,1:1)`, `xb=v(:,:,2:2)`) was then live-tested at commit
        ca26976 and worsened/failed the targeted MRSF FD checks, so the
        production change was reverted.  This source-level diagnostic keeps the
        concrete discrepancy visible without blessing the rejected spin-channel
        mapping as a fix.
        """
        source = MRSF_GRADIENT.read_text()
        self.assertRegex(source, r"call\s+tagarray_get_data\(infos%dat,\s*OQP_td_abxc,\s*td_abxc\)")
        self.assertRegex(source, r"v\(:,:,1\)\s*=\s*td_abxc")
        block = _call_block(source, "utddft_xc_gradient")
        self.assertNotRegex(block, r"xa\s*=")
        self.assertNotRegex(block, r"xb\s*=")
        self.assertNotRegex(source, r"v\(:,:,2\)\s*=\s*td_abxc")


if __name__ == "__main__":
    unittest.main()
