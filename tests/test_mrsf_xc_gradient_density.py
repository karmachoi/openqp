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


class MrsfXcGradientDensityTests(unittest.TestCase):
    def test_conventional_tddft_documents_transition_density_xc_handoff(self):
        """Conventional TDDFT enables dF_xc/dR by passing root-selected X+Y."""
        source = TDDFT_GRADIENT.read_text()
        self.assertRegex(source, r"xpy\(:,\s*infos%tddft%target_state\)")
        self.assertRegex(source, r"xpy2\(:,:,1\)")
        block = _call_block(source, "tddft_xc_gradient")
        self.assertRegex(block, r"xa\s*=\s*xpy2\(:,:,1:1\)")

    def test_utddft_xc_gradient_requires_optional_xa_xb_to_enable_fxc(self):
        """The unrestricted XC gradient computes grad_X only when xa/xb are present."""
        source = XC_GRADIENT.read_text()
        self.assertRegex(source, r"doFxc\s*=\s*present\(xa\)")
        self.assertRegex(source, r"if \(doFxc\) nxcder\s*=\s*3")
        self.assertRegex(source, r"dat%xa\s*=>\s*xa")
        self.assertRegex(source, r"dat%xb\s*=>\s*xb")

    def test_mrsf_xc_gradient_passes_td_abxc_as_transition_density(self):
        """MRSF must not drop the target-root transition density in XC gradient.

        `td_abxc` is built in the MRSF Z-vector path as the target-root
        density-like `b=A*x` handoff and copied to `v(:,:,1)` before the
        gradient's two-electron path.  For DFT/MRSF gradients the same
        transition-density information must be supplied to
        `utddft_xc_gradient` via optional `xa`/`xb`; otherwise that routine
        sets `doFxc=.false.` and skips the derivative XC (`grad_X`) term even
        though conventional TDDFT includes it.
        """
        source = MRSF_GRADIENT.read_text()
        self.assertRegex(source, r"call\s+tagarray_get_data\(infos%dat,\s*OQP_td_abxc,\s*td_abxc\)")
        self.assertRegex(source, r"v\(:,:,1\)\s*=\s*td_abxc")
        block = _call_block(source, "utddft_xc_gradient")
        self.assertRegex(block, r"xa\s*=\s*v\(:,:,1:1\)")
        self.assertRegex(block, r"xb\s*=\s*v\(:,:,2:2\)")


if __name__ == "__main__":
    unittest.main()
