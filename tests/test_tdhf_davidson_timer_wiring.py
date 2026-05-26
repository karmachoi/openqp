import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UTIL = ROOT / "source" / "util.F90"
MRSF_ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"


class TDHFDavidsonTimerWiringTests(unittest.TestCase):
    def test_util_provides_oqp_timer_log_line_helper(self):
        text = UTIL.read_text()

        self.assertIn("subroutine log_oqp_timer", text)
        self.assertIn('"OQP_TIMER label="', text)
        self.assertIn('" seconds="', text)

    def test_mrsf_davidson_loop_emits_total_sigma_and_metc_timer_labels(self):
        text = MRSF_ENERGY.read_text()

        self.assertRegex(text, r"use\s+util,\s+only:.*log_oqp_timer")
        for label in (
            "tdhf.response.total",
            "tdhf.davidson.total",
            "tdhf.davidson.sigma_build",
            "tdhf.davidson.metc_contract",
        ):
            self.assertIn(label, text)

        loop_body = re.search(r"do iter = 1, mxiter(?P<body>.*?)call rparedms", text, re.S)
        assert loop_body is not None
        body = loop_body.group("body")
        self.assertIn("timer_sigma_start", body)
        self.assertIn("timer_metc_start", body)
        self.assertIn("call log_oqp_timer(iw, \"tdhf.davidson.sigma_build\"", body)
        self.assertIn("call log_oqp_timer(iw, \"tdhf.davidson.metc_contract\"", body)


if __name__ == "__main__":
    unittest.main()
