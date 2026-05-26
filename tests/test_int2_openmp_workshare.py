from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
INT2 = ROOT / "source" / "integrals" / "int2.F90"


def _int2_twoei_source() -> str:
    text = INT2.read_text()
    start = text.index("subroutine int2_twoei")
    end = text.index("end subroutine int2_twoei", start)
    return text[start:end]


class Int2OpenmpWorkshareTest(unittest.TestCase):
    def test_int2_twoei_uses_single_flat_shell_pair_workshare(self):
        src = _int2_twoei_source().lower()
        compact = src.replace(" ", "")
        self.assertIn("ij_pair", src)
        self.assertIn("npairs=nshell*(nshell+1)/2", compact)
        self.assertIn("!$omp do schedule(dynamic,1)", src)
        self.assertRegex(src, r"do\s+ij_pair\s*=\s*1\s*,\s*npairs")

    def test_int2_twoei_avoids_nested_dynamic_k_workshare_and_nowait(self):
        src = _int2_twoei_source().lower()
        self.assertNotIn("!$omp do schedule(dynamic,2)", src)
        self.assertNotIn("!$omp end do nowait", src)


if __name__ == "__main__":
    unittest.main()
