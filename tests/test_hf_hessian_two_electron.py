import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


class TestTwoElectronHessianInfrastructure(unittest.TestCase):
    def test_two_electron_hessian_module_exports_guarded_driver(self):
        source = read("source/integrals/grd2_hessian.F90")

        self.assertIn("module grd2_hessian", source)
        self.assertIn("public :: grd2_hessian_driver", source)
        self.assertIn("subroutine grd2_hessian_driver", source)
        self.assertIn("two_electron_hessian_der2_scaffold", source)
        self.assertIn("partial_kernel", source)
        self.assertIn("WITH_ABORT", source)
        self.assertIn("grd2_rys", source)
        self.assertIn("grd2_int_data_t", source)

    def test_two_electron_hessian_scaffold_allocates_second_derivative_integral_workspace(self):
        source = read("source/integrals/grd2_hessian.F90")

        self.assertRegex(source, r"call\s+gdat%init\s*\([^\n]*maxang[^\n]*,\s*2\s*,")
        self.assertIn("two-electron second derivatives", source.lower())
        self.assertIn("finite-difference", source.lower())
        self.assertIn("grd2_driver", source)

    def test_hf_hessian_mentions_deferred_two_electron_block(self):
        source = read("source/modules/hf_hessian.F90")

        self.assertIn("two-electron Hessian block deferred", source)
        self.assertIn("grd2_hessian_driver", source)
        self.assertNotIn("call grd2_hessian_driver", source)


if __name__ == "__main__":
    unittest.main()
