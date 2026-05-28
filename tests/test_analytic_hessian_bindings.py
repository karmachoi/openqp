import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


class TestAnalyticHessianBindings(unittest.TestCase):
    def test_c_header_declares_hf_dft_hessian_entry_point_only(self):
        header = read("include/oqp.h")

        self.assertIn("void hf_hessian(struct oqp_handle_t *inf);", header)
        for symbol in ["tdhf_hessian", "tdhf_sf_hessian", "tdhf_mrsf_hessian"]:
            self.assertNotIn(f"void {symbol}(struct oqp_handle_t *inf);", header)

    def test_python_dispatch_mentions_hf_dft_native_hessian_entry_point_only(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("oqp.hf_hessian", source)
        self.assertNotIn("oqp.tdhf_hessian", source)
        self.assertNotIn("oqp.tdhf_sf_hessian", source)

    def test_hf_native_dispatch_is_explicitly_not_a_numerical_fallback(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("native_hess_func", source)
        self.assertIn("no numerical fallback", source.lower())
        self.assertIn("Native OpenQP HF/DFT analytic Hessian kernel is not available", source)
        self.assertIn("external PySCF is not used for production analytic Hessians", source)
        self.assertIn("native_hess(self.mol)", source)

    def test_response_dispatch_is_deferred_out_of_hf_dft_branch(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("def analytical_sf_hess", source)
        self.assertIn("td_type == 'sf'", source)
        self.assertIn("return self.analytical_sf_hess()", source)
        self.assertIn("Response-theory", source)
        self.assertIn("response-Hessian branch", source)

    def test_hf_hessian_fortran_scaffold_exports_c_abi_without_claiming_support(self):
        source = read("source/modules/hf_hessian.F90")

        self.assertIn("module hf_hessian_mod", source)
        self.assertIn('bind(C, name="hf_hessian")', source)
        self.assertIn("subroutine hf_hessian_C", source)
        self.assertIn("subroutine hf_hessian", source)
        self.assertIn("Analytic HF/DFT Hessian kernel scaffold reached", source)
        self.assertIn("WITH_ABORT", source)

    def test_response_hessian_fortran_scaffolds_are_not_in_hf_dft_branch(self):
        for filename in ["tdhf_hessian.F90", "tdhf_sf_hessian.F90", "tdhf_mrsf_hessian.F90"]:
            with self.subTest(filename=filename):
                self.assertFalse((ROOT / "source/modules" / filename).exists())

    def test_molecule_has_single_hessian_storage_helper_with_asymmetry_metadata(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("def set_hessian_result", source)
        self.assertIn("hessian_metadata", source)
        self.assertIn("max_asymmetry", source)
        self.assertIn("0.5 * (hessian + hessian.T)", source)
        self.assertIn("def get_hess(self):", source)
        self.assertNotIn("def get_hess(self):\n        \"\"\"\n        Get hessian results\n        \"\"\"\n\n        return []", source)

    def test_saved_hessian_json_uses_inertia_not_modes_for_inertia_field(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("'inertia': self.inertia.tolist()", source)
        self.assertNotIn("'inertia': self.modes.tolist()", source)

    def test_read_hessian_json_restores_hessian_metadata(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("self.hessian_metadata = data.get('hessian_metadata', {})", source)


if __name__ == "__main__":
    unittest.main()
