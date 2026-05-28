import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


def subroutine_body(source, name):
    start = source.index(f"SUBROUTINE {name}")
    end = source.index("END SUBROUTINE", start)
    return source[start:end]


class TestOneElectronHessianInfrastructure(unittest.TestCase):
    def test_primitive_module_exports_second_derivative_entry_points(self):
        source = read("source/integrals/mod_1e_primitives.F90")

        for name in [
            "comp_overlap_der2",
            "comp_kinetic_der2",
            "comp_coulomb_der2",
        ]:
            self.assertIn(f"PUBLIC {name}", source)
            self.assertRegex(source, rf"SUBROUTINE\s+{name}\s*\(")

    def test_second_derivative_entry_points_are_guarded_until_validated(self):
        source = read("source/integrals/mod_1e_primitives.F90")

        self.assertIn("one_electron_hessian_der2_scaffold", source)
        self.assertIn("finite_difference_validation_required", source)
        self.assertIn("no production one-electron Hessian support", source)
        self.assertIn("WITH_ABORT", source)

    def test_finite_difference_validation_hooks_reference_existing_der1_routines(self):
        source = read("source/integrals/mod_1e_primitives.F90")

        self.assertRegex(source, r"SUBROUTINE\s+validate_overlap_der2_by_finite_difference\s*\(")
        self.assertRegex(source, r"SUBROUTINE\s+validate_kinetic_der2_by_finite_difference\s*\(")
        self.assertRegex(source, r"SUBROUTINE\s+validate_coulomb_der2_by_finite_difference\s*\(")
        self.assertIn("comp_overlap_der1", source)
        self.assertIn("comp_kinetic_der1", source)
        self.assertIn("comp_coulomb_der1", source)
        self.assertIn("central finite difference", source.lower())
        self.assertIn("fd_step", source)
    def test_overlap_der2_has_native_analytic_assembly_without_abort(self):
        source = read("source/integrals/mod_1e_primitives.F90")
        body = subroutine_body(source, "comp_overlap_der2")

        self.assertIn("overlap_xyz", body)
        self.assertIn("der_kinovl_xyz", body)
        self.assertIn("der2_kinovl_xyz", body)
        self.assertIn("ovl_der2", body)
        self.assertIn("der2(1,1)", body)
        self.assertIn("der2(1,2)", body)
        self.assertIn("der2(2,3)", body)
        self.assertIn("der2(3,3)", body)
        self.assertIn("pp%expfac", body)
        self.assertNotIn("WITH_ABORT", body)
        self.assertNotIn("no production one-electron Hessian support", body)

    def test_overlap_second_derivative_1d_helper_uses_twice_differentiated_recursion(self):
        source = read("source/integrals/mod_1e_primitives.F90")

        self.assertRegex(source, r"SUBROUTINE\s+der2_kinovl_xyz\s*\(")
        self.assertIn("4.0_REAL64 * ai * ai", source)
        self.assertIn("2.0_REAL64 * ai * (2 * i + 1)", source)
        self.assertIn("i * (i - 1)", source)


if __name__ == "__main__":
    unittest.main()
