import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


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


if __name__ == "__main__":
    unittest.main()
