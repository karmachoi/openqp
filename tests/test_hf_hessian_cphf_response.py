import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


class TestHfCphfHessianResponseScaffold(unittest.TestCase):
    def test_hf_cphf_response_module_exports_guarded_solver_boundary(self):
        source = read("source/modules/hf_cphf_response.F90")

        self.assertIn("module hf_cphf_response_mod", source)
        self.assertIn("public :: hf_cphf_hessian_response", source)
        self.assertIn("subroutine hf_cphf_hessian_response", source)
        self.assertIn("hf_cphf_hessian_response_scaffold", source)
        self.assertIn("partial_kernel", source)
        self.assertIn("WITH_ABORT", source)

    def test_hf_cphf_response_scaffold_documents_required_orbital_response_invariants(self):
        source = read("source/modules/hf_cphf_response.F90")

        for token in [
            "occupied-virtual",
            "orbital-energy denominator",
            "Fock derivative RHS",
            "two-electron response coupling",
            "finite-difference gradient validation",
        ]:
            self.assertIn(token, source)

        self.assertRegex(source, r"nocc\s*=\s*infos%mol_prop%nocc")
        self.assertRegex(source, r"nvir\s*=\s*nbasis\s*-\s*nocc")
        self.assertIn("if (nocc <= 0 .or. nvir <= 0)", source)

    def test_hf_hessian_keeps_cphf_boundary_deferred(self):
        source = read("source/modules/hf_hessian.F90")

        self.assertIn("hf_cphf_hessian_response", source)
        self.assertIn("CPHF orbital-response block deferred", source)
        self.assertNotIn("call hf_cphf_hessian_response", source)


if __name__ == "__main__":
    unittest.main()
