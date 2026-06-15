import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TestMrsfReferenceScfWiring(unittest.TestCase):
    def test_density_builder_accepts_explicit_spin_occupations(self):
        source = (ROOT / "source" / "guess.F90").read_text()

        self.assertIn("alpha_occupation", source)
        self.assertIn("beta_occupation", source)
        self.assertIn(
            "call orb_to_dens(alpha_density, alpha_orbital, alpha_occupation, nbasis",
            source,
        )
        self.assertIn(
            "call orb_to_dens(beta_density, alpha_orbital, beta_occupation, nbasis",
            source,
        )

    def test_native_scf_consumes_mrsf_reference_occupation_tags(self):
        source = (ROOT / "source" / "scf.F90").read_text()

        self.assertIn("OQP_mrsf_ref_occ_a", source)
        self.assertIn("OQP_mrsf_ref_occ_b", source)
        self.assertIn("validate_mrsf_ref_occupation", source)
        self.assertIn("form_rohf_ensemble_fock", source)
        self.assertIn("alpha_occupation=mrsf_ref_occ_a", source)
        self.assertIn("beta_occupation=mrsf_ref_occ_b", source)

    def test_tagarray_defines_mrsf_reference_occupation_tags(self):
        source = (ROOT / "source" / "tagarray_driver.F90").read_text()

        self.assertIn('OQP_mrsf_ref_occ_a = OQP_prefix // "mrsf_ref_occ_a"', source)
        self.assertIn('OQP_mrsf_ref_occ_b = OQP_prefix // "mrsf_ref_occ_b"', source)

    def test_python_stages_occupation_tags_before_scf(self):
        source = (ROOT / "pyoqp" / "oqp" / "library" / "single_point.py").read_text()

        self.assertIn('self.mol.data["OQP::mrsf_ref_occ_a"] = alpha_occ', source)
        self.assertIn('self.mol.data["OQP::mrsf_ref_occ_b"] = beta_occ', source)
        self.assertIn("'applied_weights': metadata.get('weights', [])", source)
        self.assertIn("'applied_weight_model': metadata.get('weight_model', {})", source)
        self.assertLess(
            source.index("self._prepare_mrsf_reference_scf()"),
            source.index("scf_flag = self._run_scf()"),
        )

    def test_python_runs_block_response_before_legacy_guard(self):
        source = (ROOT / "pyoqp" / "oqp" / "library" / "single_point.py").read_text()

        self.assertIn("def _run_mrsf_reference_response(self):", source)
        self.assertIn("reference_mo_permutation(reference, nmo)", source)
        self.assertIn("collect_block_diagonal_response(blocks, requested_nstate)", source)
        self.assertIn("collect_state_interaction_response(", source)
        self.assertIn("'state_interaction': state_interaction", source)
        self.assertIn("block_diagonal = collect_block_diagonal_response", source)
        self.assertIn("final_response.get('model', 'block_diagonal_uncoupled')", source)
        self.assertIn("def _apply_mrsf_reference_trial_vector_policy(self", source)
        self.assertIn("_apply_mrsf_reference_trial_vector_policy(", source)
        self.assertIn("active_virtual_shift_hartree", source)
        self.assertIn("Do not use it with EKT, gradients, NAC, RT-MRSF, or transition-property workflows yet.", source)
        self.assertLess(
            source.index("if self._run_mrsf_reference_response():"),
            source.index("self._guard_mrsf_reference_mode()"),
        )
        self.assertLess(
            source.index("self.runtype != 'energy'"),
            source.index("if self.runtype == 'ekt':"),
        )


if __name__ == "__main__":
    unittest.main()
