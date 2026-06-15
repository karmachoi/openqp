import importlib.util
import math
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_input_checker():
    oqp_mod = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    utils_mod = sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    setattr(oqp_mod, "utils", utils_mod)
    setattr(utils_mod, "mpi_utils", mpi_utils)

    return load_module("openqp_mrsf_reference_checker_under_test", "pyoqp/oqp/utils/input_checker.py")


class TestMrsfReferenceParser(unittest.TestCase):
    def setUp(self):
        self.mrsf_reference = load_module(
            "mrsf_reference_under_test",
            "pyoqp/oqp/utils/mrsf_reference.py",
        )

    def test_parse_candidate_open_shell_pairs(self):
        pairs = self.mrsf_reference.parse_reference_pairs("12:13; 11,14")

        self.assertEqual(pairs, [(12, 13), (11, 14)])

    def test_explicit_weights_must_sum_to_one(self):
        with self.assertRaises(self.mrsf_reference.MrsfReferenceError):
            self.mrsf_reference.parse_weights("0.7,0.7", 2)

    def test_auto_weights_are_equal(self):
        self.assertEqual(self.mrsf_reference.parse_weights("equal", 2), [0.5, 0.5])

    def test_metadata_infers_current_rohf_open_pair_and_gap_risk(self):
        config = {
            "mrsf_ref": {
                "mode": "diagnostic",
                "open_pairs": "auto",
                "weights": "equal",
                "max_refs": 2,
                "gap_threshold": 0.01,
                "overlap_threshold": 0.85,
                "strict": False,
            }
        }
        data = {
            "nelec_A": 4,
            "nelec_B": 2,
            "OQP::E_MO_A": [-1.0, -0.500, -0.495, -0.100],
            "OQP::E_MO_B": [-1.0, -0.500, -0.495, -0.100],
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, data)

        self.assertEqual(metadata["status"], "diagnostic")
        self.assertEqual(metadata["open_pairs"], [[3, 4]])
        self.assertEqual(metadata["weights"], [1.0])
        self.assertTrue(metadata["frontier"]["ambiguous"])
        self.assertIn("closed_to_open", metadata["frontier"]["gaps_hartree"]["alpha"])

    def test_ensemble_marks_overlap_offdiagonal_response_coupling(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            }
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, {})

        self.assertTrue(metadata["implemented"])
        self.assertTrue(metadata["response_implemented"])
        self.assertEqual(metadata["open_pairs"], [[3, 4], [2, 5]])
        self.assertTrue(metadata["theory"]["inter_reference_coupling"])
        self.assertEqual(metadata["response_coupling_model"], "overlap_offdiagonal")
        self.assertEqual(metadata["theory"]["inter_reference_coupling_model"], "overlap_offdiagonal")
        self.assertEqual(metadata["theory"]["target_response_model"], "native sigma-action off-diagonal MRSF response over all reference-specific spin-flip spaces")
        self.assertTrue(any("state-interaction" in item for item in metadata["warnings"]))
        self.assertEqual(metadata["theory"]["reference_model"], "mixed_rohf_triplet_reference_ensemble")
        self.assertEqual(metadata["trial_vector_model"]["mode"], "adaptive")
        self.assertEqual(metadata["trial_vector_model"]["active_virtual_shift_hartree"], 1.0e6)

    def test_ensemble_can_request_block_diagonal_response_for_comparison(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
                "coupling": "block_diagonal",
            }
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, {})

        self.assertFalse(metadata["theory"]["inter_reference_coupling"])
        self.assertEqual(metadata["response_coupling_model"], "block_diagonal")
        self.assertTrue(any("block_diagonal" in item for item in metadata["warnings"]))

    def test_trial_vector_policy_accepts_native_mode(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
                "trial_vectors": "native",
                "trial_shift": 100.0,
            }
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, {})

        self.assertEqual(metadata["trial_vector_model"]["mode"], "native")
        self.assertEqual(metadata["trial_vector_model"]["active_virtual_shift_hartree"], 100.0)

    def test_state_average_alias_canonicalizes_to_ensemble(self):
        config = {
            "mrsf_ref": {
                "mode": "state_average",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            }
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, {})

        self.assertEqual(metadata["mode"], "ensemble")
        self.assertEqual(metadata["status"], "ensemble_requested")

    def test_ensemble_auto_selects_frontier_reference_pairs(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "auto",
                "weights": "equal",
                "max_refs": 2,
            }
        }
        data = {
            "nelec_A": 6,
            "nelec_B": 4,
            "nbf": 13,
            "OQP::E_MO_A": [
                -20.0, -1.4, -0.9, -0.445, -0.440, -0.400, -0.395,
                0.10, 0.20, 0.35, 0.50, 0.70, 0.90,
            ],
            "OQP::E_MO_B": [
                -20.0, -1.4, -0.9, -0.445, -0.440, -0.400, -0.395,
                0.10, 0.20, 0.35, 0.50, 0.70, 0.90,
            ],
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, data)
        alpha_occ, beta_occ = self.mrsf_reference.ensemble_occupation_vectors(metadata)

        self.assertEqual(metadata["pair_selection"]["mode"], "auto")
        self.assertEqual(metadata["pair_selection"]["strategy"], "frontier_window")
        self.assertEqual(metadata["open_pairs"], [[5, 6], [4, 7]])
        self.assertEqual(metadata["weights"], [0.5, 0.5])
        self.assertEqual(alpha_occ[:7], [1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5])
        self.assertEqual(beta_occ[:7], [1.0, 1.0, 1.0, 0.5, 0.5, 0.0, 0.0])

    def test_ensemble_auto_excludes_promoted_high_high_pair(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "auto",
                "weights": "equal",
            }
        }
        data = {
            "nelec_A": 6,
            "nelec_B": 4,
            "nbf": 13,
            "OQP::E_MO_A": [
                -20.0, -1.4, -0.9, -0.445, -0.440, -0.400, -0.395,
                0.10, 0.20, 0.35, 0.50, 0.70, 0.90,
            ],
            "OQP::E_MO_B": [
                -20.0, -1.4, -0.9, -0.445, -0.440, -0.400, -0.395,
                0.10, 0.20, 0.35, 0.50, 0.70, 0.90,
            ],
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, data)

        self.assertEqual(metadata["max_refs"], 6)
        self.assertEqual(metadata["pair_selection"]["active_orbitals"], [4, 5, 6, 7])
        self.assertEqual(metadata["open_pairs"], [[5, 6], [4, 7], [4, 6], [5, 7], [4, 5]])
        self.assertEqual(metadata["weights"], [1.0 / 5.0] * 5)
        self.assertEqual(metadata["pair_selection"]["excluded_pairs"][0]["pair"], [6, 7])
        self.assertEqual(
            metadata["pair_selection"]["excluded_pairs"][0]["excluded_reason"],
            "promoted_high_high_pair",
        )
        self.assertFalse(metadata["pair_selection"]["truncated"])

    def test_ensemble_auto_uses_gap_threshold_for_active_window(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "auto",
                "weights": "equal",
                "max_refs": 6,
                "gap_threshold": 0.01,
            }
        }
        data = {
            "nelec_A": 9,
            "nelec_B": 7,
            "nbf": 12,
            "OQP::E_MO_A": [
                -20.0, -10.0, -1.0, -0.8, -0.6, -0.4, -0.30,
                -0.20, -0.17, -0.08, 0.2, 0.4,
            ],
            "OQP::E_MO_B": [
                -20.0, -10.0, -1.0, -0.8, -0.6, -0.4, -0.30,
                -0.20, -0.17, -0.08, 0.2, 0.4,
            ],
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, data)

        self.assertEqual(metadata["pair_selection"]["active_orbitals"], [8, 9])
        self.assertEqual(metadata["open_pairs"], [[8, 9]])
        self.assertEqual(metadata["weights"], [1.0])

    def test_gap_softmax_weights_follow_reference_energy_proxy(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "gap_softmax",
                "weight_temperature": 0.7,
            }
        }
        data = {
            "nelec_A": 4,
            "nelec_B": 2,
            "OQP::E_MO_A": [-1.0, -0.8, -0.3, -0.2, 0.0, 0.4],
            "OQP::E_MO_B": [-1.0, -0.8, -0.3, -0.2, 0.0, 0.4],
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, data)
        alpha_occ, beta_occ = self.mrsf_reference.ensemble_occupation_vectors(metadata)

        expected_first = 1.0 / (1.0 + math.exp(-1.0))
        expected_second = 1.0 - expected_first
        self.assertEqual(metadata["weight_model"]["mode"], "gap_softmax")
        self.assertTrue(metadata["weight_model"]["resolved"])
        self.assertAlmostEqual(metadata["weights"][0], expected_first)
        self.assertAlmostEqual(metadata["weights"][1], expected_second)
        self.assertGreater(metadata["weights"][0], metadata["weights"][1])
        self.assertAlmostEqual(alpha_occ[3], expected_first)
        self.assertAlmostEqual(alpha_occ[4], expected_second)
        self.assertAlmostEqual(beta_occ[1], expected_first)
        self.assertAlmostEqual(beta_occ[2], expected_second)

    def test_reference_ensemble_builds_fractional_occupations_and_response_spaces(self):
        data = {
            "nelec_A": 4,
            "nelec_B": 2,
            "OQP::E_MO_A": [-1.0, -0.8, -0.2, -0.2, 0.1, 0.3],
            "OQP::E_MO_B": [-1.0, -0.8, -0.2, -0.2, 0.1, 0.3],
        }

        ensemble = self.mrsf_reference.build_reference_ensemble(
            [(3, 4), (2, 5)],
            [0.25, 0.75],
            data,
        )

        self.assertTrue(ensemble["available"])
        self.assertEqual(ensemble["active_open_orbitals"], [2, 3, 4, 5])
        self.assertEqual(ensemble["references"][0]["closed_orbitals"], [1, 2])
        self.assertEqual(ensemble["references"][1]["closed_orbitals"], [1, 3])
        self.assertEqual(ensemble["references"][1]["open_pair"], [2, 5])
        self.assertEqual(ensemble["references"][1]["response_space"]["raw_dimension"], 16)
        self.assertEqual(ensemble["references"][1]["response_space"]["triplet_dimension"], 13)
        self.assertEqual(ensemble["response_space"]["raw_dimension"], 32)
        self.assertEqual(ensemble["response_space"]["triplet_dimension"], 26)

        alpha = {item["mo"]: item["occupation"] for item in ensemble["ensemble_occupations"]["alpha"]}
        beta = {item["mo"]: item["occupation"] for item in ensemble["ensemble_occupations"]["beta"]}
        self.assertEqual(alpha, {1: 1.0, 2: 1.0, 3: 1.0, 4: 0.25, 5: 0.75})
        self.assertEqual(beta, {1: 1.0, 2: 0.25, 3: 0.75})
        self.assertAlmostEqual(ensemble["ensemble_occupations"]["alpha_sum"], 4.0)
        self.assertAlmostEqual(ensemble["ensemble_occupations"]["beta_sum"], 2.0)

    def test_reference_mo_permutation_places_closed_then_open_pair(self):
        reference = {
            "closed_orbitals": [1, 2, 3, 7],
            "open_pair": [8, 9],
        }

        permutation = self.mrsf_reference.reference_mo_permutation(reference, 10)

        self.assertEqual(permutation, [1, 2, 3, 7, 8, 9, 4, 5, 6, 10])

    def test_collect_block_diagonal_response_sorts_reference_blocks(self):
        blocks = [
            {
                "reference_id": 1,
                "weight": 1.0 / 6.0,
                "open_pair": [8, 9],
                "converged": True,
                "energies": [0.42, 0.35],
            },
            {
                "reference_id": 2,
                "weight": 1.0 / 6.0,
                "open_pair": [7, 10],
                "converged": True,
                "energies": [0.31, 0.50],
            },
        ]

        combined = self.mrsf_reference.collect_block_diagonal_response(blocks, 3)

        self.assertEqual(combined["model"], "block_diagonal_uncoupled")
        self.assertEqual(combined["energies"], [0.31, 0.35, 0.42])
        self.assertEqual(
            [(item["reference_id"], item["state_index"]) for item in combined["selected_states"]],
            [(2, 1), (1, 2), (1, 1)],
        )

    def test_collect_block_diagonal_response_ignores_unconverged_blocks(self):
        blocks = [
            {
                "reference_id": 1,
                "weight": 0.5,
                "open_pair": [8, 9],
                "converged": True,
                "energies": [0.42],
            },
            {
                "reference_id": 2,
                "weight": 0.5,
                "open_pair": [7, 10],
                "converged": False,
                "energies": [-2.0],
            },
        ]

        combined = self.mrsf_reference.collect_block_diagonal_response(blocks, 1)

        self.assertEqual(combined["energies"], [0.42])
        self.assertEqual(combined["candidate_count"], 1)
        self.assertEqual(combined["raw_candidate_count"], 2)
        self.assertEqual(combined["skipped_nonconverged_blocks"], [2])

    def test_mrsf_response_labels_map_local_vector_to_original_mos(self):
        reference = {
            "closed_orbitals": [1],
            "open_pair": [2, 4],
        }

        labels = self.mrsf_reference.mrsf_response_labels(reference, 4)

        self.assertEqual(labels[:4], [(1, 2), (2, 2), (4, 2), (1, 4)])
        self.assertEqual(len(labels), 9)

    def test_state_interaction_response_uses_common_basis_offdiagonal_overlap(self):
        references = [
            {"id": 1, "closed_orbitals": [1], "open_pair": [2, 3]},
            {"id": 2, "closed_orbitals": [1], "open_pair": [2, 4]},
        ]
        blocks = [
            {
                "reference_id": 1,
                "weight": 0.5,
                "open_pair": [2, 3],
                "converged": True,
                "energies": [0.30],
            },
            {
                "reference_id": 2,
                "weight": 0.5,
                "open_pair": [2, 4],
                "converged": True,
                "energies": [0.40],
            },
        ]
        vec1 = [0.0] * 9
        vec2 = [0.0] * 9
        vec1[0] = 1.0
        vec2[0] = 0.6
        vec2[2] = 0.8
        block_vectors = {
            1: [[item] for item in vec1],
            2: [[item] for item in vec2],
        }

        mixed = self.mrsf_reference.collect_state_interaction_response(
            blocks,
            block_vectors,
            references,
            nstate=1,
            nmo=4,
        )

        self.assertEqual(mixed["status"], "ready")
        self.assertEqual(mixed["model"], "state_interaction_offdiagonal")
        self.assertEqual(mixed["coupling"], "overlap_offdiagonal")
        self.assertEqual(mixed["candidate_count"], 2)
        self.assertEqual(mixed["common_dimension"], 2)
        self.assertAlmostEqual(mixed["overlap_matrix"][0][1], 0.6)
        self.assertAlmostEqual(mixed["hamiltonian_matrix"][0][1], 0.21)
        self.assertEqual(mixed["offdiagonal_count"], 1)
        self.assertAlmostEqual(mixed["max_abs_offdiagonal_hamiltonian"], 0.21)
        self.assertEqual(len(mixed["selected_states"][0]["components"]), 2)

    def test_ensemble_occupation_vectors_are_dense_spin_occupations(self):
        config = {
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.25,0.75",
            }
        }
        data = {
            "nelec_A": 4,
            "nelec_B": 2,
            "OQP::E_MO_A": [-1.0, -0.8, -0.2, -0.2, 0.1, 0.3],
            "OQP::E_MO_B": [-1.0, -0.8, -0.2, -0.2, 0.1, 0.3],
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, data)
        alpha_occ, beta_occ = self.mrsf_reference.ensemble_occupation_vectors(metadata)

        self.assertEqual(alpha_occ, [1.0, 1.0, 1.0, 0.25, 0.75, 0.0])
        self.assertEqual(beta_occ, [1.0, 0.25, 0.75, 0.0, 0.0, 0.0])

    def test_reference_ensemble_rejects_out_of_range_pair_when_mo_count_is_known(self):
        data = {
            "nelec_A": 4,
            "nelec_B": 2,
            "OQP::E_MO_A": [-1.0, -0.8, -0.2, -0.2],
        }

        ensemble = self.mrsf_reference.build_reference_ensemble([(3, 6)], [1.0], data)

        self.assertFalse(ensemble["available"])
        self.assertEqual(ensemble["status"], "invalid_reference")
        self.assertIn("exceeds available MO count", ensemble["reason"])


class TestMrsfReferenceInputChecker(unittest.TestCase):
    def setUp(self):
        self.input_checker = load_input_checker()

    def test_diagnostic_mode_accepts_rohf_triplet_mrsf(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "diagnostic",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
                "max_refs": 2,
                "gap_threshold": 0.01,
                "overlap_threshold": 0.85,
                "strict": False,
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_rejects_reference_ensemble_outside_mrsf(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "hf"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa"},
            "mrsf_ref": {"mode": "diagnostic"},
        }

        self.input_checker._check_mrsf_ref(config, report)

        errors = "\n".join(item.path for item in report.errors)
        self.assertIn("mrsf_ref.mode", errors)

    def test_ensemble_is_warned_but_kept_parseable(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertTrue(report.ok, report.to_text())
        self.assertIn("mrsf_ref.mode", "\n".join(item.path for item in report.warnings))

    def test_ensemble_accepts_auto_references(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "auto",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertTrue(report.ok, report.to_text())
        self.assertIn("mrsf_ref.mode", "\n".join(item.path for item in report.warnings))

    def test_ensemble_rejects_single_manual_reference(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("mrsf_ref.open_pairs", "\n".join(item.path for item in report.errors))

    def test_ensemble_rejects_pfon_smearing(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3, "pfon": True},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("scf.pfon", "\n".join(item.path for item in report.errors))

    def test_ensemble_rejects_non_energy_runtype(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf", "runtype": "grad"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("input.runtype", "\n".join(item.path for item in report.errors))

    def test_bad_weights_are_errors(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "diagnostic",
                "open_pairs": "3:4;2:5",
                "weights": "0.8,0.8",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("weights", report.errors[0].message)

    def test_bad_weight_temperature_is_an_error(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "ensemble",
                "open_pairs": "auto",
                "weights": "gap_softmax",
                "weight_temperature": 0,
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("weight_temperature", report.errors[0].message)


if __name__ == "__main__":
    unittest.main()
