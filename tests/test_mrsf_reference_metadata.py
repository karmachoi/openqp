import importlib.util
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

    def test_state_average_marks_unimplemented_coupled_response(self):
        config = {
            "mrsf_ref": {
                "mode": "state_average",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            }
        }

        metadata = self.mrsf_reference.build_mrsf_reference_metadata(config, {})

        self.assertFalse(metadata["implemented"])
        self.assertEqual(metadata["open_pairs"], [[3, 4], [2, 5]])
        self.assertTrue(any("coupled ensemble-response" in item for item in metadata["warnings"]))
        self.assertEqual(metadata["theory"]["reference_model"], "weighted_rohf_triplet_ensemble")

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

    def test_ensemble_occupation_vectors_are_dense_spin_occupations(self):
        config = {
            "mrsf_ref": {
                "mode": "state_average",
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

    def test_state_average_is_warned_but_kept_parseable(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "state_average",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertTrue(report.ok, report.to_text())
        self.assertIn("mrsf_ref.mode", "\n".join(item.path for item in report.warnings))

    def test_state_average_requires_multiple_explicit_references(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "state_average",
                "open_pairs": "auto",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("mrsf_ref.open_pairs", "\n".join(item.path for item in report.errors))

    def test_state_average_rejects_pfon_smearing(self):
        report = self.input_checker.CheckReport()
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rohf", "multiplicity": 3, "pfon": True},
            "tdhf": {"type": "mrsf"},
            "mrsf_ref": {
                "mode": "state_average",
                "open_pairs": "3:4;2:5",
                "weights": "0.5,0.5",
            },
        }

        self.input_checker._check_mrsf_ref(config, report)

        self.assertFalse(report.ok)
        self.assertIn("scf.pfon", "\n".join(item.path for item in report.errors))

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


if __name__ == "__main__":
    unittest.main()
