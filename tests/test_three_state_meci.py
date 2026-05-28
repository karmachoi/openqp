import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    setattr(mpi_utils, "MPIManager", type("MPIManager", (), {}))
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


def install_libscipy_stubs():
    install_minimal_oqp_stubs()
    oqp_library = sys.modules.setdefault("oqp.library", types.ModuleType("oqp.library"))
    oqp_library.__path__ = []

    single_point = types.ModuleType("oqp.library.single_point")
    for class_name in ("SinglePoint", "Gradient", "LastStep"):
        setattr(single_point, class_name, type(class_name, (), {"__init__": lambda self, *args, **kwargs: None}))
    sys.modules["oqp.library.single_point"] = single_point

    file_utils = types.ModuleType("oqp.utils.file_utils")
    setattr(file_utils, "dump_log", lambda *args, **kwargs: None)
    setattr(file_utils, "dump_data", lambda *args, **kwargs: None)
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    sys.modules["oqp.utils.file_utils"] = file_utils


class TestThreeStateMECI(unittest.TestCase):
    def test_input_checker_accepts_three_state_meci_search_and_kstate_order(self):
        install_minimal_oqp_stubs()
        input_checker = load_module(
            "input_checker_three_state_meci_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "meci", "method": "tdhf"},
            "optimize": {"lib": "scipy", "meci_search": "three_state", "istate": 0, "jstate": 1, "kstate": 2},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_rejects_three_state_meci_without_kstate_above_jstate(self):
        install_minimal_oqp_stubs()
        input_checker = load_module(
            "input_checker_three_state_meci_bad_kstate_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "meci", "method": "tdhf"},
            "optimize": {"lib": "scipy", "meci_search": "three_state", "istate": 0, "jstate": 2, "kstate": 2},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertFalse(report.ok)
        self.assertIn("optimize.kstate", report.to_text())

    def test_input_checker_rejects_three_state_meci_when_kstate_exceeds_nstate(self):
        install_minimal_oqp_stubs()
        input_checker = load_module(
            "input_checker_three_state_meci_nstate_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "meci", "method": "tdhf"},
            "tdhf": {"nstate": 2},
            "optimize": {"lib": "scipy", "meci_search": "three_state", "istate": 0, "jstate": 1, "kstate": 3},
        }

        report = input_checker.CheckReport()
        input_checker._check_requested_states(config, report)

        self.assertFalse(report.ok)
        self.assertIn("tdhf.nstate", report.to_text())
        self.assertIn(">= 3", report.to_text())

    def test_input_checker_rejects_three_state_meci_with_dlfind(self):
        install_minimal_oqp_stubs()
        input_checker = load_module(
            "input_checker_three_state_meci_dlfind_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "meci", "method": "tdhf"},
            "optimize": {"lib": "dlfind", "meci_search": "three_state", "istate": 0, "jstate": 1, "kstate": 2},
            "dlfind": {"ims": 1},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertFalse(report.ok)
        self.assertIn("optimize.lib", report.to_text())

    def test_three_state_penalty_uses_nonredundant_adjacent_gaps(self):
        self.assertIn("'kstate': self.kstate", (ROOT / "pyoqp/oqp/library/libscipy.py").read_text())
        install_libscipy_stubs()
        libscipy = load_module("libscipy_three_state_formula_under_test", "pyoqp/oqp/library/libscipy.py")
        opt = libscipy.MECIOpt.__new__(libscipy.MECIOpt)
        opt.istate = 0
        opt.jstate = 1
        opt.kstate = 2
        opt.sigma = 2.0
        opt.alpha = 0.1
        opt.incre = 0.0
        opt.weights = 1.0
        opt.itr = 1
        opt.pre_energy = 0.0
        opt.pre_coord = np.zeros(2)
        opt.atoms = np.array([[1]])
        opt.mol = types.SimpleNamespace(log_path=".")
        opt.metrics = {
            "itr": 0,
            "nstate": 3,
            "istate": 0,
            "jstate": 1,
            "kstate": 2,
            "meci_search": "three_state",
            "sigma": opt.sigma,
            "alpha": opt.alpha,
            "incre": opt.incre,
            "energy_shift": 1e-6,
            "energy_gap": 1e-4,
            "target_rmsd_step": 1e-3,
            "target_max_step": 2e-3,
            "target_rmsd_grad": 1e-4,
            "target_max_grad": 3e-4,
        }

        coordinates = np.zeros(2)
        energies = np.array([1.0, 1.2, 1.5])
        grads = np.array([
            [[1.0, 0.0]],
            [[0.0, 2.0]],
            [[3.0, 4.0]],
        ])

        f, df = opt.three_state(coordinates, energies, grads)

        gap_ij = 0.2
        gap_jk = 0.3
        expected_f = (1.0 + 1.2 + 1.5) / 3.0 + 2.0 * (
            gap_ij**2 / (gap_ij + 0.1) + gap_jk**2 / (gap_jk + 0.1)
        )
        coef_ij = (gap_ij**2 + 2 * 0.1 * gap_ij) / (gap_ij + 0.1) ** 2
        coef_jk = (gap_jk**2 + 2 * 0.1 * gap_jk) / (gap_jk + 0.1) ** 2
        avg_grad = np.array([4.0 / 3.0, 2.0])
        expected_df = avg_grad + 2.0 * (coef_ij * np.array([-1.0, 2.0]) + coef_jk * np.array([3.0, 2.0]))

        self.assertAlmostEqual(f, expected_f)
        np.testing.assert_allclose(df, expected_df)
        self.assertAlmostEqual(opt.metrics["gap"], 0.5)

    def test_three_state_adaptive_convergence_jumps_sigma_when_gap_is_not_tight(self):
        install_libscipy_stubs()
        libscipy = load_module("libscipy_three_state_convergence_under_test", "pyoqp/oqp/library/libscipy.py")
        opt = libscipy.MECIOpt.__new__(libscipy.MECIOpt)
        opt.mol = types.SimpleNamespace()
        opt.meci_search = "three_state"
        opt.itr = 4
        opt.maxit = 20
        opt.energy_shift = 1e-6
        opt.energy_gap = 1e-4
        opt.rmsd_step = 1e-3
        opt.max_step = 2e-3
        opt.rmsd_grad = 1e-4
        opt.max_grad = 3e-4
        opt.sigma = 2.0
        opt.pen_jump = 10.0
        opt.metrics = {
            "de": 0.0,
            "gap": 5e-4,
            "rmsd_step": 0.0,
            "max_step": 0.0,
            "rmsd_grad": 0.0,
            "max_grad": 0.0,
        }

        opt.check_convergence()

        self.assertEqual(opt.sigma, 12.0)
        self.assertEqual(opt.metrics["sigma"], 12.0)


if __name__ == "__main__":
    unittest.main()
