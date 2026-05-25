import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        pass

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class TestGpuXcResponseConfig(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_runtime_helper_recognizes_xc_response_target(self):
        gpu = load_module("gpu_xc_response_under_test", "pyoqp/oqp/utils/gpu.py")

        config = gpu.GpuConfig.from_config(
            {
                "gpu": {
                    "enabled": True,
                    "backend": "cuda",
                    "target": "xc_response",
                    "device": "1",
                    "precision": "float64",
                    "fallback": "cpu",
                }
            }
        )

        self.assertTrue(config.enabled)
        self.assertEqual(config.target, "xc_response")
        self.assertEqual(config.device, 1)
        self.assertTrue(config.targets_xc_response)
        self.assertFalse(config.targets_metc)

    def test_input_checker_accepts_tdhf_dft_xc_response_gpu_target(self):
        input_checker = load_module(
            "input_checker_gpu_xc_response_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"method": "tdhf", "functional": "bhhlyp"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "multiplicity": 1, "nstate": 2, "nvdav": 4},
            "gpu": {"enabled": True, "backend": "cuda", "target": "xc_response"},
        }

        report = input_checker.CheckReport()
        input_checker._check_gpu(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_rejects_xc_response_without_dft_functional(self):
        input_checker = load_module(
            "input_checker_gpu_xc_response_reject_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "multiplicity": 1, "nstate": 2, "nvdav": 4},
            "gpu": {"enabled": True, "backend": "cuda", "target": "xc_response"},
        }

        report = input_checker.CheckReport()
        input_checker._check_gpu(config, report)

        self.assertFalse(report.ok)
        self.assertIn("input.functional", report.to_text())


if __name__ == "__main__":
    unittest.main()
