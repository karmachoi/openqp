import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_cache_module():
    path = ROOT / "pyoqp/oqp/utils/tdhf_xc_response_cache.py"
    spec = importlib.util.spec_from_file_location("tdhf_xc_response_cache_under_test", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TdhfXcResponseCachePlanTest(unittest.TestCase):
    def test_cache_key_distinguishes_response_type_and_spin_channels(self):
        XcResponseCachePlan = load_cache_module().XcResponseCachePlan

        rpa = XcResponseCachePlan(
            nbf=24,
            ngrid=1024,
            functional="bhhlyp",
            basis="6-31g*",
            scf_type="rhf",
            response_type="rpa",
            spin_channels=1,
        )
        tda = XcResponseCachePlan(
            nbf=24,
            ngrid=1024,
            functional="bhhlyp",
            basis="6-31g*",
            scf_type="rhf",
            response_type="tda",
            spin_channels=2,
        )

        self.assertEqual(rpa.reuse_key(), ("bhhlyp", "6-31g*", "rhf", "rpa", 24, 1024, 1))
        self.assertNotEqual(rpa.reuse_key(), tda.reuse_key())

    def test_workspace_sizes_cover_density_potential_and_weights(self):
        XcResponseCachePlan = load_cache_module().XcResponseCachePlan

        plan = XcResponseCachePlan(
            nbf=10,
            ngrid=50,
            functional="b3lypv5",
            basis="3-21g",
            scf_type="rohf",
            response_type="rpa",
            spin_channels=2,
        )

        self.assertEqual(plan.density_values, 100)
        self.assertEqual(plan.potential_values, 100)
        self.assertEqual(plan.weight_values, 50)
        self.assertEqual(plan.ao_grid_values, 500)
        self.assertEqual(plan.total_scalar_values(), 750)

    def test_workspace_byte_count_uses_explicit_dtype_width(self):
        XcResponseCachePlan = load_cache_module().XcResponseCachePlan

        plan = XcResponseCachePlan(
            nbf=10,
            ngrid=50,
            functional="b3lypv5",
            basis="3-21g",
            scf_type="rohf",
            response_type="rpa",
            spin_channels=2,
        )

        self.assertEqual(plan.total_workspace_bytes(dtype_bytes=8), 6000)
        self.assertEqual(plan.total_workspace_bytes(dtype_bytes=4), 3000)

    def test_rejects_nonpositive_dimensions(self):
        XcResponseCachePlan = load_cache_module().XcResponseCachePlan

        with self.assertRaises(ValueError):
            XcResponseCachePlan(
                nbf=0,
                ngrid=50,
                functional="bhhlyp",
                basis="6-31g*",
                scf_type="rhf",
                response_type="rpa",
                spin_channels=1,
            )


if __name__ == "__main__":
    unittest.main()
