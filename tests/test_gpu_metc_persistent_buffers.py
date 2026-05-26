import importlib.util
import sys
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


class GpuMetcPersistentBufferTests(unittest.TestCase):
    def test_buffer_plan_accounts_for_reusable_metc_device_arrays(self):
        buffers = load_module(
            "gpu_metc_buffers_under_test", "pyoqp/oqp/utils/gpu_metc_buffers.py"
        )

        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=7,
            nf=5,
            nmatrix=3,
            max_integrals=11,
            dtype_bytes=8,
        )

        self.assertEqual(plan.reuse_key, (7, 5, 3, 11, 8))
        self.assertEqual(plan.bytes_for("ids"), 11 * 4 * 4)
        self.assertEqual(plan.bytes_for("integrals"), 11 * 8)
        self.assertEqual(plan.bytes_for("density"), 3 * 5 * 7 * 7 * 8)
        self.assertEqual(plan.bytes_for("fock"), 3 * 5 * 7 * 7 * 8)
        self.assertEqual(
            plan.total_bytes,
            plan.bytes_for("ids")
            + plan.bytes_for("integrals")
            + plan.bytes_for("density")
            + plan.bytes_for("fock"),
        )

    def test_buffer_plan_rejects_nonpositive_dimensions(self):
        buffers = load_module(
            "gpu_metc_buffers_validation_under_test", "pyoqp/oqp/utils/gpu_metc_buffers.py"
        )

        with self.assertRaisesRegex(ValueError, "nbf"):
            buffers.PersistentMetcBufferPlan.from_problem(
                nbf=0,
                nf=5,
                nmatrix=3,
                max_integrals=11,
            )


if __name__ == "__main__":
    unittest.main()
