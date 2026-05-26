import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class TDHFBenchmarkTimerTests(unittest.TestCase):
    def setUp(self):
        self.timers = load_module("tdhf_benchmark_timers_under_test", "pyoqp/oqp/utils/tdhf_benchmark_timers.py")

    def test_davidson_timer_manifest_has_stable_gpu_benchmark_labels(self):
        manifest = self.timers.davidson_timer_manifest()
        labels = [timer.label for timer in manifest]

        self.assertEqual(labels[0], "tdhf.response.total")
        self.assertIn("tdhf.davidson.total", labels)
        self.assertIn("tdhf.davidson.sigma_build", labels)
        self.assertIn("tdhf.davidson.metc_contract", labels)
        self.assertIn("tdhf.davidson.eri_buffer", labels)
        self.assertTrue(all(timer.private_gpu_group for timer in manifest))

    def test_format_timer_line_is_machine_parseable_and_unit_explicit(self):
        line = self.timers.format_timer_line(
            "tdhf.davidson.sigma_build",
            elapsed_seconds=1.25,
            metadata={"branch": "perf/tdhf-davidson-timers", "nstate": 4},
        )

        self.assertEqual(
            line,
            "OQP_TIMER label=tdhf.davidson.sigma_build seconds=1.250000 branch=perf/tdhf-davidson-timers nstate=4",
        )

    def test_parse_timer_line_round_trips_formatted_timer(self):
        line = self.timers.format_timer_line(
            "tdhf.davidson.metc_contract",
            elapsed_seconds=0.0032,
            metadata={"kernel": "cpu_baseline"},
        )

        parsed = self.timers.parse_timer_line(line)

        self.assertEqual(parsed["label"], "tdhf.davidson.metc_contract")
        self.assertAlmostEqual(parsed["seconds"], 0.0032)
        self.assertEqual(parsed["kernel"], "cpu_baseline")


if __name__ == "__main__":
    unittest.main()
