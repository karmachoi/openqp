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

    def test_parse_timer_lines_extracts_only_timer_records_from_log_text(self):
        log_text = "\n".join(
            [
                "normal OpenQP output",
                self.timers.format_timer_line("tdhf.davidson.total", 2.5, {"iter": 3}),
                "unrelated warning line",
                self.timers.format_timer_line("tdhf.response.total", 3.0),
            ]
        )

        records = self.timers.parse_timer_lines(log_text)

        self.assertEqual([record["label"] for record in records], ["tdhf.davidson.total", "tdhf.response.total"])
        self.assertEqual(records[0]["iter"], "3")
        self.assertAlmostEqual(records[1]["seconds"], 3.0)

    def test_summarize_timer_records_groups_counts_and_total_seconds_by_label(self):
        records = [
            self.timers.parse_timer_line(self.timers.format_timer_line("tdhf.davidson.total", 2.5, {"iter": 1})),
            self.timers.parse_timer_line(self.timers.format_timer_line("tdhf.davidson.total", 3.5, {"iter": 2})),
            self.timers.parse_timer_line(self.timers.format_timer_line("tdhf.response.total", 9.0)),
        ]

        summary = self.timers.summarize_timer_records(records)

        self.assertEqual(summary["tdhf.davidson.total"], {"count": 2, "seconds_total": 6.0, "seconds_mean": 3.0})
        self.assertEqual(summary["tdhf.response.total"], {"count": 1, "seconds_total": 9.0, "seconds_mean": 9.0})


if __name__ == "__main__":
    unittest.main()
