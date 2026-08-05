"""Collect the petite-list smoke scripts so CI actually runs them.

The three ``tests/smoke_petite_*.py`` scripts were written as validation for
the symmetry reduction and have never executed in CI: ``pyproject.toml`` sets
``testpaths = ["tests"]``, but pytest's default ``python_files`` is
``test_*.py`` / ``*_test.py``, which ``smoke_*.py`` does not match. Nothing
else invokes them either -- no workflow, no documentation, no Makefile.

That is not a cosmetic gap. ``smoke_petite_benchmark.py`` already asserts
``status == 'active'`` on two cc-pVTZ cases, so it encodes exactly the check
that would have caught the reduction being silently disabled on every
spherical basis. The test existed; it had simply never run.

This module is a thin collector. The scripts keep their standalone entry
points (they are useful to run by hand with a workdir argument) and keep
owning their own gates; here we only assert their exit status.

Cost: the two validation scripts together take about 6 s. The benchmark takes
several minutes because it runs acenes and cc-pVTZ, so it is opt-in via
``OQP_RUN_PETITE_BENCHMARK=1`` rather than charged to every CI run.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def _runtime_available():
    """True when the compiled OpenQP runtime can actually be driven."""
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _run_script(name, workdir):
    """Run a smoke script in its own process and return (rc, output).

    A separate process is deliberate: these scripts compare two calculations
    against each other, so any state that survives between jobs inside one
    interpreter shows up as a fake difference. That is not hypothetical --
    the Z-vector warm start does exactly this, see
    smoke_petite_gradient_validation.py.
    """
    proc = subprocess.run(
        [sys.executable, str(HERE / f"{name}.py"), str(workdir)],
        cwd=str(ROOT), capture_output=True, text=True, timeout=1800,
        env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "2")},
    )
    return proc.returncode, proc.stdout + proc.stderr


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class PetiteSmokeScripts(unittest.TestCase):

    def _check(self, name):
        import tempfile
        with tempfile.TemporaryDirectory(prefix=f"oqp_{name}_") as tmp:
            rc, out = _run_script(name, tmp)
        if rc != 0:
            # The scripts print one line per case with its own gate; surface
            # those rather than the whole config dump.
            verdicts = [l for l in out.splitlines()
                        if l.startswith(('[ok', '[FAIL')) or 'FAILED' in l]
            self.fail(f"{name} failed (rc={rc}):\n" + "\n".join(verdicts[-12:]))

    def test_petite_list_energies_match_c1(self):
        """Skeleton-Fock energies vs C1 references, gate 1e-10."""
        self._check("smoke_petite_list_validation")

    def test_petite_gradients_match_c1(self):
        """Petite gradients vs C1 references, gate 1e-9."""
        self._check("smoke_petite_gradient_validation")

    @unittest.skipUnless(os.environ.get("OQP_RUN_PETITE_BENCHMARK") == "1",
                         "set OQP_RUN_PETITE_BENCHMARK=1 to run the slow benchmark")
    def test_petite_benchmark(self):
        """Speedup/accuracy across acenes and both symmetry tiers (minutes)."""
        self._check("smoke_petite_benchmark")


if __name__ == "__main__":
    unittest.main()
