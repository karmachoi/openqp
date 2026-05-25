import ast
import os
import sys
import unittest
from pathlib import Path
from typing import Callable, cast


PYOQP_SOURCE = Path(__file__).resolve().parents[1] / "pyoqp" / "oqp" / "pyoqp.py"


class ThreadingDefaultsTest(unittest.TestCase):
    def _load_helper_from_source(self) -> Callable[[], None]:
        tree = ast.parse(PYOQP_SOURCE.read_text())
        helper_nodes = []
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "_set_threading_defaults":
                helper_nodes.append(node)
            elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                func = node.value.func
                if isinstance(func, ast.Name) and func.id == "_set_threading_defaults":
                    helper_nodes.append(node)
                    break
        module = ast.Module(body=helper_nodes, type_ignores=[])
        ast.fix_missing_locations(module)
        namespace = {"os": os, "sys": sys}
        exec(compile(module, str(PYOQP_SOURCE), "exec"), namespace)
        return cast(Callable[[], None], namespace["_set_threading_defaults"])

    def test_threading_defaults_are_set_before_loading_oqp(self):
        text = PYOQP_SOURCE.read_text()
        self.assertLess(text.index("_set_threading_defaults()"), text.index("import oqp"))

    def test_threading_defaults_preserve_user_values(self):
        helper = self._load_helper_from_source()
        keys = [
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_STACKSIZE",
            "GOMP_STACKSIZE",
        ]
        old = {key: os.environ.get(key) for key in keys}
        try:
            for key in keys:
                os.environ.pop(key, None)
            os.environ["OPENBLAS_NUM_THREADS"] = "8"
            os.environ["OMP_NUM_THREADS"] = "24"
            helper()
            self.assertEqual(os.environ["OMP_NUM_THREADS"], "24")
            self.assertEqual(os.environ["OPENBLAS_NUM_THREADS"], "8")
            self.assertEqual(os.environ["MKL_NUM_THREADS"], "1")
            self.assertEqual(os.environ["BLIS_NUM_THREADS"], "1")
            self.assertEqual(os.environ["VECLIB_MAXIMUM_THREADS"], "1")
            self.assertEqual(os.environ["OMP_STACKSIZE"], "256M")
            self.assertEqual(os.environ["GOMP_STACKSIZE"], "256M")
        finally:
            for key, value in old.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
    def test_threading_defaults_set_safe_macos_omp_count(self):
        helper = self._load_helper_from_source()
        old = os.environ.get("OMP_NUM_THREADS")
        try:
            os.environ.pop("OMP_NUM_THREADS", None)
            helper()
            if sys.platform == "darwin":
                self.assertEqual(os.environ["OMP_NUM_THREADS"], "16")
            else:
                self.assertNotIn("OMP_NUM_THREADS", os.environ)
        finally:
            if old is None:
                os.environ.pop("OMP_NUM_THREADS", None)
            else:
                os.environ["OMP_NUM_THREADS"] = old


if __name__ == "__main__":
    unittest.main()
