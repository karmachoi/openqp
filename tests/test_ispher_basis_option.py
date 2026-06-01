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


def install_oqpdata_stubs():
    oqp = types.ModuleType("oqp")
    oqp.ffi = object()
    oqp.lib = object()
    sys.modules["oqp"] = oqp

    periodic_table = types.ModuleType("oqp.periodic_table")
    periodic_table.MASSES = {}
    periodic_table.SYMBOL_MAP = {}
    sys.modules["oqp.periodic_table"] = periodic_table

    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    constants = types.ModuleType("oqp.utils.constants")
    constants.ANGSTROM_TO_BOHR = 1.8897259886
    sys.modules["oqp.utils.constants"] = constants


def load_oqpdata():
    install_oqpdata_stubs()
    return load_module("oqpdata_ispher_under_test", "pyoqp/oqp/molecule/oqpdata.py")


class TestISPHERBasisOption(unittest.TestCase):
    def test_schema_defines_ispher_default_cartesian(self):
        oqpdata = load_oqpdata()

        self.assertIn("ispher", oqpdata.OQP_CONFIG_SCHEMA["input"])
        entry = oqpdata.OQP_CONFIG_SCHEMA["input"]["ispher"]
        self.assertEqual(entry["default"], "-1")
        self.assertIs(entry["type"], oqpdata.ispher)
        self.assertEqual(oqpdata.ispher("-1"), -1)
        self.assertEqual(oqpdata.ispher("0"), 0)
        self.assertEqual(oqpdata.ispher("1"), 1)
        with self.assertRaisesRegex(ValueError, "input.ispher.*-1, 0, or 1"):
            oqpdata.ispher("2")

    def test_basis_shell_dimensions_encode_cartesian_and_pure_counts(self):
        text = (ROOT / "pyoqp/oqp/library/set_basis.py").read_text()

        self.assertIn("cartesian_shell_function_count", text)
        self.assertIn("pure_shell_function_count", text)
        self.assertIn("basis_shell_function_count", text)
        self.assertIn("ispher", text)

    def test_ispher_zero_compatibility_notice_is_user_visible(self):
        text = (ROOT / "pyoqp/oqp/library/set_basis.py").read_text()

        self.assertIn("ISPHER=0", text)
        self.assertIn("Cartesian-equivalent compatibility mode", text)
        self.assertIn("SALC", text)

    def test_c_header_exposes_ispher_control_field_for_cffi_runtime(self):
        header = (ROOT / "include/oqp.h").read_text()

        self.assertIn("struct control_parameters", header)
        self.assertIn("int64_t   ispher;", header)

    def test_native_mapping_rejects_unimplemented_pure_basis_explicitly(self):
        source = (ROOT / "source/basis_api.F90").read_text()

        self.assertIn("pure/spherical basis functions are not implemented", source)
        self.assertIn("ISPHER=1", source)
        self.assertIn("with_abort", source)

    def test_example_and_docs_include_ispher_keyword(self):
        example = ROOT / "examples/BASIS/H2O_RHF_ISPHER1_PURE_BASIS.inp"
        docs = ROOT / "docs/dev/ispher_basis_option_plan.md"

        self.assertTrue(example.exists(), "missing ISPHER=1 user-facing example")
        self.assertIn("ispher=1", example.read_text().lower())
        self.assertTrue(docs.exists(), "missing ISPHER development plan")
        docs_text = docs.read_text().lower()
        self.assertIn("ispher=-1", docs_text)
        self.assertIn("ispher=0", docs_text)
        self.assertIn("ispher=1", docs_text)
        self.assertIn("cartesian-equivalent", docs_text)

    def test_user_input_manual_documents_ispher_keyword_and_claim_boundary(self):
        readme = ROOT / "pyoqp/README.md"

        self.assertTrue(readme.exists(), "missing PyOQP input manual")
        text = readme.read_text().lower()
        self.assertIn("ispher=-1", text)
        self.assertIn("ispher=0", text)
        self.assertIn("ispher=1", text)
        self.assertIn("cartesian-equivalent", text)
        self.assertIn("pure/spherical", text)
        self.assertIn("not implemented", text)


if __name__ == "__main__":
    unittest.main()
