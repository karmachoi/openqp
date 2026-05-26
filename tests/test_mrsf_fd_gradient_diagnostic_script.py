from pathlib import Path
import importlib.util
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "fd_mrsf_gradient_check.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("fd_mrsf_gradient_check", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MrsfFdGradientDiagnosticScriptTest(unittest.TestCase):
    def test_script_exists_with_expected_variants_and_default_h2o_geometry(self):
        self.assertTrue(SCRIPT.exists(), f"missing diagnostic script: {SCRIPT}")
        mod = _load_script()
        self.assertEqual(mod.VARIANTS, ("default", "spc_ovov_0", "flip_ovov_sign"))
        self.assertEqual(len(mod.BUILTIN_GEOMETRIES["h2o"]), 3)
        self.assertEqual(mod.BUILTIN_GEOMETRIES["h2o"][0][0], "O")

    def test_input_template_keeps_mrsf_root_numbering_and_variant_controls(self):
        mod = _load_script()
        text = mod.render_input(
            geom=mod.BUILTIN_GEOMETRIES["h2o"],
            runtype="grad",
            target_root=3,
            variant="spc_ovov_0",
            basis="3-21g",
            functional="bhhlyp",
            nstate=4,
            z_solver=1,
            huckel=True,
        )
        self.assertIn("method=tdhf", text)
        self.assertIn("[scf]\n", text)
        self.assertIn("type=rohf", text)
        self.assertIn("multiplicity=3", text)
        self.assertIn("[tdhf]\n", text)
        self.assertIn("type=mrsf", text)
        self.assertIn("nstate=4", text)
        self.assertIn("z_solver=1", text)
        self.assertIn("spc_ovov=0", text)
        self.assertIn("[properties]\ngrad=3", text)

    def test_fd_gradient_conversion_uses_angstrom_displacement_to_bohr_gradient(self):
        mod = _load_script()
        # E(+h)-E(-h)=0.002 Ha at h=0.002 Ang gives 0.5 Ha/Ang.
        # Convert to Ha/Bohr by dividing by Bohr per Angstrom.
        self.assertAlmostEqual(
            mod.central_fd_ha_per_bohr(1.001, 0.999, 0.002),
            0.5 / mod.BOHR_PER_ANG,
            places=12,
        )


if __name__ == "__main__":
    unittest.main()
