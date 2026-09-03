"""Gates for the analytic MRSF transition-dipole derivative.

Two layers:

* pure-Python algebra gates, which always run: the packed/unfolded adjoint
  identity and the exactness of the bilinear dipole source;
* a calculation gate that runs H2O/MRSF-BHHLYP/6-31G* in process and checks
  the analytic derivative against a frozen reference, the analytic AO
  dipole-derivative translational sum rule, rigid-translation invariance of
  d(mu_IJ)/dR, and the same-space rotational gauge residual of the assembled
  orbital source.

The calculation gate is skipped unless OQP_HT_RUN=1 AND the imported oqp
package really carries the derivative-dipole entry point.  The skip is loud on
purpose: a silently no-op gate has cost this project time before.
"""

from __future__ import annotations

import importlib.util
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "pyoqp" / "oqp" / "library" / "htdipole_analytic.py"
_SPEC = importlib.util.spec_from_file_location("htdipole_analytic", MODULE_PATH)
HT = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(HT)


H2O_INPUT = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf
conv=1.0e-11
maxit=200

[tdhf]
type=mrsf
nstate=5
conv=1.0e-11
maxit=300
zvconv=1.0e-10
"""

# Frozen reference for the (S0, S1) pair, a.u. per bohr, rows = x/y/z of mu,
# columns = the nine Cartesian nuclear coordinates.  Reproduced by central
# finite differences of the phase-transported transition dipole at
# h = 0.01/0.005/0.0025 bohr with clean O(h**2) convergence (deviations
# 2.31e-5 / 5.81e-6 / 1.50e-6).  The overall sign is the MRSF state-phase
# convention and is not reproducible between runs, so it is projected out.
H2O_S0S1_DMU = np.array([
    [-0.0, 0.0, -0.11746507, 0.12477437, -0.03604514,
     0.05873254, -0.12477437, 0.03604514, 0.05873254],
    [-0.0, 0.0, -0.11746507, 0.03604514, -0.12477437,
     0.05873254, -0.03604514, 0.12477437, 0.05873254],
    [-0.16519604, -0.16519604, 0.0, 0.08259802, 0.08259802,
     -0.0, 0.08259802, 0.08259802, 0.0],
])
H2O_S0S1_ABS_MU = 0.17898304


def _oqp_has_dipole_derivatives() -> bool:
    try:
        import oqp
    except Exception:
        return False
    return hasattr(oqp, "electric_dipole_der_bra")


RUN_CALC = os.environ.get("OQP_HT_RUN") == "1"
HAVE_OQP = _oqp_has_dipole_derivatives()


class _StubContext(HT.MRSFTransitionDipole):
    """Fold/algebra-only context; no molecule handle is touched."""

    def __init__(self, nbf, noca, nocb, seed=0):
        rng = np.random.default_rng(seed)
        self.mol = None
        self.nbf, self.noca, self.nocb = nbf, noca, nocb
        self.nvirb = nbf - nocb
        self.nij = noca*self.nvirb
        self.mult = 1
        self.C = np.eye(nbf)
        sym = rng.normal(size=(3, nbf, nbf))
        self.R = [0.5*(s + s.T) for s in sym]
        self.m = list(self.R)
        self.energies = np.zeros(1)
        self.nstate = 1
        self.bvec = np.zeros((self.nij, 1))
        self.ilr1 = (noca - 1 - nocb - 1)*noca + noca - 2
        self.ilr2 = (noca - nocb - 1)*noca + noca - 1
        self.redundant = self.ilr2
        somo = (noca - 2, noca - 1)
        self.mask_vv = np.zeros((noca, self.nvirb))
        self.mask_oo = np.zeros((noca, self.nvirb))
        for k in somo:
            self.mask_vv[k, 0:2] = 1.0
            self.mask_oo[k, 0:2] = 1.0


class TestTransitionDipoleAlgebra(unittest.TestCase):
    """Algebra that must hold independently of any OpenQP build."""

    def setUp(self):
        self.ctx = _StubContext(nbf=19, noca=6, nocb=4)
        self.rng = np.random.default_rng(11)

    def test_missing_response_roots_flags_a_skipped_dense_root(self):
        rng = np.random.default_rng(3)
        n = 12
        q, _ = np.linalg.qr(rng.standard_normal((n, n)))
        spectrum = np.array([-0.1, 0.04, 0.2, 0.216, 0.217, 0.279, 0.291,
                             0.306, 0.313, 0.40, 0.43, 0.44])
        amat = q @ np.diag(spectrum) @ q.T
        amat[5, :] = 0.0; amat[:, 5] = 0.0      # a redundant slot
        keep = [i for i in range(n) if i != 5]
        dense = np.linalg.eigvalsh(amat[np.ix_(keep, keep)])
        # Davidson "found" the eight lowest roots except the one near 0.306
        found = np.delete(dense[:9], 7)
        missing = HT.missing_response_roots(amat, found, 5)
        self.assertEqual(len(missing), 1)
        self.assertAlmostEqual(missing[0], dense[7], places=10)
        self.assertEqual(HT.missing_response_roots(amat, dense[:9], 5), [])

    def _packed(self):
        b = self.rng.normal(size=self.ctx.nij)
        b[self.ctx.redundant] = 0.0
        return b

    def test_fold_adjoint_is_the_adjoint_of_unfold(self):
        """<fold_adjoint(g), b> == <g, unfold(b)> for every g and b."""
        for _ in range(20):
            b = self._packed()
            g = self.rng.normal(size=(self.ctx.noca, self.ctx.nvirb))
            lhs = float(self.ctx.fold_adjoint(g) @ b)
            rhs = float(np.sum(g*self.ctx.unfold(b)))
            self.assertAlmostEqual(lhs, rhs, places=12)

    def test_redundant_slot_never_leaves_the_fold(self):
        b = self._packed()
        b[self.ctx.redundant] = 7.0
        clean = np.array(b, copy=True)
        clean[self.ctx.redundant] = 0.0
        np.testing.assert_allclose(self.ctx.unfold(b), self.ctx.unfold(clean))
        g = self.rng.normal(size=(self.ctx.noca, self.ctx.nvirb))
        self.assertEqual(self.ctx.fold_adjoint(g)[self.ctx.redundant], 0.0)

    def test_transition_dipole_is_symmetric_in_the_two_states(self):
        xi = self.ctx.unfold(self._packed())
        xj = self.ctx.unfold(self._packed())
        np.testing.assert_allclose(self.ctx.transition_dipole(xi, xj),
                                   self.ctx.transition_dipole(xj, xi),
                                   atol=1e-12)

    def test_transition_density_is_traceless(self):
        """Origin independence of mu_IJ rests on this, at every geometry."""
        xi = self.ctx.unfold(self._packed())
        xj = self.ctx.unfold(self._packed())
        self.assertLess(abs(np.trace(self.ctx.gamma_mo(xi, xj))), 1e-12)

    def test_dipole_source_is_the_exact_bilinear_adjoint(self):
        """mu is bilinear, so the source must reproduce it with no step error."""
        bi, bj = self._packed(), self._packed()
        xi, xj = self.ctx.unfold(bi), self.ctx.unfold(bj)
        for k in range(3):
            w = self.ctx.fold_adjoint(self.ctx.dipole_source(xi, xj, k))
            for _ in range(5):
                d = self._packed()
                lhs = float(w @ d)
                rhs = (self.ctx.transition_dipole(self.ctx.unfold(bi + d), xj)[k]
                       - self.ctx.transition_dipole(xi, xj)[k])
                self.assertAlmostEqual(lhs, rhs, places=10)


@unittest.skipUnless(
    RUN_CALC and HAVE_OQP,
    "set OQP_HT_RUN=1 and use a build that exports electric_dipole_der_bra; "
    "this gate runs a real MRSF calculation and is not a no-op stub")
class TestTransitionDipoleDerivative(unittest.TestCase):
    """In-process calculation gates on H2O/MRSF-BHHLYP/6-31G*."""

    @classmethod
    def setUpClass(cls):
        import oqp
        from oqp.pyoqp import Runner
        cls.tmp = tempfile.TemporaryDirectory()
        inp = Path(cls.tmp.name) / "h2o_ht.inp"
        inp.write_text(H2O_INPUT)
        runner = Runner(input_file=str(inp), log=str(inp.with_suffix(".log")))
        runner.run()
        cls.mol = runner.mol
        cls.result = HT.analytic_transition_dipole_derivative(
            cls.mol, pairs=[(0, 1)])
        cls.diag = cls.result["diagnostics"][(0, 1)]
        cls.ctx = HT.MRSFTransitionDipole(cls.mol)
        oqp.electric_dipole_der_bra(cls.mol)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_dipole_derivative_integrals_obey_the_translation_identity(self):
        """sum_A dR^k/dx_c = delta_kc S -- the dipole operator has a fixed origin."""
        nbf = self.ctx.nbf
        natom = int(self.mol.data["natom"])
        dbra = HT._fortran_view(
            np.array(self.mol.data["OQP::dip_dbra"], copy=True),
            (nbf, nbf, 3, 3*natom))
        deriv = dbra + dbra.transpose(1, 0, 2, 3)
        packed = np.asarray(self.mol.data["OQP::SM"]).ravel()
        overlap = HT._unpack_lt(packed, nbf)
        for k in range(3):
            for c in range(3):
                total = sum(deriv[:, :, k, 3*a + c] for a in range(natom))
                target = overlap if k == c else np.zeros_like(overlap)
                self.assertLess(np.abs(total - target).max(), 1e-10)

    def test_transition_dipole_matches_the_frozen_reference(self):
        amps = [self.ctx.unfold(self.ctx.bvec[:, s]) for s in range(2)]
        mu = self.ctx.transition_dipole(amps[0], amps[1])
        self.assertAlmostEqual(abs(mu[0]), H2O_S0S1_ABS_MU, places=6)
        self.assertAlmostEqual(abs(mu[1]), H2O_S0S1_ABS_MU, places=6)
        self.assertLess(abs(mu[2]), 1e-8)

    def test_analytic_derivative_matches_the_finite_difference_reference(self):
        got = self.result[(0, 1)]
        sign = np.sign(float(np.sum(got*H2O_S0S1_DMU))) or 1.0
        self.assertLess(np.abs(sign*got - H2O_S0S1_DMU).max(), 5e-6)

    def test_rigid_translation_leaves_the_transition_dipole_unchanged(self):
        got = self.result[(0, 1)]
        natom = int(self.mol.data["natom"])
        total = got.reshape(3, natom, 3).sum(axis=1)
        self.assertLess(np.abs(total).max(), 1e-10)

    def test_same_space_rotational_gauge_residual_vanishes(self):
        """The orbital source must have no same-space antisymmetric part.

        The interchange machinery fixes those rotations to zero; the residual
        is therefore a direct check that the amplitude response is complete."""
        for entry in self.diag:
            self.assertLess(entry["gauge_residual"],
                            1e-5*max(entry["source_scale"], 1.0))


if __name__ == "__main__":
    unittest.main()
