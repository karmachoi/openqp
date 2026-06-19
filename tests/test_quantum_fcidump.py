"""Tests for the oqp.quantum second-quantized Hamiltonian / FCIDUMP bridge.

These exercise the pure-Python integral transforms and FCIDUMP I/O, which do
not require a compiled ``liboqp``. When the compiled extension is unavailable
(e.g. a source checkout without a build), the modules are loaded directly from
their files so the math is still covered.
"""

import importlib.util
import os

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_QDIR = os.path.join(_HERE, os.pardir, "pyoqp", "oqp", "quantum")


def _load(modname, filename):
    try:
        return __import__(f"oqp.quantum.{modname}", fromlist=[modname])
    except Exception:
        spec = importlib.util.spec_from_file_location(
            f"_oqp_quantum_{modname}", os.path.join(_QDIR, filename))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod


integrals = _load("integrals", "integrals.py")
fcidump = _load("fcidump", "fcidump.py")


# --------------------------------------------------------------------------
# Integral transforms
# --------------------------------------------------------------------------

def test_unpack_triangular_roundtrip():
    n = 4
    full = np.arange(n * n, dtype=float).reshape(n, n)
    sym = full + full.T
    packed = sym[np.tril_indices(n)]
    out = integrals.unpack_triangular(packed, n)
    assert np.allclose(out, sym)
    assert np.allclose(out, out.T)


def test_unpack_triangular_bad_size():
    with pytest.raises(ValueError):
        integrals.unpack_triangular(np.zeros(5), 4)


def test_ao_to_mo_1body_identity():
    h = np.array([[1.0, 0.2], [0.2, -0.5]])
    c = np.eye(2)
    assert np.allclose(integrals.ao_to_mo_1body(h, c), h)


def test_ao_to_mo_1body_matches_explicit():
    rng = np.random.default_rng(0)
    h = rng.standard_normal((3, 3))
    h = h + h.T
    c = rng.standard_normal((3, 3))
    ref = c.T @ h @ c
    assert np.allclose(integrals.ao_to_mo_1body(h, c), ref)


def test_ao_to_mo_2body_matches_full_einsum():
    rng = np.random.default_rng(1)
    n = 3
    eri = rng.standard_normal((n, n, n, n))
    # symmetrize to a physical chemist-notation tensor
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)
    c = rng.standard_normal((n, n))
    ref = np.einsum('pi,qj,rk,sl,ijkl->pqrs', c, c, c, c, eri, optimize=True)
    out = integrals.ao_to_mo_2body(eri, c)
    assert np.allclose(out, ref)


def test_ao_to_mo_2body_preserves_symmetry():
    rng = np.random.default_rng(2)
    n = 4
    a = rng.standard_normal((n, n))
    s = a + a.T
    eri = np.einsum('pq,rs->pqrs', s, s)  # valid (pq|rs) symmetry pattern
    c = rng.standard_normal((n, n))
    mo = integrals.ao_to_mo_2body(eri, c)
    assert np.allclose(mo, mo.transpose(1, 0, 2, 3))
    assert np.allclose(mo, mo.transpose(0, 1, 3, 2))
    assert np.allclose(mo, mo.transpose(2, 3, 0, 1))


# --------------------------------------------------------------------------
# FCIDUMP round trip
# --------------------------------------------------------------------------

def _random_hamiltonian(n, seed):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n))
    h1 = a + a.T
    g = rng.standard_normal((n, n, n, n))
    # impose full 8-fold permutation symmetry of (pq|rs)
    g = g + g.transpose(1, 0, 2, 3)
    g = g + g.transpose(0, 1, 3, 2)
    g = g + g.transpose(2, 3, 0, 1)
    return h1, g


def test_fcidump_write_read_roundtrip(tmp_path):
    n = 4
    h1, h2 = _random_hamiltonian(n, 3)
    ecore = -1.2345
    path = str(tmp_path / "test.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, ecore, n_electrons=4, ms2=0)

    data = fcidump.read_fcidump(path)
    assert data["norb"] == n
    assert data["nelec"] == 4
    assert data["ms2"] == 0
    assert np.isclose(data["ecore"], ecore)
    assert np.allclose(data["h1"], h1)
    assert np.allclose(data["h2"], h2)


def test_fcidump_header_format(tmp_path):
    n = 2
    h1, h2 = _random_hamiltonian(n, 5)
    path = str(tmp_path / "h.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, 0.7, n_electrons=2, ms2=0,
                          orbsym=[1, 1])
    with open(path) as fh:
        head = fh.read(200)
    assert "&FCI" in head
    assert "NORB=2" in head
    assert "NELEC=2" in head
    assert "ORBSYM=" in head
    assert "&END" in head


def test_fcidump_two_electron_count(tmp_path):
    # For norb=2 with all integrals non-zero, the 8-fold-unique count of
    # (pq|rs) pairs is the number of unique (P>=Q) compound-index pairs:
    # npair = 3 -> npair*(npair+1)/2 = 6 unique two-electron entries.
    n = 2
    h1, h2 = _random_hamiltonian(n, 9)
    path = str(tmp_path / "c.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, 0.0, n_electrons=2, tol=-1.0)
    twoe = 0
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 5 and int(p[3]) != 0 and int(p[4]) != 0:
                twoe += 1
    assert twoe == 6


def test_fcidump_tolerance_drops_small(tmp_path):
    n = 2
    h1 = np.array([[1.0, 1e-15], [1e-15, 2.0]])
    h2 = np.zeros((n, n, n, n))
    path = str(tmp_path / "t.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, 0.0, n_electrons=2)
    data = fcidump.read_fcidump(path)
    # tiny off-diagonal dropped -> reads back as exact zero
    assert data["h1"][0, 1] == 0.0
    assert np.isclose(data["h1"][0, 0], 1.0)
    assert np.isclose(data["h1"][1, 1], 2.0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
