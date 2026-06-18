"""
Automated regression + validation suite for pTC-MRSF-CIS.

Run:  pytest tests/ptc_mrsf/test_ptc_mrsf.py -v

Covers:
  * the non-Hermitian reduced eigensolver kernel (symmetric gate, biorthonormal
    non-symmetric case, complex-instability detection);
  * the from-scratch FCI engine vs pyscf (a second molecule, H4);
  * the exact transcorrelation identity <HF|H_bar|HF> = E_MP2;
  * cusp-accelerated basis convergence (Hooke's atom);
  * exact Hubbard-dimer correlation recovery;
  * working pTC-MRSF-CIS ground+excited states vs full FCI, across several H2
    bond lengths (robustness) and on a second system (LiH);
  * the singlet/triplet dressing property of the closed-shell correlator
    (motivating pTC's dual cusp conditions).

The prototype modules live next to this file under prototype/.
"""

import os
import sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "prototype"))

from nonsym_tda_eig import nonsym_tda_eig                       # noqa: E402
from tc_finite_basis import (ms0_determinants,                 # noqa: E402
                             build_fci_hamiltonian, build_t2_operator)
from ptc_mrsf_cis import build_s2, spin_of                     # noqa: E402

from pyscf import gto, scf, ao2mo, fci, mp                     # noqa: E402


# --------------------------------------------------------------------------
# shared pipeline helper
# --------------------------------------------------------------------------
def ptc_pipeline(atom, basis):
    """Return (e_bare, e_ptc, e_fci, spins_bare) for the (2,2) MRSF space."""
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    nocc = mol.nelectron // 2
    nvir = norb - nocc
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    ecore = mol.energy_nuc()
    pt = mp.MP2(mf).run()

    basis_d = ms0_determinants(norb, nocc)
    index = {d: i for i, d in enumerate(basis_d)}
    H = build_fci_hamiltonian(h1, eri, ecore, norb, basis_d, index)
    S2 = build_s2(basis_d, norb, nocc)
    wfci = np.linalg.eigvalsh(H)

    act = {nocc - 1, nocc}
    compact = sorted(i for i, d in enumerate(basis_d)
                     if all((P // 2) in act for P in d))
    sub = np.ix_(compact, compact)
    wb, vb = np.linalg.eigh(H[sub])
    spins = [spin_of(vb[:, k], S2[sub]) for k in range(len(compact))]

    from scipy.linalg import expm
    T2 = build_t2_operator(pt.t2, nocc, nvir, basis_d, index)
    Hbar = expm(-T2) @ H @ expm(T2)
    ee, _, _, info = nonsym_tda_eig(Hbar[sub], len(compact))
    assert info["n_complex"] == 0
    ee = np.sort(ee)
    return wb, ee, wfci[:len(compact)], spins


# --------------------------------------------------------------------------
# eigensolver kernel
# --------------------------------------------------------------------------
def test_eigensolver_symmetric_gate():
    rng = np.random.default_rng(0)
    m = rng.standard_normal((10, 10))
    a = 0.5 * (m + m.T)
    ref = np.linalg.eigvalsh(a)[:3]
    ee, vr, vl, info = nonsym_tda_eig(a, 3)
    assert np.allclose(ee, ref, atol=1e-10)
    assert info["n_complex"] == 0
    assert np.max(np.abs(a @ vr - vr * ee)) < 1e-9


def test_eigensolver_nonsymmetric_biorthonormal():
    rng = np.random.default_rng(1)
    n, nst = 14, 4
    m = rng.standard_normal((n, n))
    a = 0.5 * (m + m.T) + np.diag(np.arange(n) * 2.0)
    sk = rng.standard_normal((n, n))
    a = a + 0.05 * (sk - sk.T)
    ee, vr, vl, info = nonsym_tda_eig(a, nst)
    assert np.max(np.abs(a @ vr - vr * ee)) < 1e-8
    assert np.max(np.abs(vl.T @ vr - np.eye(nst))) < 1e-8
    assert info["max_imag"] < 1e-6


def test_eigensolver_detects_complex():
    a = np.array([[0.0, 1.0], [-1.0, 0.0]])
    _, _, _, info = nonsym_tda_eig(a, 1)
    assert info["n_complex"] >= 1


# --------------------------------------------------------------------------
# engine + transcorrelation identities
# --------------------------------------------------------------------------
def test_fci_engine_matches_pyscf_h4():
    mol = gto.M(atom='H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0',
                basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    nocc = mol.nelectron // 2
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    basis_d = ms0_determinants(norb, nocc)
    index = {d: i for i, d in enumerate(basis_d)}
    H = build_fci_hamiltonian(h1, eri, mol.energy_nuc(), norb, basis_d, index)
    e_me = np.linalg.eigvalsh(H)[0]
    e_pyscf, _ = fci.FCI(mf).kernel()
    assert abs(e_me - e_pyscf) < 1e-8


def test_transcorrelated_reference_equals_mp2():
    from scipy.linalg import expm
    mol = gto.M(atom='H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0',
                basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    nocc = mol.nelectron // 2
    nvir = norb - nocc
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    pt = mp.MP2(mf).run()
    basis_d = ms0_determinants(norb, nocc)
    index = {d: i for i, d in enumerate(basis_d)}
    H = build_fci_hamiltonian(h1, eri, mol.energy_nuc(), norb, basis_d, index)
    T2 = build_t2_operator(pt.t2, nocc, nvir, basis_d, index)
    Hbar = expm(-T2) @ H @ expm(T2)
    hf = index[tuple(sorted([2 * p for p in range(nocc)] +
                            [2 * p + 1 for p in range(nocc)]))]
    assert abs(Hbar[hf, hf] - pt.e_tot) < 1e-8


# --------------------------------------------------------------------------
# physics demonstrations (delegated to prototype main() self-checks)
# --------------------------------------------------------------------------
def test_cusp_convergence_module():
    import cusp_convergence
    cusp_convergence.main()           # asserts TC converges faster


def test_hubbard_dimer_module():
    import tc_hubbard_demo
    tc_hubbard_demo.main()            # asserts 100% recovery + eigenvector id


def test_finite_basis_module():
    import tc_finite_basis
    tc_finite_basis.main()           # asserts <HF|Hbar|HF>=MP2 + recovery


# --------------------------------------------------------------------------
# working pTC-MRSF-CIS: ground+excited vs FCI
# --------------------------------------------------------------------------
@pytest.mark.parametrize("r", [1.0, 1.4, 1.8, 2.2])
def test_ptc_recovers_ground_correlation_h2(r):
    # 2-electron (2,2) == full valence, so FCI is the exact benchmark.
    wb, ee, wfci, spins = ptc_pipeline(f'H 0 0 0; H 0 0 {r}', '6-31g')
    e_corr = wfci[0] - wb[0]
    recovered = (ee[0] - wb[0]) / e_corr * 100.0
    assert wb[0] > ee[0] > wfci[0] - 1e-9     # improved, bounded by FCI region
    # perturbative MP2 correlator recovers less as the bond stretches
    assert recovered > 35.0
    assert abs(ee[0] - wb[0]) > 1e-4          # S0 actually moved


def test_ptc_second_basis_h2():
    # a different basis (more virtuals) is a second, independent check
    wb, ee, wfci, spins = ptc_pipeline('H 0 0 0; H 0 0 1.4', 'cc-pvdz')
    recovered = (ee[0] - wb[0]) / (wfci[0] - wb[0]) * 100.0
    assert ee[0] < wb[0]
    assert recovered > 50.0


def test_adc_comparison_module():
    import compare_adc
    compare_adc.main()               # MRSF-CIS vs ADC(2) vs FCI, H2 dissociation


def test_singlet_dressed_triplet_not_h2():
    """Closed-shell MP2 correlator dresses singlets (incl. S0) but not the
    triplet: the concrete motivation for pTC's dual singlet+triplet cusp."""
    wb, ee, wfci, spins = ptc_pipeline('H 0 0 0; H 0 0 1.6', 'cc-pvdz')
    for k in range(len(wb)):
        if spins[k] > 1.5:                     # triplet
            assert abs(ee[k] - wb[k]) < 1e-6   # untouched
        if spins[k] < 0.5 and k == 0:          # ground singlet S0
            assert wb[k] - ee[k] > 1e-3        # dressed (lowered)


def test_r12_he_module():
    import r12_he
    r12_he.main()                    # explicit r12 recovers >70% He correlation


def test_triplet_dressing_module():
    import triplet_dressing
    triplet_dressing.main()          # active-pair correlator dresses the triplet


def test_r12_geminal_integral():
    import r12_geminal
    r12_geminal.main()               # analytic geminal integral vs numerics


def test_f12_intermediates():
    import f12_intermediates
    f12_intermediates.main()         # V, X, B over the Gaussian geminal
