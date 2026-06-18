"""
Validated NumPy reference for the non-Hermitian (transcorrelated) MRSF-CIS
reduced-space eigensolver.

Context
-------
pTC-MRSF-CIS replaces the bare electronic Hamiltonian with Ten-no's projective
transcorrelated effective Hamiltonian H_bar = e^{-tau} H e^{tau}. H_bar is
non-Hermitian, so the MRSF-CIS (Tamm-Dancoff) response matrix A is no longer
real-symmetric. The current OpenQP solver (`rpaeig`, TDA branch in
source/tdhf_lib.F90) assumes A = A^T and uses `diag_symm_packed`. This module
prototypes the replacement: a general (non-symmetric) reduced-space solve that
returns real eigenvalues with biorthonormal left/right Ritz vectors.

This is the *reduced-space* kernel (the small nvec x nvec subspace matrix that
the Davidson iteration builds), i.e. the drop-in replacement for the TDA branch
of `rpaeig`. The outer Davidson machinery (subspace build, residual, new-vector
generation) is unchanged in structure; only the subspace diagonalization and the
biorthonormal residual change.

Phase-1 validation gate
-----------------------
With the transcorrelation switched off (tau = 0), A is symmetric and this
non-symmetric solver MUST reproduce the symmetric eigenpairs bit-for-bit (up to
sign/degeneracy). That equivalence is the gate that lets us trust all the
non-Hermitian plumbing before any transcorrelated integral is computed.

Port target: source/modules/tdhf_mrsf_ptc.F90 :: tc_nonsym_tda_eig
LAPACK mapping: np.linalg.eig  ->  DGEEV (left+right eigenvectors).
"""

import numpy as np


def nonsym_tda_eig(amat, nstate, imag_tol=1.0e-8):
    """Lowest-`nstate` right eigenpairs of a (generally) non-symmetric A.

    Parameters
    ----------
    amat : (n, n) real array
        Reduced-space MRSF-CIS response matrix. Symmetric when tau = 0,
        non-symmetric for the transcorrelated Hamiltonian.
    nstate : int
        Number of low-lying roots requested.
    imag_tol : float
        Tolerance on |Im(eigenvalue)| / (1+|Re|) below which a root is treated
        as physically real (bound state). Complex pairs above this signal a
        non-Hermitian instability and are reported.

    Returns
    -------
    ee  : (nstate,) real      eigenvalues, ascending
    vr  : (n, nstate) real    right eigenvectors  (A   vr = vr ee)
    vl  : (n, nstate) real    left  eigenvectors  (A^T vl = vl ee)
                              biorthonormal: vl[:,i] . vr[:,j] = delta_ij
    info: dict                diagnostics (max imaginary part, n_complex)
    """
    n = amat.shape[0]
    # DGEEV analogue: right (vr) and left (wl) eigenvectors.
    w, vr_c = np.linalg.eig(amat)
    wl, vl_c = np.linalg.eig(amat.T)

    # Match left vectors to right eigenvalues by nearest eigenvalue.
    order_r = np.argsort(w.real)
    order_l = np.argsort(wl.real)
    w = w[order_r]
    vr_c = vr_c[:, order_r]
    vl_c = vl_c[:, order_l]  # same ordering -> same eigenvalue per column

    max_imag = float(np.max(np.abs(w.imag) / (1.0 + np.abs(w.real))))
    n_complex = int(np.sum(np.abs(w.imag) / (1.0 + np.abs(w.real)) > imag_tol))

    # Bound states are real; take real parts.
    ee_full = w.real
    vr_full = vr_c.real
    vl_full = vl_c.real

    # Biorthonormalize: scale so that vl_i . vr_i = 1.
    for i in range(n):
        denom = vl_full[:, i] @ vr_full[:, i]
        if abs(denom) < 1.0e-14:
            # left/right nearly orthogonal -> defective/near-degenerate; guard.
            denom = np.sign(denom) * 1.0e-14 if denom != 0 else 1.0e-14
        scal = 1.0 / np.sqrt(abs(denom))
        sgn = np.sign(denom)
        vr_full[:, i] *= scal
        vl_full[:, i] *= scal * sgn

    ee = ee_full[:nstate].copy()
    vr = vr_full[:, :nstate].copy()
    vl = vl_full[:, :nstate].copy()
    info = {"max_imag": max_imag, "n_complex": n_complex}
    return ee, vr, vl, info


# ---------------------------------------------------------------------------
# Tests (run: python3 nonsym_tda_eig.py)
# ---------------------------------------------------------------------------
def _sym_reference(amat, nstate):
    """What the existing symmetric `rpaeig` TDA branch computes."""
    ee, v = np.linalg.eigh(amat)
    return ee[:nstate], v[:, :nstate]


def test_tau0_reduces_to_symmetric():
    """Phase-1 gate: tau=0 (symmetric A) must reproduce the symmetric solve."""
    rng = np.random.default_rng(0)
    n, nstate = 12, 4
    m = rng.standard_normal((n, n))
    a_sym = 0.5 * (m + m.T)  # symmetric == tau=0 MRSF-CIS

    ee_ref, _ = _sym_reference(a_sym, nstate)
    ee, vr, vl, info = nonsym_tda_eig(a_sym, nstate)

    assert np.allclose(ee, ee_ref, atol=1e-10), (ee, ee_ref)
    assert info["n_complex"] == 0
    # eigenpair residual A vr - vr ee = 0
    res = a_sym @ vr - vr * ee
    assert np.max(np.abs(res)) < 1e-9, np.max(np.abs(res))
    print("PASS  tau=0 reduces to symmetric   max|dE|=%.2e  res=%.2e"
          % (np.max(np.abs(ee - ee_ref)), np.max(np.abs(res))))


def test_nonsymmetric_real_spectrum():
    """Transcorrelated case: small non-symmetric perturbation, real spectrum."""
    rng = np.random.default_rng(1)
    n, nstate = 16, 5
    m = rng.standard_normal((n, n))
    a_sym = 0.5 * (m + m.T)
    a_sym += np.diag(np.arange(n) * 2.0)         # well-separated diagonal
    skew = rng.standard_normal((n, n))
    a_tc = a_sym + 0.05 * (skew - skew.T)        # non-Hermitian TC perturbation

    ee, vr, vl, info = nonsym_tda_eig(a_tc, nstate)

    # right-eigenpair residual
    res_r = a_tc @ vr - vr * ee
    # left-eigenpair residual
    res_l = a_tc.T @ vl - vl * ee
    # biorthonormality vl^T vr = I (lowest block)
    gram = vl.T @ vr
    assert np.max(np.abs(res_r)) < 1e-8, np.max(np.abs(res_r))
    assert np.max(np.abs(res_l)) < 1e-8, np.max(np.abs(res_l))
    assert np.max(np.abs(gram - np.eye(nstate))) < 1e-8, gram
    assert info["max_imag"] < 1e-6, info
    print("PASS  non-symmetric real spectrum   res_R=%.2e res_L=%.2e "
          "biorth=%.2e Im=%.2e"
          % (np.max(np.abs(res_r)), np.max(np.abs(res_l)),
             np.max(np.abs(gram - np.eye(nstate))), info["max_imag"]))


def test_detects_complex_instability():
    """A strongly non-Hermitian A can have complex pairs -> must be flagged."""
    a = np.array([[0.0, 1.0], [-1.0, 0.0]])  # eigenvalues +-i
    _, _, _, info = nonsym_tda_eig(a, 1)
    assert info["n_complex"] >= 1, info
    print("PASS  complex instability detected  n_complex=%d max_imag=%.2e"
          % (info["n_complex"], info["max_imag"]))


if __name__ == "__main__":
    test_tau0_reduces_to_symmetric()
    test_nonsymmetric_real_spectrum()
    test_detects_complex_instability()
    print("\nAll Phase-1 eigensolver-kernel tests passed.")
