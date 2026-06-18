"""
Exact demonstration of what transcorrelation BUYS you.

The molecular (2,2) script shows the non-Hermitian machinery runs on real
integrals, but on a *complete* active space a similarity transform cannot change
the spectrum, so it cannot show correlation recovery. This script does, on a
system small enough to solve exactly: the half-filled Hubbard dimer (2 sites,
2 electrons, hopping t, on-site repulsion U).

The point of transcorrelation (and of pTC) is that evaluating the
*transcorrelated* Hamiltonian H_bar = J^{-1} H J in a COMPACT reference recovers
correlation that the same compact reference misses with the bare H. Here the
compact reference is a single mean-field (RHF) determinant:

    E_RHF   = <g2|   H   |g2>            (mean field, misses correlation)
    E_TC(g) = <g2| J^-1 H J |g2>/<g2|g2> (transcorrelated, projective)

With a Gutzwiller correlation factor J = g^(double occupancy), the dimer is
exactly representable, so at the optimal g the mean-field determinant becomes the
exact right-eigenvector of H_bar and E_TC recovers 100% of the correlation
energy. This is the mechanism pTC-MRSF-CIS exploits.

Run:  python3 tc_hubbard_demo.py
"""

import numpy as np
from nonsym_tda_eig import nonsym_tda_eig

# spin-orbitals: 1up=0, 1dn=1, 2up=2, 2dn=3
SITE_OF = {0: 0, 1: 0, 2: 1, 3: 1}


def apply_adag_a(det, p, q):
    """c^dag_p c_q |det>; det = sorted tuple of occupied spin-orbitals.
    Returns (sign, new_det) or (0, None)."""
    if q not in det:
        return 0, None
    occ = list(det)
    sign = (-1) ** occ.index(q)
    occ.remove(q)
    if p in occ:
        return 0, None
    # insertion sign
    pos = 0
    while pos < len(occ) and occ[pos] < p:
        pos += 1
    sign *= (-1) ** pos
    occ.insert(pos, p)
    return sign, tuple(occ)


def build_hubbard(t, U):
    """4-determinant Ms=0 Hubbard-dimer Hamiltonian and the basis."""
    basis = [(0, 1), (2, 3), (0, 3), (1, 2)]  # |1u1d>,|2u2d>,|1u2d>,|2u1d>
    idx = {d: i for i, d in enumerate(basis)}
    dim = len(basis)
    H = np.zeros((dim, dim))
    hops = [(0, 2), (2, 0), (1, 3), (3, 1)]  # c^dag_p c_q pairs for -t
    for i, det in enumerate(basis):
        # U term: double occupancy of a site
        d_occ = sum(1 for s in (0, 1)
                    if (2 * s) in det and (2 * s + 1) in det)
        H[i, i] += U * d_occ
        # hopping -t
        for (p, q) in hops:
            sign, nd = apply_adag_a(det, p, q)
            if sign != 0 and nd in idx:
                H[idx[nd], i] += -t * sign
    return H, basis


def double_occupancy_vec(basis):
    return np.array([sum(1 for s in (0, 1)
                         if (2 * s) in d and (2 * s + 1) in d)
                     for d in basis], dtype=float)


def rhf_determinant(basis):
    """|g^2> with g = (site1 + site2)/sqrt(2), in the determinant basis."""
    idx = {d: i for i, d in enumerate(basis)}
    psi = np.zeros(len(basis))
    # g^dag_up g^dag_dn = 1/2 (c1u+c2u)(c1d+c2d)
    terms = [((0, 1), +1), ((0, 3), +1), ((2, 1), +1), ((2, 3), +1)]
    for (det, _) in terms:
        d = tuple(sorted(det))
        # sign from ordering c^dag_a c^dag_b |vac> -> sorted
        a, b = det
        sign = 1.0 if a < b else -1.0
        psi[idx[d]] += 0.5 * sign
    return psi


def main():
    t, U = 1.0, 4.0
    H, basis = build_hubbard(t, U)

    # exact
    w = np.linalg.eigvalsh(H)
    e_exact = w[0]
    e_exact_closed = 0.5 * (U - np.sqrt(U**2 + 16 * t**2))
    assert abs(e_exact - e_exact_closed) < 1e-10, (e_exact, e_exact_closed)

    # mean-field reference
    g2 = rhf_determinant(basis)
    g2 /= np.linalg.norm(g2)
    e_rhf = g2 @ H @ g2
    e_rhf_closed = 0.5 * U - 2.0 * t
    assert abs(e_rhf - e_rhf_closed) < 1e-10, (e_rhf, e_rhf_closed)

    e_corr = e_exact - e_rhf  # negative

    print("=== Hubbard dimer (t=%.1f, U=%.1f) ===" % (t, U))
    print(f"exact ground energy        : {e_exact:.8f}  (closed form {e_exact_closed:.8f})")
    print(f"mean-field (RHF) energy    : {e_rhf:.8f}  (closed form {e_rhf_closed:.8f})")
    print(f"correlation energy         : {e_corr:.8f}\n")

    # The projective energy <g2|J^-1 H J|g2> is non-variational (unbounded as
    # g->0), so g is NOT fixed by minimizing it. It is fixed by the
    # transcorrelation condition that the reference satisfy the transformed
    # equations -- equivalently here, the variational Gutzwiller optimum
    #     E_G(g) = <g2| J H J |g2> / <g2| J^2 |g2>      (Hermitian, bounded).
    docc = double_occupancy_vec(basis)
    gs = np.linspace(0.05, 1.0, 1000)
    e_g = []
    for g in gs:
        jd = g ** docc
        jpsi = jd * g2
        e_g.append((jpsi @ H @ jpsi) / (jpsi @ jpsi))
    e_g = np.array(e_g)
    kmin = int(np.argmin(e_g))
    g_opt, e_g_min = gs[kmin], e_g[kmin]
    recovered = (e_g_min - e_rhf) / e_corr * 100.0

    print("Gutzwiller transcorrelation E_G(g) (bounded, variational):")
    print(f"  optimal g                : {g_opt:.4f}")
    print(f"  E_G(g_opt)               : {e_g_min:.8f}")
    print(f"  correlation recovered    : {recovered:.2f} %")
    assert e_g_min >= e_exact - 1e-9, "variational energy must bound E_exact"
    assert recovered > 99.0, recovered
    print("  VALIDATED: a single mean-field determinant recovers ~100% of the")
    print("  correlation energy through the correlation factor J.\n")

    # Key identity: at g_opt the mean-field determinant is the EXACT
    # right-eigenvector of the (non-Hermitian) transcorrelated H_bar:
    #   H |J g2> = E |J g2>  <=>  H_bar |g2> = E |g2>,  H_bar = J^-1 H J.
    jd = g_opt ** docc
    hbar = (1.0 / jd)[:, None] * H * jd[None, :]
    e_proj = (g2 @ hbar @ g2) / (g2 @ g2)
    resid = np.linalg.norm(hbar @ g2 - e_proj * g2)
    print("transcorrelated projective check at g_opt:")
    print(f"  projective energy <g2|H_bar|g2>   : {e_proj:.8f}")
    print(f"  |H_bar|g2> - E|g2>|               : {resid:.2e}")
    print(f"  matches exact ground energy       : {abs(e_proj - e_exact):.2e}")
    assert abs(e_proj - e_exact) < 5e-3, (e_proj, e_exact)
    assert resid < 5e-3, resid
    print("  VALIDATED: mean-field determinant IS the right-eigenvector of H_bar.\n")

    # the non-Hermitian solver returns the exact spectrum from H_bar
    hbar = (1.0 / jd)[:, None] * H * jd[None, :]
    assert not np.allclose(hbar, hbar.T)
    ee, vr, vl, info = nonsym_tda_eig(hbar, len(basis))
    print("non-Hermitian solve of H_bar(g_opt):")
    print(f"  spectrum real (max Im)   : {info['max_imag']:.2e}")
    print(f"  matches exact spectrum   : {np.max(np.abs(np.sort(ee)-np.sort(w))):.2e}")
    print(f"  biorthonormality         : "
          f"{np.max(np.abs(vl.T @ vr - np.eye(len(basis)))):.2e}")
    assert np.allclose(np.sort(ee), np.sort(w), atol=1e-9)


if __name__ == "__main__":
    main()
