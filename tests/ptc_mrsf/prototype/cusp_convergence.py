"""
Phase 2b: the headline transcorrelation claim, demonstrated exactly.

The whole point of an explicit correlation factor (F12 / transcorrelation, and
hence pTC) is that building the electron-electron cusp into the wavefunction
accelerates basis-set convergence: a smooth orbital basis converges slowly to the
non-smooth (cusped) exact wavefunction, whereas an explicitly-correlated basis
removes the cusp and converges much faster.

We show this exactly on Hooke's atom (two electrons in a harmonic well with
Coulomb repulsion), whose relative-motion problem separates into a 1D radial
equation with a real r12 cusp [psi'(0)/psi(0) = 1/2]:

    H_u = -lap_u + (1/4) w^2 u^2 + 1/u ,    l = 0.

A high-accuracy finite-difference solve gives the exact energy. We then expand the
relative wavefunction in n s-Gaussians, comparing:
    * bare basis  {g_k}            (smooth, no cusp),
    * correlated  {J * g_k}        (J = exp[u/(2(1+b u))], J'(0)=1/2: correct cusp).

The correlated basis reaches a given accuracy with far fewer functions -- the
basis-set-convergence acceleration that pTC-MRSF-CIS inherits.

Run:  python3 cusp_convergence.py
"""

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal

W = 1.0


def V(u):
    return 0.25 * W * W * u * u + 1.0 / u


def fd_reference(umax=20.0, N=20000):
    """Exact relative-motion ground state by finite differences (chi = u psi)."""
    u = np.linspace(umax / N, umax, N)
    h = u[1] - u[0]
    diag = 2.0 / h**2 + V(u)
    off = -1.0 / h**2 * np.ones(N - 1)
    e, _ = eigh_tridiagonal(diag, off, select='i', select_range=(0, 0))
    return e[0]


def make_grid(umax=25.0, N=40000):
    u = np.linspace(1e-6, umax, N)
    return u


def jastrow(u, b=0.5):
    """Pade-Jastrow with the exact coalescence cusp J'(0)/J(0) = 1/2."""
    s = u / (2.0 * (1.0 + b * u))
    sp = 1.0 / (2.0 * (1.0 + b * u) ** 2)
    J = np.exp(s)
    return J, J * sp


def solve(alphas, u, correlated):
    """Lowest eigenvalue of the radial generalized eigenproblem (measure u^2 du).

    Kinetic uses the gradient form  <h_i|-lap|h_j> = int u^2 h_i' h_j' du."""
    J, Jp = jastrow(u)
    funcs = []
    for a in alphas:
        g = np.exp(-a * u * u)
        gp = -2.0 * a * u * g
        if correlated:
            funcs.append((J * g, Jp * g + J * gp))
        else:
            funcs.append((g, gp))
    n = len(alphas)
    H = np.zeros((n, n))
    S = np.zeros((n, n))
    w2 = u * u
    Vu = V(u)
    for i in range(n):
        hi, hpi = funcs[i]
        for j in range(i, n):
            hj, hpj = funcs[j]
            S[i, j] = S[j, i] = np.trapezoid(w2 * hi * hj, u)
            T = np.trapezoid(w2 * hpi * hpj, u)
            Vv = np.trapezoid(w2 * Vu * hi * hj, u)
            H[i, j] = H[j, i] = T + Vv
    return eigh(H, S, eigvals_only=True)[0]


def expset(n):
    return 0.05 * (2.2 ** np.arange(n))


def main():
    e_exact = fd_reference()
    u = make_grid()
    print("=== Hooke's atom: cusp-accelerated basis convergence ===")
    print(f"exact (finite difference): {e_exact:.8f}\n")
    print("  n     bare E        TC E          err_bare    err_TC    speedup")
    err_bare, err_tc = {}, {}
    for n in range(2, 9):
        al = expset(n)
        eb = solve(al, u, correlated=False)
        et = solve(al, u, correlated=True)
        err_bare[n] = abs(eb - e_exact)
        err_tc[n] = abs(et - e_exact)
        print(f"  {n}   {eb:.6f}   {et:.6f}   {err_bare[n]:.2e}  "
              f"{err_tc[n]:.2e}   {err_bare[n]/err_tc[n]:6.1f}x")

    # Validations: the correlated basis converges markedly faster.
    assert err_tc[4] < err_bare[4] / 10.0, (err_tc[4], err_bare[4])
    # correlated basis at n=4 already beats the bare basis at n=8
    assert err_tc[4] < err_bare[8], (err_tc[4], err_bare[8])
    print("\nVALIDATED: the explicit cusp factor accelerates basis-set")
    print("convergence -- the correlated basis at n=4 is more accurate than")
    print("the bare basis at n=8. This is the convergence benefit pTC-MRSF-CIS")
    print("inherits from the transcorrelated Hamiltonian.")


if __name__ == "__main__":
    main()
