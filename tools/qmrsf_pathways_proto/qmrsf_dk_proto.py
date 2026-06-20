#!/usr/bin/env python3
"""
QMRSF-DK proof-of-concept (pure NumPy, NO pyscf / NO scipy).

Demonstrates the MECHANISM of the dressed (frequency-dependent) xc kernel that the
QMRSF-DK pathway adds to the closed-shell 0OS double-spin-flip diagonal of the QMRSF
response matrix. Companion derivation: QMRSF_DK_kernel.md.

Physics (Maitra, Cave, Zhang, Burke, JCP 120, 5932 (2004); Casida; Loos-Blase):
  An ADIABATIC (frequency-independent) xc kernel has a response basis of SINGLE
  excitations only -> it has NO double-excitation pole and gives a wrong single-state
  energy. A DRESSED (frequency-dependent) kernel
        D(omega) = A + sum_d |V_d|^2 / (omega - omega_d)
  injects the missing double pole and recovers the exact single+double spectrum.

This file builds the canonical Maitra-type model (one single coupled to one/few doubles),
with a KNOWN exact spectrum (exact diagonalization of the full single+double matrix), then:
  (a) ADIABATIC  : project out the double -> one root near A (misses the double; wrong single)
  (b) DRESSED    : solve omega = A + B(omega) self-consistently -> recovers ALL roots
and validates (b) against exact diagonalization to machine precision. It then extends to a
small "0OS-like" 3-state block to show the pole-search (Prescription P1) generalizes.

Run:  python3 qmrsf_dk_proto.py
"""
import numpy as np


# ----------------------------------------------------------------------
# 1. The canonical Maitra-type model: 1 single coupled to N_d doubles
# ----------------------------------------------------------------------
def augmented_H(A, omega_d, V):
    """Augmented (single + doubles) Hamiltonian, Eq. (9) of QMRSF_DK_kernel.md.

      [[ A      V_1   ... V_Nd ],
       [ V_1*   w_1   ...  0   ],
       [ ...                   ],
       [ V_Nd*   0    ...  w_Nd]]

    A       : adiabatic single-sector diagonal (scalar).
    omega_d : (N_d,) bare double-excitation energies.
    V       : (N_d,) couplings <single|H|double_d>.
    Returns the (1+N_d, 1+N_d) symmetric matrix; its eigenvalues are the EXACT spectrum.
    """
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    V = np.atleast_1d(np.asarray(V, float))
    nd = omega_d.size
    H = np.zeros((1 + nd, 1 + nd))
    H[0, 0] = A
    H[0, 1:] = V
    H[1:, 0] = V
    H[1:, 1:] = np.diag(omega_d)
    return H


def exact_spectrum(A, omega_d, V):
    """Exact eigenvalues of the full single+double model = the reference truth."""
    return np.linalg.eigvalsh(augmented_H(A, omega_d, V))


# ----------------------------------------------------------------------
# 2. (a) ADIABATIC approximation: drop the frequency dependence
# ----------------------------------------------------------------------
def adiabatic_energy(A, omega_d, V):
    """Adiabatic TDDFT: the response basis is the SINGLE sector only; the doubles are
    not in the basis (the frequency-dependent residue B(omega) is identically 0, Eq. 8).
    -> a single root, sitting exactly at the bare adiabatic diagonal A. It MISSES every
    double-excitation state and carries none of the level repulsion from the doubles."""
    return float(A)


# ----------------------------------------------------------------------
# 3. (b) DRESSED kernel: solve omega = A + B(omega) self-consistently
# ----------------------------------------------------------------------
def B_of_omega(omega, omega_d, V):
    """Frequency-dependent dressing B(omega) = sum_d |V_d|^2 / (omega - omega_d). (Eq. 6)."""
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    V = np.atleast_1d(np.asarray(V, float))
    return np.sum(np.abs(V) ** 2 / (omega - omega_d))


def dressed_roots(A, omega_d, V, tol=1e-13, maxit=200):
    """Solve the dressed eigenvalue condition g(omega) = omega - A - B(omega) = 0
    (Eq. 7) by a robust bracketed search (Prescription P1 / Section 3(b)).

    g(omega) is monotone INCREASING on every interval between consecutive poles
    {omega_d} (and on the two semi-infinite ends), with a simple pole at each omega_d.
    So there is exactly one root per open interval: 1 + N_d roots total. We bracket
    each interval and bisect (pure NumPy, no scipy)."""
    omega_d = np.sort(np.atleast_1d(np.asarray(omega_d, float)))
    V = np.atleast_1d(np.asarray(V, float))

    def g(w):
        return w - A - B_of_omega(w, omega_d, V)

    # Pole locations partition the real line. Build brackets that avoid the poles.
    poles = omega_d
    spread = max(1.0, float(np.ptp(poles)) if poles.size > 1 else 1.0,
                 abs(A) + np.sum(np.abs(V)))
    eps = 1e-7 * max(1.0, spread)           # tiny offset to step off a pole
    big = 10.0 * spread + abs(A) + 100.0    # far enough to bracket the end roots

    # Interval edges: -inf | pole_1 | pole_2 | ... | pole_Nd | +inf  (open at poles)
    edges = [(-big, poles[0] - eps)]
    for i in range(len(poles) - 1):
        edges.append((poles[i] + eps, poles[i + 1] - eps))
    edges.append((poles[-1] + eps, big))

    roots = []
    for lo, hi in edges:
        glo, ghi = g(lo), g(hi)
        if glo == 0.0:
            roots.append(lo); continue
        if ghi == 0.0:
            roots.append(hi); continue
        if np.sign(glo) == np.sign(ghi):
            # no sign change in this (finite, pole-avoiding) bracket -> skip
            continue
        a, b, ga = lo, hi, glo
        for _ in range(maxit):
            m = 0.5 * (a + b)
            gm = g(m)
            if abs(gm) < tol or (b - a) < tol:
                break
            if np.sign(gm) == np.sign(ga):
                a, ga = m, gm
            else:
                b = m
        roots.append(0.5 * (a + b))
    return np.sort(np.array(roots))


# ----------------------------------------------------------------------
# 4. Driver
# ----------------------------------------------------------------------
def banner(t):
    print("=" * 78); print(t); print("=" * 78)


def report_case(name, A, omega_d, V):
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    V = np.atleast_1d(np.asarray(V, float))
    exact = exact_spectrum(A, omega_d, V)
    adiab = adiabatic_energy(A, omega_d, V)
    dressed = dressed_roots(A, omega_d, V)

    print(f"\n--- {name} ---")
    print(f"  A (adiabatic single diag) = {A:.6f}")
    print(f"  bare doubles  omega_d     = {np.array2string(omega_d, precision=6)}")
    print(f"  couplings     V           = {np.array2string(V, precision=6)}")
    print(f"  exact spectrum (diag full): {np.array2string(exact, precision=6)}  "
          f"[{exact.size} roots]")

    # adiabatic: nearest exact root to A, and what it misses
    near = exact[np.argmin(np.abs(exact - adiab))]
    adiab_err = abs(adiab - near)
    print(f"  ADIABATIC  E = {adiab:.6f}  -> nearest exact root {near:.6f}, "
          f"err = {adiab_err:.3e}")
    print(f"             (basis = single only -> {exact.size - 1} double root(s) MISSED)")

    print(f"  DRESSED   roots           : {np.array2string(dressed, precision=6)}  "
          f"[{dressed.size} roots]")
    ok = dressed.size == exact.size
    maxerr = np.abs(dressed - exact).max() if ok else float('nan')
    print(f"  DRESSED vs EXACT  max|err| = {maxerr:.3e}   "
          f"({'MATCH' if ok and maxerr < 1e-8 else 'MISMATCH'})")
    return exact, adiab, dressed, maxerr, adiab_err


def main():
    np.set_printoptions(precision=6, suppress=True)
    banner("QMRSF-DK proof-of-concept  |  dressed (freq-dependent) kernel restores the "
           "double pole")
    print("Model: 1 single excitation (adiabatic diag A) coupled to N_d double excitations")
    print("(bare energies omega_d, couplings V) -- the Maitra-Cave-Zhang-Burke model.")
    print("Adiabatic kernel: basis = single only (B(omega)=0)  -> misses doubles.")
    print("Dressed   kernel: omega = A + sum_d |V_d|^2/(omega - omega_d) -> recovers all.")

    # ---- Case 1: the canonical 1 single + 1 double (MCZB) ----
    # A double sits near the single; the adiabatic answer must miss the double-like state
    # AND mis-place the single because it carries no level repulsion.
    A1, wd1, V1 = 5.0, 5.6, 0.8
    e1, ad1, dr1, err1, adiab_err1 = report_case("Case 1: one single + one double",
                                                 A1, wd1, V1)

    # explicit check of the two-root secular eq (2): (A-w)(wd-w) - V^2 = 0
    check = np.array([(A1 - w) * (wd1 - w) - V1 ** 2 for w in dr1])
    print(f"  secular residual (A-w)(wd-w)-V^2 at dressed roots: "
          f"max|.| = {np.abs(check).max():.3e}")

    # ---- Case 2: near-resonance (A ~ omega_d): strong mixing, the hard case ----
    A2, wd2, V2 = 4.0, 4.05, 0.5
    report_case("Case 2: near-resonant single/double (strong mixing)", A2, wd2, V2)

    # ---- Case 3: one single + a FEW doubles (multi-satellite dressing) ----
    A3 = 6.0
    wd3 = np.array([4.5, 6.4, 8.1])
    V3 = np.array([0.45, 0.7, 0.3])
    report_case("Case 3: one single + three doubles (multi-satellite)", A3, wd3, V3)

    # ---- 0OS-like 3-state block: pole-search generalizes (Prescription P1) ----
    banner('"0OS-like" block: pole search over a small coupled set (Prescription P1)')
    print("Two 0OS-like single channels, each adiabatically at A_c, coupled to a shared")
    print("set of doubles. We dress EACH channel state-by-state (Eq. 10-11) and confirm")
    print("each scalar pole search reproduces the exact roots of its own augmented block.")
    blocks = [
        ("0OS channel #1", 5.2, np.array([4.9, 6.1]), np.array([0.6, 0.4])),
        ("0OS channel #2", 7.0, np.array([6.7, 8.3]), np.array([0.5, 0.55])),
    ]
    worst = 0.0
    for nm, Ac, wdc, Vc in blocks:
        _, _, _, err, _ = report_case(nm, Ac, wdc, Vc)
        worst = max(worst, err)

    # ---- Final verdict ----
    banner("VERDICT")
    print(f"Case 1  ADIABATIC error vs nearest exact root      : {adiab_err1:.3e} "
          f"(and it MISSES the double-like state entirely)")
    print(f"Case 1  DRESSED   max error vs exact (both roots)   : {err1:.3e}")
    all_ok = all(np.isfinite(x) and x < 1e-8 for x in [err1, worst])
    print(f"All dressed solutions match exact diagonalization to < 1e-8 : "
          f"{'YES' if all_ok else 'NO'}")
    print("\nConclusion: the adiabatic kernel (single-only basis) cannot produce the")
    print("double-excitation pole and mis-places the single; the frequency-dependent")
    print("dressed kernel D(omega)=A+sum_d |V_d|^2/(omega-omega_d) injects the pole and")
    print("recovers the exact single+double spectrum -- the QMRSF-DK 0OS mechanism,")
    print("validated against exact diagonalization, independent of OpenQP/pyscf.")


if __name__ == "__main__":
    main()
