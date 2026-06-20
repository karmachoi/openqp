#!/usr/bin/env python3
"""
QMRSF-DK LIVE pathway prototype (pure NumPy, NO pyscf / NO scipy).

This is the "live-shaped" prototype for the density-functional / dressed-kernel
realization of QMRSF-DK -- the object the eventual live module
source/modules/tdhf_qmrsf_dk.F90 must reproduce. It takes EXACTLY the inputs the
live module will assemble from OpenQP tagarray quantities:

    A0       : the (Ns x Ns) ADIABATIC single-spin-flip response block over the
               quintet reference's frontier (2OS/4OS sectors). In the live module
               this is the converged backbone singles block restricted to the
               single-spin-flip response space (real, symmetric).
    omega_d  : the (Nd,) bare 0OS closed-shell double-spin-flip energies (the
               diagonal of the 0OS sector; from the active-space integrals).
    V        : the (Ns x Nd) couplings <single c | H | 0OS double d> between the
               single-spin-flip block and the injected 0OS doubles.

and produces ALL physical roots: the DRESSED single-like states AND the INJECTED
0OS double-like states, via the frequency-dependent dressed kernel
    A(omega) = A0 + g_xc(omega),   g_xc(omega)_{c c'} = sum_d V_{c,d} V_{c',d}/(omega - omega_d).
The dressed problem is solved by the secular determinant / pole search
    fsec(omega) = det[ omega I - A0 - g_xc(omega) ] = 0,
which is the multi-channel form of the scalar secular function
    fsec(omega) = omega - A - sum_d V_d^2/(omega - omega_d)
solved by bracketed bisection in qmrsf_dk_proto.py / qmrsf_dk_block_proto.py.

KEY POINT vs the icPT2 pathway: there is NO determinant CI here and NO external-Q
PT2. The 0OS doubles are injected as EXTRA ROOTS purely through the
frequency-dependent kernel poles; the adiabatic (omega-independent) kernel
A(omega)->A0 has ONLY Ns roots and structurally MISSES the Nd doubles.

GATES (printed PASS/FAIL):
  (1) dressed-kernel roots  == exact eigenvalues of the explicit augmented
      [[A0, V],[V^T, diag(omega_d)]] matrix, to < 1e-9  (ALL Ns+Nd roots);
  (2) adiabatic kernel (omega-independent) MISSES the doubles: it yields only Ns
      roots and the most-doubly-excited exact state is absent (gap >> 1e-3).

This file also writes nothing; the Fortran cross-check model is dumped separately
by dump_dk_ref.py (mirroring dump_icpt2_downfold_ref.py).

Run:  python3 qmrsf_dk_live_proto.py
"""
import numpy as np


# ======================================================================
# 0. A minimal CAS(4,4)-like LIVE model
# ======================================================================
def build_live_model(seed=20260620):
    """A small, faithful stand-in for what the live module assembles.

    Ns single-spin-flip channels (the 2OS/4OS frontier response block A0) and
    Nd 0OS closed-shell double-spin-flip configurations (energies omega_d,
    couplings V into the singles). Sizes echo the QMRSF CAS(4,4) structure:
    a handful of single-spin-flip channels and the SIX 0OS doubles -- but the
    DK mechanism is size-agnostic, so we keep it small and well-conditioned.

    Returns A0 (Ns,Ns) symmetric, omega_d (Nd,), V (Ns,Nd).
    """
    rng = np.random.default_rng(seed)
    Ns, Nd = 5, 6                      # 5 single-spin-flip channels, 6 0OS doubles

    # --- adiabatic single-spin-flip block A0: symmetric, diagonally dominant ---
    # diagonals ~ frontier single-spin-flip excitation energies (eV-like scale),
    # small off-diagonal couplings among the single channels.
    diagA = np.array([3.2, 4.1, 4.8, 5.9, 6.6])
    A0 = np.diag(diagA)
    off = 0.15 * (rng.standard_normal((Ns, Ns)))
    off = 0.5 * (off + off.T)
    np.fill_diagonal(off, 0.0)
    A0 = A0 + off                      # real symmetric singles block

    # --- bare 0OS double energies: spread across and above the singles window ---
    # deliberately place some doubles INSIDE the singles band (the hard case the
    # adiabatic kernel cannot reach) and some above it.
    omega_d = np.array([4.45, 5.35, 6.2, 7.1, 8.0, 9.3])

    # --- singles<->0OS couplings V (Ns x Nd): moderate, all channels coupled ---
    V = 0.30 + 0.20 * rng.standard_normal((Ns, Nd))
    V *= 0.6                            # keep mixing moderate so P-search is clean

    return A0, omega_d, V


# ======================================================================
# 1. Explicit augmented matrix (the reference truth) + exact spectrum
# ======================================================================
def augmented_matrix(A0, omega_d, V):
    """The EXPLICIT (singles + 0OS doubles) matrix whose eigenvalues are the
    reference truth:
        [[ A0 ,  V            ],
         [ V^T,  diag(omega_d) ]].
    The DK dressed kernel is the EXACT Feshbach downfold of this matrix onto the
    singles sector -- so dressed roots must reproduce eigvalsh(this) exactly."""
    A0 = np.asarray(A0, float)
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    V = np.asarray(V, float).reshape(A0.shape[0], omega_d.size)
    ns, nd = A0.shape[0], omega_d.size
    M = np.zeros((ns + nd, ns + nd))
    M[:ns, :ns] = A0
    M[:ns, ns:] = V
    M[ns:, :ns] = V.T
    M[ns:, ns:] = np.diag(omega_d)
    return M


def exact_spectrum(A0, omega_d, V):
    return np.linalg.eigvalsh(augmented_matrix(A0, omega_d, V))


def doubles_weight(A0, omega_d, V):
    """Doubles-sector weight of every exact eigenvector (to identify which exact
    states are 'predominantly double' -- the ones the adiabatic kernel misses)."""
    ns = A0.shape[0]
    w, vecs = np.linalg.eigh(augmented_matrix(A0, omega_d, V))
    return w, (vecs[ns:, :] ** 2).sum(axis=0)


# ======================================================================
# 2. The dressed (frequency-dependent) kernel and the secular determinant
# ======================================================================
def gxc(omega, omega_d, V):
    """Frequency-dependent quadratic kernel g_xc(omega) in the singles sector:
        g_xc(omega)_{c c'} = sum_d V_{c,d} V_{c',d} / (omega - omega_d).
    This is the matrix generalization of B(omega) = |V|^2/(omega-omega_d). The
    poles at {omega_d} are what inject the missing 0OS double-excitation roots."""
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    V = np.asarray(V, float)
    d = omega - omega_d
    return (V / d) @ V.T                  # sum_d V[:,d] V[:,d]^T / (omega-omega_d)


def fsec(omega, A0, omega_d, V):
    """Pole-cancelled secular function of the dressed eigenvalue condition.

    The raw condition is det[ omega I - A0 - g_xc(omega) ] = 0. Multiplying through
    by prod_d (omega - omega_d) clears every pole and leaves a SMOOTH polynomial-
    like function with the SAME roots (the augmented spectrum) and NO singularities:

        fsec(omega) = det[ omega I - A0 - g_xc(omega) ] * prod_d (omega - omega_d).

    This is exactly the characteristic polynomial of the augmented matrix
    [[A0, V],[V^T, diag(omega_d)]] (Schur-complement / Feshbach identity), so its
    Ns+Nd real roots ARE the exact eigenvalues. Being pole-free and smooth, it is
    trivially and robustly bracketed by a sign scan -- and it ports verbatim to
    Fortran (no slogdet, no eigvals needed in the search). Returned scaled by a
    fixed normalization to keep magnitudes tame."""
    A0 = np.asarray(A0, float)
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    V = np.asarray(V, float)
    ns = A0.shape[0]
    M = omega * np.eye(ns) - A0 - gxc(omega, omega_d, V)
    detM = np.linalg.det(M)
    return detM * np.prod(omega - omega_d)


# ======================================================================
# 3. DRESSED roots via the (pole-cancelled) secular function: the EXACT
#    multi-channel Feshbach downfold -- Prescription P0.
# ======================================================================
def dressed_roots(A0, omega_d, V, tol=1e-13, maxit=300):
    """Solve the dressed eigenvalue condition for ALL Ns+Nd roots.

    Uses the pole-cancelled secular function fsec(omega) (= the augmented
    characteristic polynomial), which is smooth with no singularities, so every
    root is a clean sign change. We bracket the whole spectral window with a fine
    uniform scan (the spectrum lives in the Gershgorin range of the augmented
    matrix) and bisect each sign change. This is the EXACT downfold (no
    off-diagonal-B approximation) and must reproduce eigvalsh(augmented_matrix)
    to machine precision.

    Returns (sorted roots, expected count Ns+Nd)."""
    A0 = np.asarray(A0, float)
    omega_d = np.sort(np.atleast_1d(np.asarray(omega_d, float)))
    V = np.asarray(V, float)
    ns, nd = A0.shape[0], omega_d.size
    ntot = ns + nd

    def f(w):
        return fsec(w, A0, omega_d, V)

    # Gershgorin bounds of the augmented matrix bracket the entire real spectrum.
    aug = augmented_matrix(A0, omega_d, V)
    radii = np.abs(aug).sum(axis=1) - np.abs(np.diag(aug))
    lo = float((np.diag(aug) - radii).min()) - 1.0
    hi = float((np.diag(aug) + radii).max()) + 1.0

    # Fine uniform scan over the bounded window: resolution must be finer than the
    # closest root spacing. ntot roots in [lo,hi] -> use >> ntot samples.
    nscan = max(20000, 2000 * ntot)
    xs = np.linspace(lo, hi, nscan)
    fv = np.array([f(x) for x in xs])

    roots = []
    for i in range(nscan - 1):
        if fv[i] == 0.0:
            roots.append(xs[i]); continue
        if fv[i] * fv[i + 1] < 0.0:
            a, b, fa = xs[i], xs[i + 1], fv[i]
            for _ in range(maxit):
                m = 0.5 * (a + b)
                fm = f(m)
                if (b - a) < tol or fm == 0.0:
                    break
                if fm * fa > 0.0:
                    a, fa = m, fm
                else:
                    b = m
            roots.append(0.5 * (a + b))

    roots = np.sort(np.array(roots))
    # numerical de-dup: collapse roots closer than 1e-10
    if roots.size > 1:
        keep = [roots[0]]
        for r in roots[1:]:
            if r - keep[-1] > 1e-10:
                keep.append(r)
        roots = np.array(keep)
    return roots, ntot


# ======================================================================
# 4. ADIABATIC kernel (omega-independent): g_xc -> g_xc(0-limit dropped) = 0
# ======================================================================
def adiabatic_roots(A0):
    """Adiabatic TDDFT/MRSF: the kernel carries NO frequency dependence, so the
    response basis is the SINGLES sector only and the dressing is identically
    zero (g_xc -> 0). The spectrum is just eigvalsh(A0): Ns single roots, and the
    Nd 0OS double-excitation states are structurally MISSING."""
    return np.linalg.eigvalsh(np.asarray(A0, float))


# ======================================================================
# 5. Driver / gates
# ======================================================================
def banner(t):
    print("=" * 80); print(t); print("=" * 80)


def main():
    np.set_printoptions(precision=6, suppress=True, linewidth=120)
    banner("QMRSF-DK LIVE prototype | dressed kernel A(omega)=A0+g_xc(omega) "
           "injects the 0OS doubles")
    print("Inputs mirror the live module: A0 (single-spin-flip block), omega_d (0OS")
    print("doubles), V (single<->0OS couplings). Dressed kernel poles inject the")
    print("0OS double-excitation roots the adiabatic kernel structurally misses.")
    print("NO determinant CI, NO external-Q PT2 -- pure frequency-dependent kernel.\n")

    A0, omega_d, V = build_live_model()
    Ns, Nd = A0.shape[0], omega_d.size

    exact = exact_spectrum(A0, omega_d, V)
    w_all, dbl_w = doubles_weight(A0, omega_d, V)
    dressed, ntot = dressed_roots(A0, omega_d, V)
    adiab = adiabatic_roots(A0)

    print(f"  Ns single-spin-flip channels = {Ns}")
    print(f"  Nd 0OS doubles               = {Nd}")
    print(f"  A0 eigenvalues (adiabatic)   = {np.array2string(np.sort(adiab), precision=5)}")
    print(f"  bare 0OS doubles omega_d     = {np.array2string(np.sort(omega_d), precision=5)}")
    print(f"  EXACT spectrum (explicit)    = {np.array2string(exact, precision=5)}  "
          f"[{exact.size} roots]")
    print(f"  DRESSED roots (pole search)  = {np.array2string(dressed, precision=5)}  "
          f"[{dressed.size} roots]")
    print(f"  ADIABATIC roots (g_xc=0)     = {np.array2string(np.sort(adiab), precision=5)}  "
          f"[{adiab.size} roots]")

    # ---- GATE 1: dressed roots == exact, ALL Ns+Nd, to < 1e-9 ----
    count_ok = (dressed.size == exact.size == ntot)
    maxerr = float(np.abs(dressed - exact).max()) if count_ok else float('nan')
    gate1 = count_ok and maxerr < 1e-9
    print(f"\n  [GATE 1] dressed vs exact: count {dressed.size}=={exact.size}? "
          f"{count_ok}; max|err| = {maxerr:.3e}  -> {'PASS' if gate1 else 'FAIL'}")

    # ---- GATE 2: adiabatic MISSES the doubles ----
    # The adiabatic kernel has only Ns roots; the most-doubly-excited exact states
    # (largest doubles-weight) have no adiabatic partner within 1e-3.
    order = np.argsort(dbl_w)[::-1]
    most_dbl_E = w_all[order[:Nd]]                 # the Nd most-double exact states
    # nearest adiabatic root to each most-double exact state
    near_gap = np.array([np.min(np.abs(adiab - E)) for E in most_dbl_E])
    n_missed_by_adiab = exact.size - adiab.size    # = Nd by construction
    worst_dbl_gap = float(near_gap.max())
    # the dressed kernel, by contrast, reproduces these double states:
    dressed_dbl_gap = float(np.array([np.min(np.abs(dressed - E)) for E in most_dbl_E]).max())
    gate2 = (n_missed_by_adiab == Nd) and (worst_dbl_gap > 1e-3) and (dressed_dbl_gap < 1e-9)
    print(f"  [GATE 2] adiabatic has {adiab.size} roots, EXACT has {exact.size} "
          f"-> {n_missed_by_adiab} double(s) absent.")
    print(f"           most-doubly-excited exact states E = "
          f"{np.array2string(np.sort(most_dbl_E), precision=5)} "
          f"(doubles-weight {np.array2string(np.sort(dbl_w[order[:Nd]])[::-1], precision=3)})")
    print(f"           worst gap nearest-ADIABATIC-root = {worst_dbl_gap:.3e}  "
          f"(>> 1e-3: adiabatic MISSES them)")
    print(f"           worst gap nearest-DRESSED-root   = {dressed_dbl_gap:.3e}  "
          f"(<1e-9: dressed INJECTS them)  -> {'PASS' if gate2 else 'FAIL'}")

    # ---- secular consistency: |fsec(omega)| (pole-cancelled) at each root ----
    detres = [abs(fsec(r, A0, omega_d, V)) for r in dressed]
    print(f"\n  secular residual |fsec(omega)| (pole-cancelled char. poly) at "
          f"dressed roots: max = {max(detres):.3e}")

    banner("VERDICT")
    all_ok = gate1 and gate2
    print(f"  GATE 1 (dressed == explicit-exact, all {ntot} roots, <1e-9): "
          f"{'PASS' if gate1 else 'FAIL'}  (max|err| {maxerr:.3e})")
    print(f"  GATE 2 (adiabatic misses the {Nd} 0OS doubles; dressed injects them): "
          f"{'PASS' if gate2 else 'FAIL'}")
    print(f"\n  OVERALL: {'PASS' if all_ok else 'FAIL'}")
    print("\nConclusion: the frequency-dependent dressed kernel A(omega)=A0+g_xc(omega),")
    print("solved by the secular-determinant pole search, reproduces the EXACT explicit")
    print("(single-spin-flip + 0OS-doubles) spectrum to machine precision and INJECTS the")
    print("0OS double-spin-flip states that the adiabatic kernel structurally cannot --")
    print("the QMRSF-DK live mechanism, validated independent of OpenQP/pyscf, with NO")
    print("determinant CI and NO external-Q PT2.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
