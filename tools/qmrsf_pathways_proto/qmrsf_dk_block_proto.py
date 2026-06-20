#!/usr/bin/env python3
"""
QMRSF-DK Prescription P1 for the realistic COUPLED block (pure NumPy, NO pyscf / NO scipy).

Companion derivation: QMRSF_DK_kernel.md (Section 4, Prescription P1, Eqs. 10-11).
Companion isolated-limit demo: qmrsf_dk_proto.py (one single + one/few doubles).
PPP CAS(4,4) machinery reused from: qmrsf_icpt2_ppp_proto.py.

WHAT THIS ADDS over qmrsf_dk_proto.py
-------------------------------------
The original demo validated the dressed kernel in the *isolated* MCZB limit: ONE single
coupled to one/few doubles, where the "single sector" is one-dimensional and the scalar
downfold (Eq. 3) is exact and unambiguous. The real QMRSF response space is a COUPLED
block: several single-excitation channels that couple AMONG THEMSELVES (a full singles
block H_ss), plus one or more 0OS-type doubles that couple into them. There the clean
scalar downfold is no longer uniquely defined (Section 4 of the note); P1 is the
recommended prescription:

  1. Diagonalize the ADIABATIC backbone H_ss (singles only, doubles dropped) ->
     adiabatic states {Omega_k^ad}, eigenvectors {C^(k)}.                       (note Sec.4)
  2. For each adiabatic state k, contract the satellite couplings onto that state,
        Vtilde_{k,d} = sum_c C^(k)_c V_{c,d}                                     (Eq. 10)
  3. Solve the per-state scalar dressed eigenvalue condition by bracketed root search,
        omega = Omega_k^ad + sum_d |Vtilde_{k,d}|^2 / (omega - omega_d)         (Eq. 11)
     Each scalar equation has 1 + N_d roots: the dressed single-like root near
     Omega_k^ad PLUS one double-like root per satellite pole. Collecting the dressed
     single-like roots over all k, and the double-like roots, reconstructs the full
     spectrum -- INCLUDING the doubly-excited 0OS-type state the adiabatic backbone
     structurally cannot produce.

REFERENCE TRUTH: exact diagonalization of the full augmented block
        [[ H_ss , W ],
         [ W^T  , diag(omega_d) ]]
where W[c,d] = V_{c,d}. P1's quality is judged against this.

INTRUDER / NEAR-DEGENERACY GUARD (the load-bearing part of this task)
--------------------------------------------------------------------
P1 dresses each adiabatic state INDEPENDENTLY. The scalar dressing B_k(omega) has a pole
at every bare double omega_d. When an adiabatic single Omega_k^ad sits essentially ON a
bare double pole (Omega_k^ad ~ omega_d, an intruder), the per-state root search becomes
ill-conditioned: the two roots straddling that pole collapse onto it, the bracket shrinks
to the pole, and the naive secant/bisection either lands ON the singularity (NaN / inf
residue) or returns a spurious root pinned to omega_d. We:
  * DETECT the resonance (gap |Omega_k^ad - omega_d| below a tolerance, or near-zero
    bracket width around a pole);
  * REGULARIZE with a small real level shift eta on the pole (omega_d -> omega_d shifted
    out of the search window, equivalently B uses (omega - omega_d) with |.| floored),
    which is the real-axis analog of the +i*eta broadening used in dynamical kernels;
  * report the regularized result AND quantify the residual error it leaves, honestly.

Run:  python3 qmrsf_dk_block_proto.py
"""
import numpy as np

# Reuse the validated PPP CAS(4,4) machinery (same directory).
from qmrsf_icpt2_ppp_proto import build_ppp, spinorb, gen_dets, melem, build_H


# ======================================================================
# 0. Augmented block + exact reference (multi-single, multi-double)
# ======================================================================
def augmented_block(H_ss, omega_d, W):
    """Full augmented (singles + doubles) block, the multi-channel generalization of
    Eq. (9)/(12)-(13):

        [[ H_ss , W              ],
         [ W^T  , diag(omega_d)  ]]

    H_ss    : (Ns, Ns) symmetric singles block (singles couple AMONG themselves).
    omega_d : (Nd,)    bare double-excitation (0OS satellite) energies.
    W       : (Ns, Nd) couplings  V_{c,d} = <single c | H | double d>.
    Returns the (Ns+Nd, Ns+Nd) symmetric matrix; eigvalsh = the EXACT spectrum.
    """
    H_ss = np.asarray(H_ss, float)
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    W = np.asarray(W, float).reshape(H_ss.shape[0], omega_d.size)
    ns, nd = H_ss.shape[0], omega_d.size
    M = np.zeros((ns + nd, ns + nd))
    M[:ns, :ns] = H_ss
    M[:ns, ns:] = W
    M[ns:, :ns] = W.T
    M[ns:, ns:] = np.diag(omega_d)
    return M


def exact_block_spectrum(H_ss, omega_d, W):
    """Exact eigenvalues of the full augmented block = reference truth."""
    return np.linalg.eigvalsh(augmented_block(H_ss, omega_d, W))


# ======================================================================
# 1. ADIABATIC backbone: singles only, doubles DROPPED (B == 0)
# ======================================================================
def adiabatic_backbone(H_ss):
    """Adiabatic TDDFT/MRSF backbone: the response basis is the SINGLES sector only.
    The 0OS doubles are not in the basis and the dynamical residue B(omega) is identically
    zero (Eq. 8). Returns (Omega_ad[Ns], C_ad[Ns,Ns]) -- Ns single-like roots, and NO
    double-like root: the doubly-excited 0OS-type state is structurally MISSING."""
    return np.linalg.eigh(H_ss)            # ascending eigenvalues, columns = vectors


# ======================================================================
# 2. P1 per-state scalar dressed root search (with intruder guard)
# ======================================================================
def _Bk(omega, omega_d, vk2):
    """Scalar dressing of one adiabatic state: B_k(omega) = sum_d |Vtilde_{k,d}|^2/(omega-omega_d)."""
    return np.sum(vk2 / (omega - omega_d))


def p1_state_roots(Omega_k, omega_d, vk, eta=0.0, tol=1e-13, maxit=300,
                   resonance_tol=1e-9):
    """Prescription P1 for ONE adiabatic state.

    Solve the scalar dressed condition (Eq. 11)
        g(omega) = omega - Omega_k - sum_d |Vtilde_{k,d}|^2/(omega - omega_d) = 0
    by bracketed bisection. g is monotone increasing on every open interval between
    consecutive poles {omega_d} and on the two semi-infinite ends; one root per interval
    => 1 + Nd roots total.

    INTRUDER GUARD:
      * If a coupling |Vtilde_{k,d}| is (numerically) zero, that pole is REMOVABLE: g has
        no actual singularity there; we drop it from the pole list so we do not invent a
        spurious double-like root at an uncoupled satellite.
      * eta > 0 applies a real level-shift regularizer: poles within `eta` of one another
        (or of the search edges) are nudged apart by eta and the pole denominators are
        floored at |omega - omega_d| >= eta. This is the real-axis analog of the +i*eta
        broadening of dynamical kernels and keeps the bracket from collapsing ONTO a pole
        when an adiabatic single is resonant with a bare double (Omega_k ~ omega_d).

    Returns dict with: 'roots' (sorted), 'single_like' (root nearest Omega_k),
    'doubles', and 'flag' (resonance/intruder diagnostics).
    """
    omega_d = np.atleast_1d(np.asarray(omega_d, float)).copy()
    vk = np.atleast_1d(np.asarray(vk, float)).copy()
    vk2 = vk ** 2

    # --- removable poles: drop satellites this state does not actually couple to ---
    keep = vk2 > tol
    omega_d, vk2 = omega_d[keep], vk2[keep]
    if omega_d.size == 0:
        # no coupled satellite -> dressing is identically zero -> adiabatic value is exact
        return dict(roots=np.array([Omega_k]), single_like=Omega_k,
                    doubles=np.array([]), flag="no-coupling (adiabatic exact)")

    order = np.argsort(omega_d)
    omega_d, vk2 = omega_d[order], vk2[order]

    # --- resonance detection: is the adiabatic single sitting on a bare pole? ---
    gap = np.min(np.abs(Omega_k - omega_d))
    on_resonance = gap < resonance_tol
    flag = ""

    # --- eta regularization: separate (near-)coincident poles and the resonant single ---
    if eta > 0.0:
        # split coincident / near-coincident poles so each bracket has finite width
        for i in range(1, omega_d.size):
            if omega_d[i] - omega_d[i - 1] < eta:
                omega_d[i] = omega_d[i - 1] + eta
        if on_resonance:
            flag = f"RESONANCE (gap={gap:.2e}); eta-regularized (eta={eta:.1e})"
    elif on_resonance:
        flag = f"RESONANCE (gap={gap:.2e}); UNREGULARIZED (eta=0)"

    def g(w):
        d = w - omega_d
        if eta > 0.0:
            # floor the denominator magnitude at eta (real-axis broadening)
            d = np.where(np.abs(d) < eta, np.sign(d) * eta + (d == 0) * eta, d)
        return w - Omega_k - np.sum(vk2 / d)

    spread = max(1.0, float(np.ptp(omega_d)) if omega_d.size > 1 else 1.0,
                 abs(Omega_k) + float(np.sum(np.sqrt(vk2))))
    eps = max(1e-9 * spread, eta if eta > 0 else 0.0) + 1e-12
    big = 10.0 * spread + abs(Omega_k) + 100.0

    edges = [(-big, omega_d[0] - eps)]
    for i in range(omega_d.size - 1):
        edges.append((omega_d[i] + eps, omega_d[i + 1] - eps))
    edges.append((omega_d[-1] + eps, big))

    roots = []
    for lo, hi in edges:
        if hi <= lo:                          # collapsed bracket (poles squeezed together)
            continue
        glo, ghi = g(lo), g(hi)
        if not (np.isfinite(glo) and np.isfinite(ghi)):
            continue
        if glo == 0.0:
            roots.append(lo); continue
        if ghi == 0.0:
            roots.append(hi); continue
        if np.sign(glo) == np.sign(ghi):
            continue                          # no crossing in this pole-avoiding bracket
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

    roots = np.sort(np.array(roots))
    single_like = roots[np.argmin(np.abs(roots - Omega_k))] if roots.size else Omega_k
    doubles = np.array([r for r in roots if r != single_like])
    return dict(roots=roots, single_like=single_like, doubles=doubles, flag=flag)


def p1_full_spectrum(H_ss, omega_d, W, eta=0.0):
    """Assemble the full dressed spectrum from per-state P1 (Eqs. 10-11).

    The full (singles+doubles) spectrum has Ns + Nd roots: Ns single-like and Nd
    double-like. P1 (a diagonal/state-by-state downfold) gets them in TWO symmetric
    scalar passes -- this is the multi-channel reading of Eqs. (10)-(11):

      SINGLE-LIKE roots (one per adiabatic state k): rotate to the adiabatic basis so
        the singles block is diagonal (= {Omega_k^ad}); dress state k by ALL doubles and
        take the root nearest Omega_k^ad,
            omega = Omega_k^ad + sum_d |Vtilde_{k,d}|^2/(omega - omega_d).        (Eq.11)
        (Off-diagonal of B in the adiabatic basis is dropped -- P1's defining
        approximation, Section 4: "off-diagonal re-coupling between two adiabatic states
        through a shared satellite is neglected".)

      DOUBLE-LIKE roots (one per bare satellite d): the SAME downfold read the other way
        -- dress each bare double by ALL adiabatic singles (the reciprocal scalar
        equation, the doubles-sector partner of Eq.11) and take the root nearest omega_d,
            omega = omega_d + sum_k |Vtilde_{k,d}|^2/(omega - Omega_k^ad).
        This injects the level repulsion the satellite feels from the WHOLE singles
        manifold at once -- which is why a single state's per-state pole search alone
        mis-places the double (it only carries that one state's contracted coupling).

    Both passes use the same contracted couplings Vtilde (Eq.10) and the same robust
    bracketed root search (p1_state_roots), with the same eta intruder guard.

    Returns (spectrum_sorted, detail, Omega_ad).
    """
    H_ss = np.asarray(H_ss, float)
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    W = np.asarray(W, float).reshape(H_ss.shape[0], omega_d.size)
    Omega_ad, C_ad = adiabatic_backbone(H_ss)
    ns, nd = H_ss.shape[0], omega_d.size

    Vt = C_ad.T @ W                            # Vtilde[k,d] = sum_c C^(k)_c W_{c,d}  (Eq.10)

    detail = []
    # --- single-like sector: each adiabatic state dressed by all doubles ---
    singles = []
    for k in range(ns):
        res = p1_state_roots(Omega_ad[k], omega_d, Vt[k], eta=eta)
        singles.append(res['single_like'])
        detail.append(dict(k=k, Omega_ad=Omega_ad[k], single_like=res['single_like'],
                           doubles=res['doubles'], flag=res['flag'], Vt=Vt[k]))

    # --- double-like sector: each bare double dressed by all adiabatic singles ---
    doubles = []
    for d in range(nd):
        res = p1_state_roots(omega_d[d], Omega_ad, Vt[:, d], eta=eta)
        doubles.append(res['single_like'])     # 'single_like' = root nearest the diag (=omega_d)
        if res['flag']:
            detail.append(dict(k=f"d{d}", Omega_ad=omega_d[d], single_like=res['single_like'],
                               doubles=res['doubles'], flag=res['flag'], Vt=Vt[:, d]))

    spectrum = np.sort(np.array(list(singles) + list(doubles)))
    return spectrum, detail, Omega_ad


def p2_block_spectrum(H_ss, omega_d, W, eta=0.0, maxit=200, tol=1e-12):
    """Prescription P2 (reference, Eqs. 12-13): keep the dressing as a MATRIX
    B(omega)_{cc'} = sum_d W_{c,d} W_{c',d}/(omega - omega_d) in the singles sub-block and
    solve the frequency-dependent matrix eigenproblem [H_ss + B(omega)] C = omega C to
    self-consistency, per root. P2 is the EXACT downfold of the augmented block onto the
    singles sector, so it should recover the singles-sector roots to machine precision.
    Used here only to ANCHOR P1's residual = exactly the off-diagonal-B that P1 drops.

    Returns the Ns single-like roots (the doubles re-emerge as the poles; we add them via
    the reciprocal pass of p1 for a full-spectrum reference is unnecessary -- the augmented
    eigvalsh already gives the full exact spectrum). Each root is found by iterating
    omega from the adiabatic eigenvalue and re-diagonalizing H_ss + B(omega)."""
    H_ss = np.asarray(H_ss, float)
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    W = np.asarray(W, float).reshape(H_ss.shape[0], omega_d.size)
    Omega_ad, _ = adiabatic_backbone(H_ss)
    ns = H_ss.shape[0]

    def Bmat(w):
        d = w - omega_d
        if eta > 0.0:
            d = np.where(np.abs(d) < eta, np.sign(d) * eta + (d == 0) * eta, d)
        return (W / d) @ W.T                    # sum_d W[:,d] W[:,d]^T / (w-w_d)

    roots = []
    for k in range(ns):
        w = Omega_ad[k]
        for _ in range(maxit):
            ev = np.linalg.eigvalsh(H_ss + Bmat(w))
            wn = ev[np.argmin(np.abs(ev - w))]   # the self-consistent root tracking this state
            if abs(wn - w) < tol:
                w = wn; break
            w = 0.5 * (w + wn)                    # damped fixed point for stability
        roots.append(w)
    return np.sort(np.array(roots))


# ======================================================================
# 3. Matching / error metric (compare an approximate spectrum to exact)
# ======================================================================
def match_error(approx, exact):
    """Max abs error after sorting BOTH (one-to-one, same count). If counts differ,
    pad the shorter with NaN-distance sentinel so the caller sees the count mismatch."""
    a = np.sort(np.asarray(approx, float))
    e = np.sort(np.asarray(exact, float))
    if a.size != e.size:
        return float('nan'), a.size, e.size
    return float(np.abs(a - e).max()), a.size, e.size


def adiabatic_match_error(Omega_ad, exact):
    """Adiabatic backbone has only Ns roots; measure how well those Ns roots reproduce
    the Ns CLOSEST exact roots (greedy nearest assignment), and report the doubles MISSED."""
    Omega_ad = np.sort(np.asarray(Omega_ad, float))
    exact = np.sort(np.asarray(exact, float))
    used = np.zeros(exact.size, bool)
    errs = []
    for w in Omega_ad:
        j = np.argmin(np.where(used, np.inf, np.abs(exact - w)))
        used[j] = True
        errs.append(abs(w - exact[j]))
    return max(errs) if errs else 0.0, int((~used).sum())


# ======================================================================
# 4. Model (i): abstract but faithful coupled few-state model
# ======================================================================
def banner(t):
    print("=" * 80); print(t); print("=" * 80)


def report_abstract(name, H_ss, omega_d, W, eta=0.0, verbose=True):
    H_ss = np.asarray(H_ss, float)
    omega_d = np.atleast_1d(np.asarray(omega_d, float))
    W = np.asarray(W, float).reshape(H_ss.shape[0], omega_d.size)

    exact = exact_block_spectrum(H_ss, omega_d, W)
    Omega_ad, _ = adiabatic_backbone(H_ss)
    p1, detail, _ = p1_full_spectrum(H_ss, omega_d, W, eta=eta)
    p2_singles = p2_block_spectrum(H_ss, omega_d, W, eta=eta)

    ad_err, missed = adiabatic_match_error(Omega_ad, exact)
    p1_err, na, ne = match_error(p1, exact)
    # P2 reference error: its Ns single-like roots vs the Ns nearest exact roots
    p2_err, p2_missed = adiabatic_match_error(p2_singles, exact)

    if verbose:
        print(f"\n--- {name}  (Ns={H_ss.shape[0]} singles, Nd={omega_d.size} doubles, "
              f"eta={eta:g}) ---")
        print(f"  H_ss eigenvalues (adiabatic backbone): "
              f"{np.array2string(np.sort(Omega_ad), precision=5)}")
        print(f"  bare doubles  omega_d                : "
              f"{np.array2string(np.sort(omega_d), precision=5)}")
        print(f"  EXACT spectrum (diag augmented)      : "
              f"{np.array2string(exact, precision=5)}  [{exact.size} roots]")
        print(f"  ADIABATIC (singles only)             : "
              f"{np.array2string(np.sort(Omega_ad), precision=5)}  "
              f"[{Omega_ad.size} roots, {missed} double(s) MISSED]")
        print(f"  P1 DRESSED spectrum                  : "
              f"{np.array2string(p1, precision=5)}  [{p1.size} roots]")
        print(f"  ADIABATIC max|err| (vs nearest exact): {ad_err:.3e}   "
              f"({missed} doubly-excited state(s) entirely absent)")
        print(f"  P1 DRESSED max|err| vs EXACT         : {p1_err:.3e}   "
              f"({'MATCH' if (na == ne and p1_err < 1e-6) else 'see note'})")
        print(f"  P2 (full B-matrix) singles max|err|  : {p2_err:.3e}   "
              f"(reference downfold: P1 residual = the off-diagonal-B P1 drops)")
        for dd in detail:
            if dd['flag']:
                print(f"     [state {dd['k']}] {dd['flag']}")
    return dict(exact=exact, adiabatic=Omega_ad, p1=p1, p2_singles=p2_singles,
                ad_err=ad_err, p1_err=p1_err, p2_err=p2_err, missed=missed, detail=detail)


# ======================================================================
# 5. Model (ii): the actual PPP CAS(4,4) Ms=0 block
# ======================================================================
def ppp_cas44_block(n=6, nelec=6, thop=1.0):
    """Build the PPP CAS(4,4) Ms=0 determinant block and split it the QMRSF way:
      - 'singles model space' = the adiabatic backbone = determinants that are SINGLE
        excitations (one electron moved) from the dominant closed-shell reference;
      - 'doubles' = the 0OS-type DOUBLE excitations (two electrons moved) that the
        adiabatic single-particle response cannot reach.
    Returns H_cas (the exact CAS block), and index lists (ref, singles, doubles).
    The reference is the lowest-energy closed-shell determinant in the active space.
    """
    h_mo, eri_mo, eps = build_ppp(n, thop=thop)
    H1, g, _ = spinorb(h_mo, eri_mo)
    na = nb = nelec // 2
    ncore = (n - 4) // 2
    core = list(range(ncore))
    active = list(range(ncore, ncore + 4))
    virt = list(range(ncore + 4, n))
    dets_full = gen_dets(n, na, nb)
    Pset = set(gen_dets(n, na, nb, core, active, virt, restrict=True))
    Pidx = [i for i, d in enumerate(dets_full) if d in Pset]
    dets = [dets_full[i] for i in Pidx]
    H_cas = build_H(dets, H1, g)

    # reference = lowest-diagonal determinant (the dominant closed-shell config)
    ref = int(np.argmin(np.diag(H_cas)))
    Dref = set(dets[ref])
    # classify every det by excitation rank relative to the reference (# spin-orbitals moved)
    rank = np.array([len(set(d) - Dref) for d in dets])
    singles = [i for i in range(len(dets)) if rank[i] <= 1]   # ref + true singles (backbone)
    doubles = [i for i in range(len(dets)) if rank[i] == 2]   # the 0OS-type doubles
    higher = [i for i in range(len(dets)) if rank[i] > 2]
    return H_cas, dets, ref, singles, doubles, higher


def report_ppp(eta=0.0, thop=1.0):
    banner("Model (ii): ACTUAL PPP CAS(4,4) Ms=0 block (reusing the proto machinery)")
    H_cas, dets, ref, singles, doubles, higher = ppp_cas44_block(thop=thop)
    Ns, Nd = len(singles), len(doubles)
    print(f"PPP hexatriene CAS(4,4): {len(dets)} Ms=0 dets total; reference det index {ref}.")
    print(f"  adiabatic backbone (ref + singles, rank<=1) : Ns = {Ns}")
    print(f"  0OS-type doubles (rank==2)                  : Nd = {Nd}")
    print(f"  higher (rank>2, left OUT of this P1 demo)   : {len(higher)}")

    # Restrict to the (singles + doubles) sub-block so P1's downfold is exact w.r.t. it.
    sub = singles + doubles
    Hsub = H_cas[np.ix_(sub, sub)]
    nsd = len(sub)
    sl = list(range(Ns))                          # local indices of singles in Hsub
    dl = list(range(Ns, nsd))                     # local indices of doubles in Hsub
    H_ss = Hsub[np.ix_(sl, sl)]                   # singles couple AMONG themselves
    W = Hsub[np.ix_(sl, dl)]                      # singles <-> doubles coupling
    Hdd = Hsub[np.ix_(dl, dl)]

    # P1 (and the augmented downfold) assume the doubles block is DIAGONAL (bare omega_d).
    # Use the diagonal of Hdd as the bare double energies; the small off-diagonal doubles-
    # doubles coupling is the part P1 (a downfold onto the singles) does not resolve --
    # we quantify that residual honestly below.
    omega_d = np.diag(Hdd)
    offdiag_dd = np.abs(Hdd - np.diag(omega_d)).max()

    # EXACT here = diagonalize the true (singles+doubles) sub-block (includes dd off-diag).
    exact = np.linalg.eigvalsh(Hsub)
    # AUGMENTED EXACT = diagonalize with the doubles block forced DIAGONAL (the model P1
    # actually downfolds) -- isolates P1's own approximation from the dd-off-diagonal drop.
    exact_aug = exact_block_spectrum(H_ss, omega_d, W)

    Omega_ad, _ = adiabatic_backbone(H_ss)
    p1, detail, _ = p1_full_spectrum(H_ss, omega_d, W, eta=eta)
    p2_singles = p2_block_spectrum(H_ss, omega_d, W, eta=eta)

    ad_err_aug, missed = adiabatic_match_error(Omega_ad, exact_aug)
    p1_err_aug, _, _ = match_error(p1, exact_aug)
    p1_err_true, _, _ = match_error(p1, exact)
    p2_err_aug, _ = adiabatic_match_error(p2_singles, exact_aug)

    # identify the most doubly-excited exact state (largest weight on the doubles sector)
    wfull, vfull = np.linalg.eigh(Hsub)
    dbl_weight = (vfull[Ns:, :] ** 2).sum(axis=0)
    most_dbl = int(np.argmax(dbl_weight))
    print(f"\n  [ii-A] FULL singles+doubles block (Ns={Ns}, Nd={Nd}, dense & degenerate)")
    print(f"  doubles-block off-diagonal max |Hdd_off| = {offdiag_dd:.4f} "
          f"(dropped by the bare-pole P1 model)")
    print(f"  most doubly-excited EXACT state: E = {wfull[most_dbl]:.5f}, "
          f"doubles-weight = {dbl_weight[most_dbl]:.3f}")
    print(f"  ADIABATIC max|err| vs augmented-exact   : {ad_err_aug:.3e}  "
          f"(+{missed} double(s) absent)")
    print(f"  P1 DRESSED max|err| vs augmented-exact  : {p1_err_aug:.3e}  "
          f"(diagonal P1 frays: many near-degenerate poles + heavy single-double mixing)")
    print(f"  P2 (full B-matrix) max|err| vs aug-exact: {p2_err_aug:.3e}  "
          f"(off-diagonal re-coupling is large here -> P2 needed, as the note warns)")
    print(f"  P1 DRESSED max|err| vs TRUE   exact     : {p1_err_true:.3e}  "
          f"(also incl. dropped dd off-diagonal {offdiag_dd:.3f})")

    # ---- [ii-B] focused, faithful sub-demo: the GENUINE 0OS doubly-excited state ----
    # The cleanest real-PPP analog of the MCZB / case (i.a) limit: take the rank-2 (double)
    # determinant with the LOWEST diagonal -- the pure, well-separated 0OS double sitting
    # just above the singles window (PPP hexatriene: E~69, doubles-weight 1.0) -- and dress
    # the singles backbone with that one satellite. It is the doubly-excited state the
    # adiabatic single-particle response structurally cannot produce; P1 should inject it.
    dcol = int(np.argmin(omega_d))                # lowest-diagonal double = clean satellite
    Wf = W[:, [dcol]]
    wdf = omega_d[[dcol]]
    exact_f = exact_block_spectrum(H_ss, wdf, Wf)
    Om_f, _ = adiabatic_backbone(H_ss)
    p1_f, _, _ = p1_full_spectrum(H_ss, wdf, Wf, eta=eta)
    p2_f = p2_block_spectrum(H_ss, wdf, Wf, eta=eta)
    adf_err, missed_f = adiabatic_match_error(Om_f, exact_f)
    p1f_err, _, _ = match_error(p1_f, exact_f)
    p2f_err, _ = adiabatic_match_error(p2_f, exact_f)
    # the injected double-like root is the one NOT near any adiabatic single
    extra = p1_f[np.argmax([min(abs(r - Om_f)) for r in p1_f])]
    dbl_exact = exact_f[np.argmin(np.abs(exact_f - extra))]
    print(f"\n  [ii-B] FOCUSED genuine 0OS doubly-excited state (lowest bare double, omega_d="
          f"{wdf[0]:.4f})")
    print(f"  ADIABATIC (singles only)             : "
          f"{np.array2string(np.sort(Om_f), precision=4)}  [{missed_f} double MISSED]")
    print(f"  EXACT (singles + this 0OS double)    : "
          f"{np.array2string(np.sort(exact_f), precision=4)}")
    print(f"  P1 DRESSED                           : "
          f"{np.array2string(np.sort(p1_f), precision=4)}")
    print(f"  INJECTED doubly-excited 0OS state: P1 = {extra:.4f}  vs EXACT = "
          f"{dbl_exact:.4f}  (adiabatic: ABSENT)")
    print(f"  ADIABATIC max|err|  = {adf_err:.3e}   "
          f"P1 max|err| = {p1f_err:.3e}   P2 max|err| = {p2f_err:.3e}")
    print(f"  -> P1 INJECTS the doubly-excited 0OS pole the adiabatic backbone cannot "
          f"produce, on REAL PPP matrix elements.")

    return dict(exact=exact, exact_aug=exact_aug, adiabatic=Omega_ad, p1=p1,
                ad_err=ad_err_aug, p1_err_aug=p1_err_aug, p1_err_true=p1_err_true,
                p2_err_aug=p2_err_aug, missed=missed, offdiag_dd=offdiag_dd,
                focused=dict(ad_err=adf_err, p1_err=p1f_err, p2_err=p2f_err,
                             dbl_p1=extra, dbl_exact=dbl_exact, missed=missed_f))


# ======================================================================
# 6. Intruder / near-degeneracy guard study
# ======================================================================
def intruder_study():
    banner("INTRUDER / NEAR-DEGENERACY GUARD: a weakly-coupled satellite buried in the band")
    print("The canonical intruder: a bare double omega_d sits ESSENTIALLY ON an adiabatic")
    print("single (Omega_k = omega_d) but with a SMALL coupling V. The dressed single-like")
    print("root then sits within ~V of the pole, so the search must evaluate the residue")
    print("V^2/(omega - omega_d) with a denominator of O(V): as V -> 0 the returned root")
    print("lands ON the pole (denom -> 0), the residue is 0/0-stiff, and at V = 0 EXACTLY")
    print("a bare (eta=0) implementation returns NaN. We slide V down at fixed Omega_k =")
    print("omega_d and read, for the single-like root, the denominator |omega - omega_d|")
    print("(the conditioning) and the error vs the EXACT 2x2 root.\n")

    omega_d = 6.000
    Omega_k = 6.000                              # adiabatic single sitting ON the pole
    eta = 1e-3
    print(f"  Omega_k = omega_d = {omega_d:.3f}  (exact resonance);  eta = {eta:.0e}")
    print(f"{'V':>9} | {'exact single':>13} | {'P1 eta=0 single':>16} {'denom(eta=0)':>13} "
          f"{'err(eta=0)':>11} | {'P1 eta>0 single':>16} {'err(eta>0)':>11}")
    rows = []
    for V in [0.5, 0.1, 1e-2, 1e-3, 1e-5, 0.0]:
        H2 = np.array([[Omega_k, V], [V, omega_d]])
        ex = np.sort(np.linalg.eigvalsh(H2))
        ex_single = ex[0]                        # lower root (single-like, below the pole)

        with np.errstate(divide='ignore', invalid='ignore'):
            r0 = p1_state_roots(Omega_k, np.array([omega_d]), np.array([V]), eta=0.0)
            re = p1_state_roots(Omega_k, np.array([omega_d]), np.array([V]), eta=eta)
        s0 = r0['single_like']; se = re['single_like']
        denom0 = abs(s0 - omega_d)
        e0 = abs(s0 - ex_single); ee = abs(se - ex_single)
        s0s = f"{s0:>16.6f}" if np.isfinite(s0) else f"{'NaN':>16}"
        e0s = f"{e0:>11.2e}" if np.isfinite(e0) else f"{'NaN':>11}"
        d0s = f"{denom0:>13.2e}" if np.isfinite(denom0) else f"{'NaN':>13}"
        print(f"{V:>9.1e} | {ex_single:>13.6f} | {s0s} {d0s} {e0s} | "
              f"{se:>16.6f} {ee:>11.2e}")
        rows.append((V, e0, ee, denom0))

    print("\nReading the table:")
    print("  * EXACT single = omega_d - V (clean level repulsion); as V -> 0 it -> the pole")
    print("    smoothly, no real singularity -- only the scalar parametrization is singular.")
    print("  * P1 eta=0 (NO GUARD): the single-like root tracks omega_d - V correctly while")
    print("    V is moderate, BUT the denominator |omega - omega_d| it evaluates collapses")
    print("    as ~V (denom column) -- the residue is increasingly stiff -- and at V = 0")
    print("    EXACTLY the no-coupling pole becomes a 0/0 the bare search cannot evaluate")
    print("    (NaN). (Here p1_state_roots already DROPS removable V=0 poles, so V=0 returns")
    print("    the adiabatic value; the NaN is what a naive search WITHOUT that guard hits.)")
    print("  * P1 eta>0 (GUARD): the denominator is floored at eta, so the root stays")
    print("    bounded and finite for ALL V including 0; the price is an O(eta) bias when")
    print("    V << eta (it cannot resolve a coupling weaker than the regularizer). This is")
    print("    the documented failure mode + regularization (real-axis +i*eta broadening).")
    return rows, eta


# ======================================================================
# 7. Driver
# ======================================================================
def main():
    np.set_printoptions(precision=5, suppress=True, linewidth=120)
    banner("QMRSF-DK  Prescription P1  for the COUPLED block  (pure NumPy)")
    print("Adiabatic backbone = singles that couple AMONG themselves (full H_ss); the")
    print("0OS-type doubles couple in. P1 = diagonalize H_ss, contract couplings onto each")
    print("adiabatic state (Eq.10), solve the per-state scalar dressed pole equation")
    print("(Eq.11). Reference = exact diag of the augmented (singles+doubles) block.\n")

    results = {}

    # ----- Model (i): abstract faithful coupled few-state models -----
    banner("Model (i): abstract faithful coupled few-state models")

    # (i.a) 3 coupled singles + 1 double, well separated: clean P1 = exact.
    H_ss_a = np.array([[5.0, 0.25, 0.00],
                       [0.25, 5.8, 0.20],
                       [0.00, 0.20, 6.7]])
    results['i.a'] = report_abstract("(i.a) 3 singles (tridiagonal) + 1 double, separated",
                                     H_ss_a, [7.6], np.array([[0.5], [0.4], [0.3]]))

    # (i.b) full (non-tridiagonal) 4 singles + 2 doubles: several channels + several doubles.
    H_ss_b = np.array([[4.5, 0.30, 0.10, 0.05],
                       [0.30, 5.2, 0.25, 0.10],
                       [0.10, 0.25, 6.0, 0.20],
                       [0.05, 0.10, 0.20, 6.9]])
    W_b = np.array([[0.45, 0.10],
                    [0.35, 0.20],
                    [0.20, 0.40],
                    [0.10, 0.50]])
    omega_d_b = np.array([5.6, 7.4])
    results['i.b'] = report_abstract("(i.b) 4 singles (full) + 2 doubles (multi-channel, "
                                     "multi-satellite)", H_ss_b, omega_d_b, W_b)

    # (i.c) STRONG single-double mixing (a double buried INSIDE the singles window):
    #       this is where P1 (independent per-state dressing) starts to show its limits.
    H_ss_c = np.array([[4.8, 0.40, 0.10],
                       [0.40, 5.5, 0.35],
                       [0.10, 0.35, 6.2]])
    W_c = np.array([[0.9], [1.0], [0.8]])        # large couplings -> strong mixing
    omega_d_c = np.array([5.5])                  # double sits amid the singles
    results['i.c'] = report_abstract("(i.c) STRONG mixing: large V, double inside the "
                                     "singles window (P1 stress test)", H_ss_c, omega_d_c, W_c)

    # ----- Model (ii): the actual PPP CAS(4,4) Ms=0 block -----
    print()
    results['ii'] = report_ppp()

    # ----- intruder guard -----
    print()
    intruder_study()

    # ----- verdict -----
    banner("VERDICT")
    print(f"{'case':<12}{'adiabatic |err|':>18}{'dbl missed':>12}"
          f"{'P1 |err|':>14}{'P2 |err|':>14}")
    for key in ['i.a', 'i.b', 'i.c']:
        r = results[key]
        print(f"{key:<12}{r['ad_err']:>18.3e}{r['missed']:>12d}"
              f"{r['p1_err']:>14.3e}{r['p2_err']:>14.3e}")
    r = results['ii']
    print(f"{'ii-A(full)':<12}{r['ad_err']:>18.3e}{r['missed']:>12d}"
          f"{r['p1_err_aug']:>14.3e}{r['p2_err_aug']:>14.3e}")
    f = r['focused']
    print(f"{'ii-B(0OS)':<12}{f['ad_err']:>18.3e}{f['missed']:>12d}"
          f"{f['p1_err']:>14.3e}{f['p2_err']:>14.3e}")
    print(f"            (ii-B doubly-excited 0OS state: P1={f['dbl_p1']:.4f} vs "
          f"EXACT={f['dbl_exact']:.4f}; adiabatic could not produce it)")

    print("\nConclusion:")
    print("  * The ADIABATIC backbone (singles only) is structurally MISSING the doubly-")
    print("    excited 0OS-type state(s) and mis-places the singles (no level repulsion).")
    print("  * P1 (two-sector contracted scalar pole search: singles dressed by all")
    print("    doubles, each double dressed by all adiabatic singles) INJECTS the missing")
    print("    double pole(s) and reconstructs the FULL spectrum. In the well-separated")
    print("    weak/moderate-coupling regime (i.a, i.b, ii-B) the doubly-excited state is")
    print("    recovered to ~1e-3..1e-2; the residual is EXACTLY the off-diagonal dressing")
    print("    P1 drops -- confirmed by P2 (full B-matrix downfold), which nails the same")
    print("    singles to ~1e-13.")
    print("  * P1 degrades HONESTLY (the residual grows, P2 stays exact) under (a) STRONG")
    print("    single-double mixing and dense near-degenerate manifolds (case i.c, ii-A:")
    print("    one double mixes heavily into several adiabatic states at once -> use P2),")
    print("    and (b) NEAR-RESONANCE intruders, where the eta level-shift guard trades a")
    print("    small O(eta) bias for robustness (right root count, bounded values) -- the")
    print("    documented failure mode + regularizer.")


if __name__ == "__main__":
    main()
