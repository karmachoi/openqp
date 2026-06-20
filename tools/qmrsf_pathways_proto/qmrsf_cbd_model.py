#!/usr/bin/env python3
"""
QMRSF cyclobutadiene (CBD) model -- the manuscript's smoking gun, in a
controllable pure-NumPy CAS(4,4) toy (NO pyscf, NO scipy).

WHAT THIS ADDS over the hexatriene proto / multistate file
----------------------------------------------------------
The hexatriene toy (qmrsf_icpt2_multistate.py) could NOT induce a hard
discontinuity in the naive glued/patched QMRSF scheme: its ground state stays
robustly closed-shell-dominated, the two leading diabats run nearly parallel,
and the failure of the naive hybrid is SYSTEMATIC (a wrong barrier / mis-tuned
gap) rather than a kink.  The honest note at the bottom of that file says so.

The reason, physically, is that hexatriene's ground state never CHANGES
CHARACTER.  Cyclobutadiene does.  This file builds the regime that breaks a
glued scheme:

  *  4 pi centers on a RING (4-site cyclic PPP/Hubbard-Ohno), 4 electrons,
     CAS(4,4) = the FULL pi space (all four pi orbitals active, no pi core).
  *  A bond-length-alternation parameter `delta` interpolates
       - square  (delta=0, D4h): cyclic hopping uniform -> Hueckel levels
         -2t, 0, 0, +2t.  The two frontier orbitals are EXACTLY DEGENERATE
         -> strong static correlation, the COVALENT open-shell (4OS-type)
         singlet biradical is the ground state (antiaromatic D4h point).
       - rectangular (delta>0, D2h): alternating hopping t*(1+/-delta) splits
         the frontier pair -> a CLOSED-SHELL (0OS-type) IONIC determinant drops
         and becomes the ground state.
     So the ground state's character CROSSES from covalent/open-shell (4OS) at
     and below the square point to ionic/closed-shell (0OS) on the rectangular
     side -- exactly the 4OS<->0OS crossover the manuscript says a glued QMRSF
     patches per-system and therefore renders DISCONTINUOUS.

  *  A handful of higher "sigma-like" external virtual orbitals are appended so
     the icPT2 external-Q space is non-empty (a bare 4-orbital pi model has no
     perturbers).  The pi CAS(4,4) is unchanged by them at zeroth order; they
     only supply the uniform dynamic-correlation layer.

WHAT WE COMPUTE (the three curves the manuscript contrasts)
-----------------------------------------------------------
For the lowest 2-3 states, scanned in `delta` THROUGH the square point:
  (a) EXACT FCI                              -- the truth.
  (b) ONE-OPERATOR multistate icPT2 downfold -- the determinant-union backbone
        (one consistent CAS-CI operator) + ONE uniform Hermitian des-Cloizeaux
        external-Q self-energy.  Should track FCI SMOOTHLY (continuous, correct
        ordering, accurate automerization barrier).
  (c) The NAIVE GLUED HYBRID that mimics the original QMRSF defect: the 1SF
        (2OS/4OS) diagonals are DFT-dressed/shifted, but the closed-shell 0OS
        diagonal is referenced INDEPENDENTLY with a per-system constant shift
        fit at ONE geometry (the rectangular end, where 0OS is the ground
        state -- exactly how a real code would calibrate it).  As the 0OS
        character crosses IN toward the square point, the frozen, separately-
        referenced 0OS block sits at the wrong height relative to the dressed
        1SF block -> a KINK / wrong-sign automerization energy.

QUANTIFIED
----------
  * Automerization barrier analogue  E(square) - E(rectangular minimum)
    for FCI vs one-operator vs naive hybrid (sign + magnitude).
  * max|2nd difference| (smoothness) of each S0 curve over the fine scan.
  * whether the naive hybrid shows a kink / wrong-order / wrong-sign result.

Run:  python3 qmrsf_cbd_model.py
"""
import os
import numpy as np
from itertools import combinations

# Reuse the validated Slater-Condon / determinant-CI core verbatim.
import qmrsf_icpt2_ppp_proto as proto
from qmrsf_icpt2_ppp_proto import spinorb, gen_dets, build_H


# ======================================================================
# 1. 4-site cyclic PPP/Hubbard-Ohno model with bond-length alternation
# ======================================================================
def build_cbd(t=2.40, U=8.0, delta=0.0, n_ext=2, ext_gap=6.0, R=1.40,
              kappa=1.0):
    """Cyclobutadiene-like 4-center pi ring, CAS(4,4) = full pi space, with a
    few external virtual orbitals appended for the dynamic-correlation layer.

    Ring sites 0-1-2-3-0.  Bond-length alternation `delta` modulates the four
    ring bonds as  t*(1+delta), t*(1-delta), t*(1+delta), t*(1-delta)  so that
        delta = 0  -> square D4h  (all bonds equal; degenerate frontier pair)
        delta > 0  -> rectangular D2h (alternating short/long bonds; split pair)
    Hopping is the off-diagonal Hueckel element  h = -t_bond.

    Ohno 2e (ZDO): on-site U on the diagonal, screened 1/r off-site using ring
    geometry (square of side R; pos of the 4 sites placed on a square so that
    1,3 are para (distance R*sqrt2) and 0,2 are para).  The screened-Coulomb
    form matches proto.build_ppp (gamma_ij = 14.397/sqrt((14.397/U)^2 + r^2)).

    External orbitals: `n_ext` higher virtuals at ~ +ext_gap above the pi
    manifold, weakly Coulomb-coupled (a compact "sigma*" set).  They give the
    icPT2 a real external-Q space without touching the pi CAS at zeroth order.

    Returns h_mo (m,m), eri_mo (m,m,m,m) chemist (pq|rs), eps, where m = 4+n_ext.
    The first 4 MOs are the pi space (the CAS active set); the rest are external.
    """
    npi = 4
    # ----- ring geometry: a square of side R -----
    # place the four centers on a square; bond-alternation only enters hopping.
    sq = np.array([[0.0, 0.0],
                   [R,   0.0],
                   [R,   R  ],
                   [0.0, R  ]])
    # ----- one-electron (Hueckel) ring matrix with alternating hopping -----
    h_site = np.zeros((npi, npi))
    tb = [t * (1.0 + delta), t * (1.0 - delta),
          t * (1.0 + delta), t * (1.0 - delta)]            # bonds 0-1,1-2,2-3,3-0
    bonds = [(0, 1), (1, 2), (2, 3), (3, 0)]
    for (i, j), tij in zip(bonds, tb):
        h_site[i, j] = h_site[j, i] = -tij
    eps_pi, Cpi = np.linalg.eigh(h_site)                    # pi MOs ascending

    # ----- Ohno gamma over the ring sites (U=0 -> no 2e at all) -----
    gam = np.zeros((npi, npi))
    if U > 0.0:
        for i in range(npi):
            for j in range(npi):
                if i == j:
                    gam[i, j] = U
                else:
                    r = np.linalg.norm(sq[i] - sq[j])
                    gam[i, j] = 14.397 / np.sqrt((14.397 / U) ** 2 + r ** 2)

    # pi-block MO integrals
    if n_ext == 0:
        h_pi = Cpi.T @ h_site @ Cpi
        Mpi = np.einsum('ip,iq->pqi', Cpi, Cpi)
        eri_pi = np.einsum('pqi,ij,rsj->pqrs', Mpi, gam, Mpi)
        return h_pi, eri_pi, eps_pi

    # ----- append n_ext external "auxiliary" virtual orbitals -----
    # The externals are high-lying virtuals (a compact sigma*-like set) that are
    # DECOUPLED from the ring in the ONE-electron Hamiltonian (block diagonal, so
    # they never mix into the pi CAS at zeroth order) but coupled to the pi space
    # in the TWO-electron part by a weak transfer integral (below).  That weak
    # transfer makes the external Q space a genuine, smooth dynamic-correlation
    # channel for the pi CAS -- a pure block-diagonal ZDO would give only Coulomb
    # (pq|ab) and no transfer (pa|qr), leaving CAS==FCI (verified).
    m = npi + n_ext

    # ----- ONE-electron part (block diagonal: aux are pure external virtuals) --
    # The aux on-site energy must clear not just the top pi level but also the
    # large pi-pi Ohno Coulomb repulsion an added electron would feel in the pi
    # manifold (~2U here).  Otherwise an "external" orbital at +ext_gap above the
    # bare Hueckel LUMO is actually BELOW a doubly-occupied pi orbital once the
    # ~2U Coulomb is included, and the FCI collapses into aux occupation
    # (CAS weight -> 0).  We place the aux at eps_pi[-1] + 2U + ext_gap.
    h_mo = np.zeros((m, m))
    h_mo[:npi, :npi] = Cpi.T @ h_site @ Cpi                  # pi MO 1e (block-diag)
    base = eps_pi[-1] + 2.0 * U + ext_gap
    eps_aux = np.array([base + 2.0 * a for a in range(n_ext)])
    for a in range(n_ext):
        h_mo[npi + a, npi + a] = eps_aux[a]                  # high-lying virtuals

    # ----- TWO-electron part ------------------------------------------------
    # (i)  pi-pi-pi-pi block: the EXACT ZDO Ohno tensor (reproduced bit-for-bit).
    # (ii) a WEAK, uniform pi<->aux EXCHANGE/transfer coupling that drives modest
    #      pi->aux DOUBLE excitations -- the dynamic-correlation channel the icPT2
    #      downfolds.  A bare block-diagonal ZDO gives only Coulomb (pq|ab) and no
    #      transfer (pa|qr), so the external Q space would NOT correlate the CAS
    #      (verified: CAS==FCI in that case).  We add the transfer explicitly via
    #      one small density-fitting factor per aux orbital -- each aux "borrows"
    #      a little frontier-pi charge in the FITTING FACTORS ONLY (never in the
    #      1e, so aux stay external).  Built as eri += sum_L b_L b_L^T so the
    #      tensor stays 8-fold symmetric and positive-semidefinite.  The aux
    #      self-density is kept modest (and the aux 1e energies high) so electrons
    #      are NOT pulled into the aux orbitals: the FCI stays CAS-dominated and
    #      the correlation is a small, smooth perturbation.
    eri_mo = np.zeros((m, m, m, m))
    if U > 0.0:
        Mpi = np.einsum('ip,iq->pqi', Cpi, Cpi)             # (npi,npi,npi)
        eri_mo[:npi, :npi, :npi, :npi] = np.einsum('pqi,ij,rsj->pqrs',
                                                   Mpi, gam, Mpi)
        scale = kappa                                       # weak coupling scale
        # frontier weighting: couple mainly the two middle (frontier) pi MOs,
        # so the channel is the physical HOMO/LUMO -> aux double excitation.
        wpi = np.array([0.3, 1.0, 1.0, 0.3])[:npi]
        for a in range(n_ext):
            ao = npi + a
            bvec = np.zeros((m, m))
            for p in range(npi):
                bvec[ao, p] = bvec[p, ao] = scale * wpi[p] / np.sqrt(n_ext)
            # NOTE: deliberately NO aux self-density term (bvec[ao,ao]=0).
            # An aux self-Coulomb (aa|aa) would lower the doubly-aux-occupied
            # configuration toward the CAS and create an INTRUDER (the FCI
            # collapses into aux occupation, CAS weight -> 0).  With only the
            # transfer (pa|qa) the aux orbitals contribute purely as virtual,
            # perturbative double-excitation channels -> smooth dynamic
            # correlation with the FCI staying CAS-dominated.
            eri_mo += np.einsum('pq,rs->pqrs', bvec, bvec)
    eps_all = np.concatenate([eps_pi, eps_aux])
    return h_mo, eri_mo, eps_all


# ======================================================================
# 2. Build full FCI + CAS(4,4) partition for one geometry
# ======================================================================
def build_case(delta, t=2.40, U=8.0, n_ext=2, ext_gap=6.0, kappa=1.0):
    """Everything the FCI / downfolds need for one bond-alternation `delta`.

    CAS(4,4): the FOUR pi orbitals are active, 4 pi electrons.  There is NO pi
    core (the full pi space is the CAS, as it must be for CBD).  The external
    virtuals (indices 4..) form the external-Q space; the FCI here is the FCI in
    (pi + external) for the 4 electrons -- the same 4 electrons can excite into
    the externals (dynamic correlation), which is what the icPT2 downfolds.
    """
    h_mo, eri_mo, eps_all = build_cbd(t=t, U=U, delta=delta, n_ext=n_ext,
                                      ext_gap=ext_gap, kappa=kappa)
    eps_pi = eps_all[:4]                          # the four pi MO energies
    H1, g, m = spinorb(h_mo, eri_mo)
    nelec = 4
    na = nb = nelec // 2
    dets = gen_dets(m, na, nb)
    Hfull = build_H(dets, H1, g)

    core = []                                   # no pi core: full pi space active
    active = [0, 1, 2, 3]                        # the four pi MOs
    virt = list(range(4, m))                     # external sigma*-like virtuals
    Pset = set(gen_dets(m, na, nb, core, active, virt, restrict=True))
    Pidx = [i for i, d in enumerate(dets) if d in Pset]
    Qidx = [i for i in range(len(dets)) if i not in set(Pidx)]

    # Fock (MP) spin-orbital energies for a Dyall-type denominator option.
    # Reference = the closed-shell aufbau on the lowest two pi MOs.
    nso = H1.shape[0]
    occ_spatial = [0, 1]                          # lowest-2 pi (closed shell ref)
    ref_occ = sorted(occ_spatial + [x + m for x in occ_spatial])
    fock_eps = np.empty(nso)
    for p in range(nso):
        e = H1[p, p]
        for i in ref_occ:
            e += g[p, i, p, i]
        fock_eps[p] = e

    return dict(h_mo=h_mo, eri_mo=eri_mo, eps_pi=eps_pi, H1=H1, g=g, m=m,
                dets=dets, Hfull=Hfull, Pidx=Pidx, Qidx=Qidx,
                part=(core, active, virt), fock_eps=fock_eps,
                ref_occ=set(ref_occ))


# ======================================================================
# 3. One-operator multistate icPT2 (des-Cloizeaux symmetric Hermitian)
# ======================================================================
def icpt2_multistate(case, nP=3, level_shift=0.3):
    """ONE consistent operator: diagonalize the CAS(4,4) backbone H_PP (one CI
    operator over the determinant union), then dress ALL its low roots with ONE
    uniform Hermitian external-Q self-energy (Epstein-Nesbet denominators,
    des-Cloizeaux symmetric average):

        H_eff[k,l] = E_k0 delta_kl
                   + (1/2) sum_q coup_qk coup_ql (1/(E_k0-H_qq)+1/(E_l0-H_qq))

    coup_qk = <q|H|Psi_P^k>.  A small real `level_shift` regularizes intruders
    (standard CASPT2 IPEA-free shift).  Returns (Edressed, eP) ascending.
    """
    Hfull, Pidx, Qidx = case['Hfull'], case['Pidx'], case['Qidx']
    HPP = Hfull[np.ix_(Pidx, Pidx)]
    eP, cP = np.linalg.eigh(HPP)
    nP = min(nP, len(eP))
    HQP = Hfull[np.ix_(Qidx, Pidx)]
    Hqq = np.diag(Hfull)[Qidx]
    coup = HQP @ cP[:, :nP]                                  # (nQ, nP)
    inv = np.empty((len(Qidx), nP))
    for k in range(nP):
        d = eP[k] - Hqq
        # regularized denominator: d/(d^2+ls^2) (real, finite through d=0)
        inv[:, k] = d / (d * d + level_shift * level_shift)
    Heff = np.diag(eP[:nP]).astype(float)
    for k in range(nP):
        for l in range(nP):
            Heff[k, l] += 0.5 * np.sum(coup[:, k] * coup[:, l]
                                       * (inv[:, k] + inv[:, l]))
    herm = np.abs(Heff - Heff.T).max()
    Heff = 0.5 * (Heff + Heff.T)
    Edr = np.linalg.eigvalsh(Heff)
    return Edr, eP[:nP], herm


# ======================================================================
# 4. The naive GLUED hybrid (mimics the original QMRSF defect)
# ======================================================================
def classify_P_dets(case):
    """Split the P-space (CAS) determinants into the QMRSF response classes by
    their pi-occupation pattern:

      0OS : closed-shell (every occupied pi spatial orbital doubly occupied).
            These are the |2200>, |0220>, ... ionic/covalent-closed configs --
            the DOUBLE spin-flip references.  In the original QMRSF these get a
            BARE (H-E0) diagonal, separately referenced.
      1SF : single-spin-flip (2OS/4OS) -- one or more singly-occupied pi
            orbitals.  In the original QMRSF these get the DFT-dressed orbital-
            Hessian diagonal.

    Returns (idx0OS, idx1SF) as lists of indices INTO case['Pidx'].
    """
    m = case['m']
    Pdets = [case['dets'][i] for i in case['Pidx']]
    idx0, idx1 = [], []
    for j, d in enumerate(Pdets):
        alpha = set(x for x in d if x < m)
        beta = set(x - m for x in d if x >= m)
        # doubly-occupied spatial = alpha & beta ; singly = symmetric difference
        singly = alpha ^ beta
        if len(singly) == 0:
            idx0.append(j)                       # closed shell -> 0OS
        else:
            idx1.append(j)                       # has open shells -> 1SF (2OS/4OS)
    return idx0, idx1


def block_diabats(case, d_dft, c0):
    """The two GLUED diabatic energies, computed in SEPARATE block operators with
    the inter-block coupling DROPPED -- this is the structural defect of the
    naive QMRSF scheme (it builds the 1SF response and the 0OS configs from two
    different operators and patches them per-system, instead of one consistent
    operator that carries the 1SF<->0OS coupling).

      * 1SF diabat E_1SF: lowest eigenvalue of the 1SF (2OS/4OS) sub-block of
        H_PP, each 1SF diagonal DFT-dressed by `d_dft` (the orbital-Hessian
        dynamic-correlation dressing).
      * 0OS diabat E_0OS: lowest eigenvalue of the 0OS (closed-shell) sub-block
        of H_PP, each 0OS diagonal shifted by the single per-system constant
        `c0` (the bare-referenced closed-shell block).

    The 1SF<->0OS off-diagonal block of H_PP is NEVER used.  In the ONE consistent
    operator (the determinant-union backbone) that coupling turns the diabatic
    crossing into a smooth avoided crossing; dropping it leaves two diabats that
    physically CROSS as the ground-state character changes.

    Returns (E_1SF, E_0OS).
    """
    Hfull, Pidx = case['Hfull'], case['Pidx']
    HPP = Hfull[np.ix_(Pidx, Pidx)]
    idx0, idx1 = classify_P_dets(case)
    H1 = HPP[np.ix_(idx1, idx1)].copy()
    for a in range(len(idx1)):
        H1[a, a] += d_dft
    H0 = HPP[np.ix_(idx0, idx0)].copy()
    for a in range(len(idx0)):
        H0[a, a] += c0
    e1 = np.linalg.eigvalsh(H1)[0] if len(idx1) else np.inf
    e0 = np.linalg.eigvalsh(H0)[0] if len(idx0) else np.inf
    return e1, e0


def naive_glued(case, d_dft, c0, nroots=3):
    """Naive GLUED S0/S1: the two block diabats (1SF DFT-dressed, 0OS bare +
    per-system constant), SELECTED by energy (min for S0).  Because the
    inter-block coupling is dropped, the diabats cross as the ground-state
    character changes -> the energy-selected S0 has a derivative KINK at the
    crossing (and the automerization barrier comes out wrong).  Returns the two
    energy-sorted block-diabat energies (S0, S1); higher roots padded with S1."""
    e1, e0 = block_diabats(case, d_dft, c0)
    es = sorted([e1, e0])
    out = (es + [es[-1]] * nroots)[:nroots]
    return np.array(out)


# ======================================================================
# 5. Ground-state character diagnostic (covalent/4OS <-> ionic/0OS)
# ======================================================================
def gs_character(case):
    """Closed-shell (0OS) weight of the FCI ground state, measured WITHIN the
    CAS(4,4) projection (renormalized over the P-space determinants).  ~0 =>
    covalent open-shell (square/4OS-dominated); ~1 => ionic closed-shell (0OS).

    We project onto the CAS and renormalize so the diagnostic is insensitive to
    the small dynamic-correlation leakage into the external orbitals (which would
    otherwise dilute a full-vector weight).  Returns (w0OS, cas_weight)."""
    w, V = np.linalg.eigh(case['Hfull'])
    v = V[:, 0]
    m = case['m']
    Pidx = case['Pidx']
    vp = v[Pidx]
    cas_weight = float(np.dot(vp, vp))
    vpn = vp / (np.linalg.norm(vp) + 1e-30)
    Pdets = [case['dets'][i] for i in Pidx]
    w0 = 0.0
    for j, d in enumerate(Pdets):
        alpha = set(x for x in d if x < m)
        beta = set(x - m for x in d if x >= m)
        if len(alpha ^ beta) == 0:               # closed shell (0OS)
            w0 += vpn[j] ** 2
    return float(w0), cas_weight


# ======================================================================
# 6. Driver
# ======================================================================
def banner(s):
    print("=" * 78); print(s); print("=" * 78)


def main():
    np.set_printoptions(precision=5, suppress=True)
    t, U, n_ext, ext_gap, kappa = 2.40, 6.0, 2, 4.0, 1.0

    banner("QMRSF CBD model | 4-site cyclic PPP-Ohno, CAS(4,4)=full pi, +ext virtuals")

    # ---------- correctness gates ----------
    c0 = build_case(0.0, t=t, U=U, n_ext=n_ext, ext_gap=ext_gap, kappa=kappa)
    print(f"\n[gate 1] Hermiticity max|Hfull-Hfull^T| = "
          f"{np.abs(c0['Hfull'] - c0['Hfull'].T).max():.2e}")
    print(f"         m={c0['m']} orbitals | active(pi)={c0['part'][1]} "
          f"ext={c0['part'][2]} | dets full={len(c0['dets'])} "
          f"P(CAS)={len(c0['Pidx'])} Q(ext)={len(c0['Qidx'])}")
    # U=0 gate: FCI == sum of lowest-4 spin-orbital (Hueckel/ext) energies
    h0, e0, _ = build_cbd(t=t, U=0.0, delta=0.0, n_ext=n_ext, ext_gap=ext_gap)
    H1z, gz, mz = spinorb(h0, np.zeros((c0['m'],) * 4))
    detsz = gen_dets(mz, 2, 2)
    Hz = build_H(detsz, H1z, gz)
    e_ci0 = np.linalg.eigvalsh(Hz)[0]
    so = np.sort(np.concatenate([np.diag(h0), np.diag(h0)]))
    print(f"[gate 2] U=0  FCI={e_ci0:.6f}  sum-eps={so[:4].sum():.6f}  "
          f"diff={abs(e_ci0 - so[:4].sum()):.2e}")
    # square-point Hueckel levels: should be (-2t, 0, 0, +2t) -> degenerate pair
    print(f"[gate 3] square (delta=0) pi Hueckel levels = {np.round(c0['eps_pi'],4)} "
          f"(expect -2t,0,0,+2t = {-2*t:.2f},0,0,{2*t:.2f}; degenerate frontier pair)")

    # ---------- ground-state CHARACTER across the square point ----------
    banner("CHARACTER CROSSOVER | FCI ground-state closed-shell (0OS) weight vs delta")
    print(f"\n  {'delta':>7} {'E_FCI(S0)':>12} {'gap S1-S0':>10} "
          f"{'0OS weight':>11}  character")
    for dl in [0.00, 0.02, 0.05, 0.10, 0.20, 0.30, 0.45]:
        c = build_case(dl, t=t, U=U, n_ext=n_ext, ext_gap=ext_gap, kappa=kappa)
        ev = np.linalg.eigvalsh(c['Hfull'])
        w0, _ = gs_character(c)
        ch = "covalent/4OS (open-shell)" if w0 < 0.5 else "ionic/0OS (closed-shell)"
        print(f"  {dl:>7.2f} {ev[0]:>12.5f} {ev[1]-ev[0]:>10.4f} {w0:>11.3f}  {ch}")
    print("\n  => the ground state CHANGES CHARACTER across the square point: covalent")
    print("     open-shell (4OS) at/near D4h -> ionic closed-shell (0OS) for D2h.")
    print("     This is the 4OS<->0OS crossover that a glued/patched QMRSF cannot")
    print("     follow with a single per-system 0OS reference.")

    # ================================================================
    # SCAN through the square point: FCI vs one-operator vs naive glued
    # ================================================================
    banner("SCAN | square<->rectangular: FCI vs ONE-OPERATOR icPT2 vs NAIVE GLUED")

    deltas = np.linspace(0.0, 0.45, 46)          # square (0) -> rectangular
    ls = 0.3

    # --- naive glued calibration, the way a real code calibrates each block at
    #     its OWN home geometry (each block is individually correct there; the
    #     GLUE -- the dropped 1SF<->0OS coupling -- is what fails in between):
    #   (1) DFT dressing d_dft of the 1SF block fit at the SQUARE point, where
    #       the open-shell (4OS/1SF) covalent state is the ground state, so the
    #       1SF diabat matches FCI there.
    #   (2) the SEPARATE bare-0OS constant c0 fit at the RECTANGULAR end, where
    #       the closed-shell (0OS) state is the ground state, so the 0OS diabat
    #       matches FCI there.
    case_sq = build_case(deltas[0], t=t, U=U, n_ext=n_ext, ext_gap=ext_gap, kappa=kappa)
    e_fci_sq = np.linalg.eigvalsh(case_sq['Hfull'])[0]
    e1_bare, _ = block_diabats(case_sq, 0.0, 0.0)
    d_dft = e_fci_sq - e1_bare                      # 1SF diabat -> FCI at square

    case_rc = build_case(deltas[-1], t=t, U=U, n_ext=n_ext, ext_gap=ext_gap, kappa=kappa)
    e_fci_rc = np.linalg.eigvalsh(case_rc['Hfull'])[0]
    _, e0_bare = block_diabats(case_rc, 0.0, 0.0)
    c0fit = e_fci_rc - e0_bare                      # 0OS diabat -> FCI at rect
    d1 = d_dft
    print(f"\n  naive glued calibration (each block at its OWN home geometry):")
    print(f"     1SF DFT dressing d_dft = {d1:+.5f} (fit at SQUARE delta={deltas[0]:.3f},"
          f" 4OS=GS)")
    print(f"     separate bare-0OS constant c0 = {c0fit:+.5f} (fit at RECT delta="
          f"{deltas[-1]:.3f}, 0OS=GS)")
    print(f"  one-operator: multistate des-Cloizeaux icPT2, level_shift={ls}\n")

    E_fci = np.zeros((len(deltas), 3))
    E_op = np.zeros((len(deltas), 3))            # one-operator icPT2
    E_gl = np.zeros((len(deltas), 3))            # naive glued hybrid
    W0 = np.zeros(len(deltas))
    for i, dl in enumerate(deltas):
        c = build_case(dl, t=t, U=U, n_ext=n_ext, ext_gap=ext_gap, kappa=kappa)
        E_fci[i] = np.linalg.eigvalsh(c['Hfull'])[:3]
        Edr, _, _ = icpt2_multistate(c, nP=3, level_shift=ls)
        E_op[i] = Edr[:3]
        E_gl[i] = naive_glued(c, d1, c0fit, nroots=3)
        W0[i], _ = gs_character(c)

    # ---- automerization barrier analogue: E(square) - E(rect. minimum) ----
    # The CBD automerization barrier ~ E(D4h square TS) - E(D2h rectangular min).
    # Here delta=0 is the square; the rectangular "minimum" of each S0 curve is
    # its lowest point over delta>0.  Barrier = E(square) - min_{delta>0} E(S0).
    def barrier(Es0):
        i_sq = 0                                  # delta = 0 is the square point
        i_rmin = 1 + int(np.argmin(Es0[1:]))      # rectangular minimum (delta>0)
        return Es0[i_sq] - Es0[i_rmin], deltas[i_rmin]

    b_fci, dr_fci = barrier(E_fci[:, 0])
    b_op, dr_op = barrier(E_op[:, 0])
    b_gl, dr_gl = barrier(E_gl[:, 0])
    print(f"  AUTOMERIZATION BARRIER analogue  E(square) - E(rectangular min):")
    print(f"     FCI            : {b_fci:+.5f}   (rect. min at delta={dr_fci:.3f})")
    print(f"     one-operator   : {b_op:+.5f}   (rect. min at delta={dr_op:.3f})  "
          f"err {b_op-b_fci:+.5f}")
    print(f"     naive glued    : {b_gl:+.5f}   (rect. min at delta={dr_gl:.3f})  "
          f"err {b_gl-b_fci:+.5f}")
    sign_ok_op = np.sign(b_op) == np.sign(b_fci)
    sign_ok_gl = np.sign(b_gl) == np.sign(b_fci)
    print(f"     barrier sign vs FCI: one-operator {'OK' if sign_ok_op else 'WRONG-SIGN'}"
          f" | naive glued {'OK' if sign_ok_gl else 'WRONG-SIGN'}")

    # ---- smoothness ----
    # IMPORTANT NUANCE: the EXACT FCI S0 surface is itself NON-smooth at the
    # antiaromatic D4h square point -- the covalent<->ionic transition there is a
    # genuine physical near-cusp, so a RAW max|2nd diff| of S0 is dominated by
    # real physics and cannot, on its own, separate a method's spurious kink from
    # the true curvature.  The correct discriminator is the smoothness of each
    # method's ERROR relative to FCI: a consistent (one-operator) method has a
    # small, smooth error through the transition; a glued/patched method's error
    # KINKS at the diabatic crossing where its two operators meet.
    def maxd2(arr):
        return np.abs(np.diff(arr, 2)).max()
    def slope_jump(arr):                                # max |d/ddelta jump| (kink)
        sl = np.diff(arr) / np.diff(deltas)
        return np.abs(np.diff(sl)).max()
    d2_fci = maxd2(E_fci[:, 0]); d2_op = maxd2(E_op[:, 0]); d2_gl = maxd2(E_gl[:, 0])
    err_op = E_op[:, 0] - E_fci[:, 0]
    err_gl = E_gl[:, 0] - E_fci[:, 0]
    sj_op = slope_jump(err_op); sj_gl = slope_jump(err_gl)
    print(f"\n  SMOOTHNESS:")
    print(f"   raw max|2nd diff| of S0 (dominated by the PHYSICAL square-point cusp):")
    print(f"     FCI={d2_fci:.3e}  one-operator={d2_op:.3e}  naive glued={d2_gl:.3e}")
    print(f"   *** kink-of-ERROR metric (physical FCI curvature subtracted) ***")
    print(f"     one-operator: max|err|={np.abs(err_op).max():.5f}  "
          f"max|slope-jump of err|={sj_op:.4f}")
    print(f"     naive glued : max|err|={np.abs(err_gl).max():.5f}  "
          f"max|slope-jump of err|={sj_gl:.4f}  (ratio to one-op {sj_gl/sj_op:.1f}x)")
    i_spike = int(np.argmax(np.abs(np.diff(err_gl, 2))))
    print(f"     naive-glued error kink at delta={deltas[i_spike+1]:.3f} "
          f"(0OS weight there = {W0[i_spike+1]:.3f}) -- i.e. where the 1SF and 0OS")
    print(f"     diabats cross (the glue seam).")

    # ================================================================
    # VERDICT
    # ================================================================
    banner("VERDICT")
    op_smooth = sj_op < 2.0 * d2_fci / np.diff(deltas)[0]   # error stays smooth
    # glued is discontinuous if its error kinks materially more than one-op's,
    # and the error is materially larger.
    gl_kink = (sj_gl > 2.0 * sj_op) and (np.abs(err_gl).max() > 3.0 * np.abs(err_op).max())
    gl_wrongsign = not sign_ok_gl
    print(f"  one-operator icPT2 tracks FCI smoothly      : "
          f"{'YES' if op_smooth else 'NO'}  (max|err|={np.abs(err_op).max():.4f}, "
          f"slope-jump of err {sj_op:.3f})")
    print(f"  one-operator barrier sign/magnitude OK      : "
          f"{'YES' if sign_ok_op else 'NO'}  ({b_op:+.4f} vs FCI {b_fci:+.4f})")
    print(f"  naive glued shows a KINK (discontinuous err): "
          f"{'YES' if gl_kink else 'NO'}  (max|err|={np.abs(err_gl).max():.4f}, "
          f"slope-jump of err {sj_gl:.3f} = {sj_gl/sj_op:.1f}x one-op)")
    print(f"  naive glued gives WRONG-SIGN barrier        : "
          f"{'YES' if gl_wrongsign else 'NO'}  ({b_gl:+.4f} vs FCI {b_fci:+.4f})")
    if gl_kink or gl_wrongsign:
        print("\n  ==> VALIDATES the manuscript's CBD claim: the naive glued QMRSF")
        print("      (DFT-dressed 1SF diabat + independently-referenced, per-system-")
        print("      pinned bare-0OS diabat, inter-block coupling DROPPED) produces a")
        print("      KINKED automerization surface precisely where the ground state")
        print("      crosses 4OS<->0OS.  It is forced exact at the two calibration")
        print("      geometries, yet its ERROR develops a derivative discontinuity at")
        print("      the diabatic crossing in between.  The ONE consistent operator")
        print("      (determinant-union backbone, which KEEPS the 1SF<->0OS coupling,")
        print("      + uniform Hermitian dynamic correlation) stays smooth and tracks")
        print("      FCI through the same transition.")
    else:
        print("\n  ==> HONEST NOTE: the naive glued hybrid did not produce a hard kink")
        print("      here; report shows its systematic error instead (see barrier/err).")

    # ---------- optional plot ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
        # panel A: S0 surfaces
        ax = axes[0]
        ax.plot(deltas, E_fci[:, 0], 'k-', lw=2.2, label='FCI S0')
        ax.plot(deltas, E_op[:, 0], 'C0o-', ms=3, label='one-operator icPT2 S0')
        ax.plot(deltas, E_gl[:, 0], 'C3s-', ms=3, label='naive glued S0')
        ax.axvline(0.0, color='gray', ls=':', lw=1)
        ax.set_xlabel('bond alternation delta (0=square D4h)')
        ax.set_ylabel('energy'); ax.set_title('S0: automerization surface')
        ax.legend(fontsize=8)
        # panel B: lowest 3 states, FCI vs one-operator
        ax = axes[1]
        for k, (c, lab) in enumerate(zip('kkk', ['S0', 'S1', 'S2'])):
            ax.plot(deltas, E_fci[:, k], color='k',
                    ls=['-', '--', ':'][k], lw=2, label=f'FCI {lab}')
            ax.plot(deltas, E_op[:, k], color=f'C{k}', ls='-', marker='.',
                    ms=2, label=f'one-op {lab}')
        ax.set_xlabel('delta'); ax.set_ylabel('energy')
        ax.set_title('lowest 3 states: one-operator tracks FCI'); ax.legend(fontsize=7)
        # panel C: error-vs-FCI (the kink discriminator) + character crossover
        ax = axes[2]
        ax.plot(deltas, (E_op[:, 0] - E_fci[:, 0]) * 1e3, 'C0o-', ms=3,
                label='one-operator S0 err (mEh-like)')
        ax.plot(deltas, (E_gl[:, 0] - E_fci[:, 0]) * 1e3, 'C3s-', ms=3,
                label='naive glued S0 err (KINK)')
        ax.axhline(0.0, color='gray', lw=0.8)
        ax.set_xlabel('delta'); ax.set_ylabel('S0 error vs FCI (x1e-3)')
        ax.set_title('error vs FCI: glued kinks at the crossing')
        axb = ax.twinx()
        axb.plot(deltas, W0, 'C2:', lw=1.5, label='FCI 0OS weight (character)')
        axb.set_ylabel('0OS weight', color='C2')
        ax.legend(loc='upper right', fontsize=7); axb.legend(loc='center right', fontsize=7)
        fig.tight_layout()
        png = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "qmrsf_cbd_model.png")
        fig.savefig(png, dpi=130)
        print(f"\n  [plot] saved {png}")
    except Exception as exc:
        print(f"\n  [plot] skipped ({type(exc).__name__}: {exc})")

    banner("DONE")


if __name__ == "__main__":
    main()
