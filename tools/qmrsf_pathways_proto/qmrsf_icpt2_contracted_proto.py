#!/usr/bin/env python3
"""
QMRSF-icPT2  CONTRACTED  external-Q perturber generation (pure NumPy, NO pyscf).

================================================================================
WHAT THIS REPLACES
================================================================================
The validated brute-force engine (qmrsf_icpt2_ppp_proto.py + the live gate
stageB/gate_icpt2_full.py) computes the QMRSF-icPT2 dressed spectra by:

  1. enumerating the FULL determinant space  (na alpha + nb beta electrons in
     `norb` window orbitals  ->  C(norb,na) * C(norb,nb) determinants),
  2. partitioning P (all electrons confined to the first NACT active orbitals =
     the CAS(4,4) backbone)  from  external Q (everything else),
  3. building the dense H_QP block and diag(H)_Q over ALL Q determinants,
  4. des-Cloizeaux symmetric multistate downfold onto the CAS roots.

That is O(ndet^2) and, in the number of virtual orbitals nvirt = norb - NACT,
the external-Q count grows as ~ nvirt^4 (4 active electrons can be promoted into
up to 4 distinct virtual spin-orbitals).  Correct, but it materialises the full
FCI determinant list.

================================================================================
WHAT THIS DOES INSTEAD  (the contraction)
================================================================================
The frozen-core window has NO core orbitals inside it; the reference has
exactly NACT-derived active electrons (here 4: na=nb=2) all in active orbitals,
and virtuals = orbitals NACT .. norb-1.  Therefore EVERY external-Q determinant
is reached from the active space purely by promoting 1..(#active electrons)
electrons into virtual spin-orbitals.  We can enumerate Q *intrinsically* as

      |q>  =  (occupied VIRTUAL spin-orbital subset  V, size 1..nelec)
              x
              (residual occupied ACTIVE spin-orbital subset  A)

with the per-spin electron counts (na,nb) conserved.  The virtual subset V is
the only part that grows with the basis; the active residual A is enumerated
inside the fixed NACT-orbital active space (a small, basis-independent loop).

For each such q we form the INTERNALLY-CONTRACTED coupling against each CAS
eigenvector |Psi_P^k> = sum_{dP} c_dP^k |dP> :

      coup_qk = <q|H|Psi_P^k> = sum_{dP in CAS} c_dP^k <q|H|dP>

using the SAME Slater-Condon melem() as the brute force.  Only CAS determinants
dP that differ from q by <= 2 spin-orbitals contribute (Slater-Condon), so the
inner sum is cheap and, crucially, we never build or store the FCI list.

The two zeroth-order Q diagonals needed by the downfold are computed in closed
form from q's own occupation (no extra determinants):

  * EN    :  <q|H|q>  via the diagonal branch of melem (Slater-Condon nd==0),
  * Dyall :  frozen-core Dyall denominator  d_q = -(sum of eps over q's occupied
             VIRTUAL spin-orbitals), eps indexed by spatial = so % norb.

Then we apply the IDENTICAL des-Cloizeaux symmetric multistate downfold as the
brute force (stageB/gate_icpt2_full.py:downfold).

================================================================================
VALIDATION GATE
================================================================================
Reads stageB/qmrsf_icpt2_full_live.dat, builds the CAS(4,4) (the 36 P
determinants, H_PP diagonalised -> eP, cP), runs the contracted perturber
generation + downfold, and checks the EN and Dyall dressed spectra against the
dumped edr_en / edr_dy oracle (and, independently, against the brute-force
re-run) to < 1e-9 (sorted).  Prints PASS/FAIL and max abs differences.

Run:  python3 qmrsf_icpt2_contracted_proto.py     (from tools/qmrsf_pathways_proto/)
"""
import os
import sys
import numpy as np
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Reuse the VALIDATED Slater-Condon core verbatim (never duplicate the signs).
import qmrsf_icpt2_ppp_proto as P
from qmrsf_icpt2_ppp_proto import spinorb, gen_dets, melem, build_H

DATFILE = os.path.join(HERE, "stageB", "qmrsf_icpt2_full_live.dat")
NACT = 4                       # active spatial orbitals  (CAS(4,4))


# ----------------------------------------------------------------------
# 0. Read the live window fixture (same parser as stageB/gate_icpt2_full.py)
# ----------------------------------------------------------------------
def read_fixture(path=DATFILE):
    with open(path) as f:
        norb, nPd = map(int, f.readline().split())
        h = np.array([[float(x) for x in f.readline().split()] for _ in range(norb)])
        eri = np.zeros((norb,) * 4)
        for p in range(norb):
            for q in range(norb):
                for r in range(norb):
                    eri[p, q, r, :] = [float(x) for x in f.readline().split()]
        ecore = float(f.readline())
        eps = np.array([float(x) for x in f.readline().split()])
        eP_live = np.array([float(x) for x in f.readline().split()])
        en_live = np.array([float(x) for x in f.readline().split()])
        dy_live = np.array([float(x) for x in f.readline().split()])
    return dict(norb=norb, nPd=nPd, h=h, eri=eri, ecore=ecore, eps=eps,
                eP_live=eP_live, en_live=en_live, dy_live=dy_live)


# ----------------------------------------------------------------------
# 1. CAS(4,4) backbone:  the 36 P determinants and their H_PP eigenpairs
# ----------------------------------------------------------------------
def build_cas(norb, na, nb, H1, g):
    """P determinants = all na alpha + nb beta electrons confined to the first
    NACT active spatial orbitals.  Returns (Pdets, eP, cP).

    This is a SMALL, basis-INDEPENDENT enumeration: C(NACT,na)*C(NACT,nb) dets
    regardless of norb.  It is the only determinant list the contracted method
    ever materialises.
    """
    Pdets = []
    for a in combinations(range(NACT), na):
        for b in combinations(range(NACT), nb):
            Pdets.append(tuple(sorted(list(a) + [x + norb for x in b])))
    Pdets.sort()
    HPP = build_H(Pdets, H1, g)
    eP, cP = np.linalg.eigh(HPP)
    return Pdets, eP, cP


# ----------------------------------------------------------------------
# 2. CONTRACTED external-Q enumeration  (the scalable core)
# ----------------------------------------------------------------------
def gen_external_Q(norb, na, nb):
    """Enumerate external-Q determinants INTRINSICALLY by virtual-occupation
    pattern crossed with the residual active occupation.

    A spin-orbital so has spatial = so % norb, spin = so // norb (alpha < norb,
    beta >= norb).  Active spatial = 0..NACT-1, virtual spatial = NACT..norb-1.

    For each spin sigma in {alpha, beta} independently we choose:
        n_v^sigma  electrons in virtual spin-orbitals  (0..n_sigma),
        n_a^sigma = n_sigma - n_v^sigma  electrons in active spin-orbitals,
    and take all subsets of the appropriate size.  A determinant is in Q iff it
    has AT LEAST ONE virtual-occupied spin-orbital (n_v^alpha + n_v^beta >= 1);
    the all-active case (n_v=0) is exactly the CAS / P space and is excluded.

    The virtual loop ranges over subsets of the 2*nvirt virtual spin-orbitals
    -> this is the only part that grows with the basis.  The active residual
    loop is inside the fixed NACT active space.

    Yields tuples det (sorted spin-orbital indices).  No FCI list is built.
    """
    act_a = list(range(NACT))                       # alpha active spin-orbitals
    vir_a = list(range(NACT, norb))                 # alpha virtual spin-orbitals
    act_b = [x + norb for x in range(NACT)]         # beta  active spin-orbitals
    vir_b = [x + norb for x in range(NACT, norb)]   # beta  virtual spin-orbitals

    # precompute the per-spin (active-subset, virtual-subset) building blocks
    def spin_blocks(nelec_sigma, act_so, vir_so):
        """all ways to place nelec_sigma electrons, labelled by #virtual used."""
        blocks = []   # list of (nv, occ_tuple)
        for nv in range(0, nelec_sigma + 1):
            na_act = nelec_sigma - nv
            if na_act > len(act_so) or nv > len(vir_so):
                continue
            for vsub in combinations(vir_so, nv):
                for asub in combinations(act_so, na_act):
                    blocks.append((nv, tuple(asub) + tuple(vsub)))
        return blocks

    a_blocks = spin_blocks(na, act_a, vir_a)
    b_blocks = spin_blocks(nb, act_b, vir_b)

    for nva, aocc in a_blocks:
        for nvb, bocc in b_blocks:
            if nva + nvb == 0:
                continue                            # all-active = P, skip
            yield tuple(sorted(aocc + bocc))


# ----------------------------------------------------------------------
# 3. Contracted couplings + closed-form Q diagonals
# ----------------------------------------------------------------------
def contracted_perturbers(norb, na, nb, Pdets, cP, eps, H1, g):
    """For every external-Q determinant q produced by gen_external_Q, compute

        coup[q, k] = <q|H|Psi_P^k> = sum_{dP} c_dP^k <q|H|dP>   (contracted)
        Hqq[q]     = <q|H|q>                                    (EN diagonal)
        sumv[q]    = sum of eps over q's occupied VIRTUAL spin-orbitals
                       -> Dyall denominator d_q = -sumv[q]

    Only CAS dP within 2 spin-orbital substitutions of q contribute (Slater-
    Condon), so the inner loop is sparse.  Returns (coup, Hqq, sumv, nQ).
    """
    nPd = cP.shape[1]
    Pset_index = {d: j for j, d in enumerate(Pdets)}   # for completeness/debug
    # Pre-store CAS det occupation sets to prune the inner Slater-Condon loop.
    Pocc = [set(d) for d in Pdets]

    coup_rows = []
    Hqq_list = []
    sumv_list = []
    nQ = 0
    for q in gen_external_Q(norb, na, nb):
        nQ += 1
        sq = set(q)
        # ---- contracted coupling against all CAS dets (prune by excitation rank) ----
        # <q|H|dP> nonzero only if |q \ dP| <= 2  (Slater-Condon).
        cvec = np.zeros(nPd)
        for j, dP in enumerate(Pdets):
            ndiff = len(sq - Pocc[j])               # # spin-orbitals in q not in dP
            if ndiff > 2:
                continue
            mel = melem(q, dP, H1, g)               # <q|H|dP>, validated signs
            if mel != 0.0:
                cvec += mel * cP[j, :]              # accumulate c_dP^k * <q|H|dP>
        coup_rows.append(cvec)
        # ---- EN diagonal  <q|H|q>  (Slater-Condon diagonal branch) ----
        Hqq_list.append(melem(q, q, H1, g))
        # ---- Dyall: sum of eps over occupied VIRTUAL spin-orbitals of q ----
        sumv_list.append(sum(eps[so % norb] for so in q if (so % norb) >= NACT))

    coup = np.array(coup_rows)                       # (nQ, nPd)
    Hqq = np.array(Hqq_list)                         # (nQ,)
    sumv = np.array(sumv_list)                       # (nQ,)
    return coup, Hqq, sumv, nQ


# ----------------------------------------------------------------------
# 4. des-Cloizeaux symmetric multistate downfold  (IDENTICAL to brute force)
# ----------------------------------------------------------------------
def guard(d):
    return np.where(np.abs(d) < 1e-6, np.sign(d) * 1e-6 + 1e-30, d)


def downfold(eP, coup, invd):
    """H_eff[k,l] = delta_kl eP_k + 1/2 sum_q coup_qk coup_ql (invd_qk + invd_ql),
    diagonalised.  invd is the per-(q,k) inverse zeroth-order denominator matrix.
    This is verbatim the brute-force stageB/gate_icpt2_full.py downfold.
    """
    nPd = len(eP)
    Heff = np.diag(eP).astype(float)
    for k in range(nPd):
        for l in range(nPd):
            Heff[k, l] += 0.5 * np.sum(coup[:, k] * coup[:, l] * (invd[:, k] + invd[:, l]))
    return np.linalg.eigvalsh(0.5 * (Heff + Heff.T))


def run_contracted(fix):
    """Full contracted pipeline -> (eP, edr_en, edr_dy, nQ, coup, Hqq, sumv)."""
    norb = fix['norb']; nPd = fix['nPd']
    H1, g, _ = spinorb(fix['h'], fix['eri'])
    na = nb = 2
    Pdets, eP, cP = build_cas(norb, na, nb, H1, g)
    assert len(eP) == nPd, f"CAS dim {len(eP)} != nPd {nPd}"

    coup, Hqq, sumv, nQ = contracted_perturbers(norb, na, nb, Pdets, cP,
                                                fix['eps'], H1, g)

    # EN:    invd_qk = 1/(eP_k - Hqq)             (state-specific)
    inv_en = np.column_stack([1.0 / guard(eP[k] - Hqq) for k in range(nPd)])
    # Dyall: invd_qk = 1/(-sumv)                  (state-independent denominator)
    inv_dy = np.column_stack([1.0 / guard(-sumv) for _ in range(nPd)])

    edr_en = downfold(eP, coup, inv_en)
    edr_dy = downfold(eP, coup, inv_dy)
    return dict(eP=eP, edr_en=edr_en, edr_dy=edr_dy, nQ=nQ,
                coup=coup, Hqq=Hqq, sumv=sumv, Pdets=Pdets, cP=cP)


# ----------------------------------------------------------------------
# 5. Independent brute-force re-run (cross-check, reusing the proto core)
# ----------------------------------------------------------------------
def run_bruteforce(fix):
    """Reproduce the brute-force pipeline of stageB/gate_icpt2_full.py exactly,
    so we can cross-validate the contracted result independent of the dumped
    oracle.  Returns (eP, edr_en, edr_dy, nQ).
    """
    norb = fix['norb']; nPd = fix['nPd']
    H1, g, _ = spinorb(fix['h'], fix['eri'])
    na = nb = 2
    dets = gen_dets(norb, na, nb)
    Hfull = build_H(dets, H1, g)
    Pidx = [i for i, d in enumerate(dets) if all((so % norb) < NACT for so in d)]
    Qidx = [i for i in range(len(dets)) if i not in set(Pidx)]
    HPP = Hfull[np.ix_(Pidx, Pidx)]
    eP, cP = np.linalg.eigh(HPP)
    HQP = Hfull[np.ix_(Qidx, Pidx)]
    Hqq = np.diag(Hfull)[Qidx]
    coup = HQP @ cP[:, :nPd]
    eps = fix['eps']
    sumv = np.array([sum(eps[so % norb] for so in dets[qi] if (so % norb) >= NACT)
                     for qi in Qidx])
    inv_en = np.column_stack([1.0 / guard(eP[k] - Hqq) for k in range(nPd)])
    inv_dy = np.column_stack([1.0 / guard(-sumv) for _ in range(nPd)])
    edr_en = downfold(eP, coup, inv_en)
    edr_dy = downfold(eP, coup, inv_dy)
    return dict(eP=eP, edr_en=edr_en, edr_dy=edr_dy, nQ=len(Qidx),
                Pidx=Pidx, Qidx=Qidx, dets=dets)


# ----------------------------------------------------------------------
# 6. Asymptotic perturber-count model (informational)
# ----------------------------------------------------------------------
def count_model(nvirt, na=2, nb=2):
    """Closed-form external-Q count for na alpha + nb beta active electrons
    promoted into nvirt virtual orbitals (frozen-core window, NACT=4 active).

    Per spin sigma with n_sigma electrons: choices =
        sum_{nv=0..n_sigma} C(NACT, n_sigma-nv) * C(nvirt, nv)   (spatial)
    but virtual spin-orbitals = nvirt, active spin-orbitals = NACT.  Total =
    (alpha choices)*(beta choices) - (all-active = P count).
    """
    from math import comb

    def spin_choices(nsig):
        return sum(comb(NACT, nsig - nv) * comb(nvirt, nv)
                   for nv in range(0, nsig + 1)
                   if nsig - nv <= NACT and nv <= nvirt)

    total = spin_choices(na) * spin_choices(nb)
    p = comb(NACT, na) * comb(NACT, nb)
    return total - p


# ----------------------------------------------------------------------
# 7. Driver / gate
# ----------------------------------------------------------------------
def main():
    np.set_printoptions(precision=6, suppress=True)
    fix = read_fixture()
    norb, nPd = fix['norb'], fix['nPd']
    nvirt = norb - NACT

    print("=" * 78)
    print("QMRSF-icPT2  CONTRACTED  external-Q proto  |  H4/6-31G live window")
    print("=" * 78)
    print(f"  norb={norb}  NACT={NACT}  nvirt={nvirt}  nPd(CAS)={nPd}  (na=nb=2)")

    con = run_contracted(fix)
    bru = run_bruteforce(fix)

    # ---- bare CAS sanity (contracted CAS == brute CAS == live) ----
    dcas_live = np.max(np.abs(np.sort(con['eP']) - np.sort(fix['eP_live'])))
    dcas_bru = np.max(np.abs(np.sort(con['eP']) - np.sort(bru['eP'])))

    # ---- GATE vs dumped oracle ----
    den_live = np.max(np.abs(np.sort(con['edr_en']) - np.sort(fix['en_live'])))
    ddy_live = np.max(np.abs(np.sort(con['edr_dy']) - np.sort(fix['dy_live'])))
    # ---- cross-check vs independent brute-force re-run ----
    den_bru = np.max(np.abs(np.sort(con['edr_en']) - np.sort(bru['edr_en'])))
    ddy_bru = np.max(np.abs(np.sort(con['edr_dy']) - np.sort(bru['edr_dy'])))

    TOL = 1e-9
    print("\n  -- perturber count --")
    print(f"     brute-force nQ (full-det partition)      = {bru['nQ']}")
    print(f"     contracted  nQ (virtual x active blocks) = {con['nQ']}")
    print(f"     closed-form count_model(nvirt={nvirt})         = {count_model(nvirt)}")
    print(f"     match: {'YES' if con['nQ'] == bru['nQ'] == count_model(nvirt) else 'NO'}")

    print("\n  -- bare CAS backbone --")
    print(f"     contracted vs live   max|d| = {dcas_live:.3e}  "
          f"-> {'PASS' if dcas_live < 1e-9 else 'FAIL'}")
    print(f"     contracted vs brute  max|d| = {dcas_bru:.3e}  "
          f"-> {'PASS' if dcas_bru < 1e-9 else 'FAIL'}")

    print("\n  -- icPT2 EN dressed spectrum --")
    print(f"     contracted vs live   max|d| = {den_live:.3e}  "
          f"-> {'PASS' if den_live < TOL else 'FAIL'}")
    print(f"     contracted vs brute  max|d| = {den_bru:.3e}  "
          f"-> {'PASS' if den_bru < TOL else 'FAIL'}")

    print("\n  -- icPT2 Dyall dressed spectrum --")
    print(f"     contracted vs live   max|d| = {ddy_live:.3e}  "
          f"-> {'PASS' if ddy_live < TOL else 'FAIL'}")
    print(f"     contracted vs brute  max|d| = {ddy_bru:.3e}  "
          f"-> {'PASS' if ddy_bru < TOL else 'FAIL'}")

    ecore = fix['ecore']
    print("\n  -- ground-state totals (electronic + ecore) --")
    print(f"     CAS   = {con['eP'].min() + ecore:.8f}")
    print(f"     EN    = {con['edr_en'].min() + ecore:.8f}")
    print(f"     Dyall = {con['edr_dy'].min() + ecore:.8f}")

    ok = (dcas_live < 1e-9 and den_live < TOL and ddy_live < TOL and
          den_bru < TOL and ddy_bru < TOL and
          con['nQ'] == bru['nQ'] == count_model(nvirt))
    print("\n" + "=" * 78)
    print("  RESULT:", "PASS  (contracted EN + Dyall reproduce oracle AND brute force)"
          if ok else "FAIL")
    print("=" * 78)

    # ---- asymptotic scaling table ----
    print("\n  -- perturber-count scaling in nvirt (na=nb=2, NACT=4) --")
    print(f"     {'nvirt':>6} {'nQ_contracted':>14} {'~nvirt^4 ref':>14}")
    for nv in (1, 2, 4, 8, 16, 32, 64):
        print(f"     {nv:>6} {count_model(nv):>14} {nv**4:>14}")
    print("     (leading term: the 2-electrons-per-spin-in-2-distinct-virtuals")
    print("      sector ~ C(nvirt,2)^2 ~ nvirt^4 / 4  -> quartic, as expected.)")

    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
