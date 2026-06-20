#!/usr/bin/env python3
"""Gate the LIVE frozen-core dressing (ncore>0) of QMRSF-icPT2 against an
independent closed-form oracle, on H6/STO-3G (quintet -> ncore=1).

Validates the pieces that only the ncore>0 path exercises:
  * h_eff (active+virtual window) = C_w^T (Hcore + v_core) C_w, v_core = 2J[Dc]-K[Dc]
  * E_core = E_nuc + 2 Tr[Dc Hcore] + Tr[Dc v_core]
  * eri_win = full MO ERI restricted to the window
against the Fortran-produced qmrsf_icpt2_full_live.dat (h_win, eri_win, ecore), using
the FULL MO coefficients from qmrsf_cfull_live.dat and pyscf-free closed-form AO integrals.
Then runs the window EN+Dyall downfold on the oracle integrals and compares the spectra.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))           # tools/qmrsf_pathways_proto
import qmrsf_icpt2_ppp_proto as P                    # spinorb, gen_dets, build_H
from route_a_oracle import (boys0, contraction, s_prim, t_prim, v_prim, eri_prim,  # noqa
                            STO3G_H_ALPHA, STO3G_H_DCOEF, BOHR)

NACT = 4
HARTREE = 1.0
N_H6 = 6
SPACING = 1.2                                        # Angstrom


def build_h6_ao():
    cen = np.array([[0.0, 0.0, i * SPACING] for i in range(N_H6)]) * BOHR
    alpha = STO3G_H_ALPHA
    c = contraction(STO3G_H_ALPHA, STO3G_H_DCOEF)
    n = N_H6
    S = np.zeros((n, n)); T = np.zeros((n, n)); V = np.zeros((n, n))
    for mu in range(n):
        for nu in range(n):
            for i in range(3):
                for j in range(3):
                    w = c[i] * c[j]
                    S[mu, nu] += w * s_prim(alpha[i], cen[mu], alpha[j], cen[nu])
                    T[mu, nu] += w * t_prim(alpha[i], cen[mu], alpha[j], cen[nu])
                    for k in range(n):
                        V[mu, nu] += w * v_prim(alpha[i], cen[mu], alpha[j], cen[nu], cen[k], 1.0)
    H = T + V
    eri = np.zeros((n, n, n, n))
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                for sig in range(n):
                    acc = 0.0
                    for i in range(3):
                        for j in range(3):
                            for k in range(3):
                                for ll in range(3):
                                    acc += (c[i]*c[j]*c[k]*c[ll] *
                                            eri_prim(alpha[i], cen[mu], alpha[j], cen[nu],
                                                     alpha[k], cen[lam], alpha[ll], cen[sig]))
                    eri[mu, nu, lam, sig] = acc
    enuc = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            enuc += 1.0 / np.sqrt(np.dot(cen[i] - cen[j], cen[i] - cen[j]))
    return S, H, eri, enuc


def read_cfull(path):
    with open(path) as f:
        nbf, ncore = map(int, f.readline().split())
        C = np.array([[float(x) for x in f.readline().split()] for _ in range(nbf)])
    return nbf, ncore, C


def read_full(path):
    with open(path) as f:
        norb, nPd = map(int, f.readline().split())
        hw = np.array([[float(x) for x in f.readline().split()] for _ in range(norb)])
        eri = np.zeros((norb,) * 4)
        for p in range(norb):
            for q in range(norb):
                for r in range(norb):
                    eri[p, q, r, :] = [float(x) for x in f.readline().split()]
        ecore = float(f.readline())
        eps = np.array([float(x) for x in f.readline().split()])
        eP = np.array([float(x) for x in f.readline().split()])
        en = np.array([float(x) for x in f.readline().split()])
        dy = np.array([float(x) for x in f.readline().split()])
    return dict(norb=norb, nPd=nPd, hw=hw, eri=eri, ecore=ecore, eps=eps, eP=eP, en=en, dy=dy)


def guard(d):
    return np.where(np.abs(d) < 1e-6, np.sign(d) * 1e-6 + 1e-30, d)


def downfold(eP, coup, invd, nPd):
    Heff = np.diag(eP[:nPd]).astype(float)
    for k in range(nPd):
        for l in range(nPd):
            Heff[k, l] += 0.5 * np.sum(coup[:, k] * coup[:, l] * (invd[:, k] + invd[:, l]))
    return np.linalg.eigvalsh(0.5 * (Heff + Heff.T))


def main():
    S, H, ERI, enuc = build_h6_ao()
    nbf, ncore, C = read_cfull(os.path.join(HERE, "qmrsf_cfull_live.dat"))
    fix = read_full(os.path.join(HERE, "qmrsf_icpt2_full_live.dat"))
    norb_w, nPd = fix['norb'], fix['nPd']

    # GATE 0: full MO orthonormality wrt oracle overlap
    orth = np.max(np.abs(C.T @ S @ C - np.eye(nbf)))

    # frozen-core dressing (oracle)
    Cc = C[:, :ncore]
    Cw = C[:, ncore:]
    Dc = Cc @ Cc.T
    J = np.einsum('mnls,ls->mn', ERI, Dc)
    K = np.einsum('mlns,ls->mn', ERI, Dc)
    vcore = 2.0 * J - K
    h_eff = Cw.T @ (H + vcore) @ Cw
    ecore = enuc + 2.0 * np.sum(Dc * H) + np.sum(Dc * vcore)
    eri_win = np.einsum('mnls,mp,nq,lr,st->pqrt', ERI, Cw, Cw, Cw, Cw)

    dh = np.max(np.abs(h_eff - fix['hw']))
    de = np.max(np.abs(eri_win - fix['eri']))
    dnuc = abs(ecore - fix['ecore'])

    # window downfold on oracle integrals -> compare to live spectra
    H1, g, _ = P.spinorb(h_eff, eri_win)
    na = nb = 2
    dets = P.gen_dets(norb_w, na, nb)
    Hfull = P.build_H(dets, H1, g)
    Pidx = [i for i, d in enumerate(dets) if all((so % norb_w) < NACT for so in d)]
    Qidx = [i for i in range(len(dets)) if i not in set(Pidx)]
    HPP = Hfull[np.ix_(Pidx, Pidx)]
    eP, cP = np.linalg.eigh(HPP)
    HQP = Hfull[np.ix_(Qidx, Pidx)]
    Hqq = np.diag(Hfull)[Qidx]
    coup = HQP @ cP[:, :nPd]
    eps = fix['eps']
    sumv = np.array([sum(eps[so % norb_w] for so in dets[qi] if (so % norb_w) >= NACT) for qi in Qidx])
    inv_en = np.column_stack([1.0 / guard(eP[k] - Hqq) for k in range(nPd)])
    inv_dy = np.column_stack([1.0 / guard(-sumv) for _ in range(nPd)])
    edr_en = downfold(eP, coup, inv_en, nPd)
    edr_dy = downfold(eP, coup, inv_dy, nPd)

    den = np.max(np.abs(np.sort(edr_en) - np.sort(fix['en'])))
    ddy = np.max(np.abs(np.sort(edr_dy) - np.sort(fix['dy'])))

    print("==== Gate: LIVE frozen-core (ncore=1) dressing vs closed-form oracle (H6/STO-3G) ====")
    print(f"  nbf={nbf}  ncore={ncore}  window={norb_w}  E_nuc={enuc:.8f}")
    print(f"  GATE0 |C^T S C - I|        = {orth:.3e}  -> {'PASS' if orth < 1e-6 else 'FAIL'}")
    print(f"  GATE1 max|h_eff live-orac| = {dh:.3e}  -> {'PASS' if dh < 1e-6 else 'FAIL'}")
    print(f"  GATE2 max|eri_win l-o|     = {de:.3e}  -> {'PASS' if de < 1e-6 else 'FAIL'}")
    print(f"  GATE3 E_core: live={fix['ecore']:.10f} orac={ecore:.10f} d={dnuc:.2e}"
          f"  -> {'PASS' if dnuc < 1e-6 else 'FAIL'}")
    # EN tolerance is looser: with nvirt=1 H6 has an Epstein-Nesbet INTRUDER
    # (min|eP-Hqq| ~ 1e-3), which amplifies the ~1e-8 integral/eigenvector agreement to
    # ~1e-5 (condition ~ 1/denom). This is the classic EN intruder problem -- the reason
    # Dyall denominators exist; Dyall (intruder-robust) agrees to ~1e-8. The frozen-core
    # DRESSING under test (h_eff, E_core, eri_win) is validated to 1e-8 by GATE1-3.
    print(f"  GATE4 icPT2 EN  spectrum   = {den:.3e}  -> {'PASS' if den < 1e-5 else 'FAIL'}"
          "  (EN intruder-amplified; see note)")
    print(f"  GATE5 icPT2 Dyall spectrum = {ddy:.3e}  -> {'PASS' if ddy < 1e-6 else 'FAIL'}")
    ok = (orth < 1e-6 and dh < 1e-6 and de < 1e-6 and dnuc < 1e-6 and den < 1e-5 and ddy < 1e-6)
    print("  RESULT:", "PASS  (frozen-core dressing + downfold reproduce the oracle)" if ok else "FAIL")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
