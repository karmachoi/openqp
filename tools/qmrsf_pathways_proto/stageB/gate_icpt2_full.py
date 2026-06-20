#!/usr/bin/env python3
"""Gate the LIVE QMRSF-icPT2 downfold against the NumPy prototype on identical
window integrals. Reads qmrsf_icpt2_full_live.dat (norb, h, chemist eri, ecore,
live bare-CAS eP, live icPT2 edr), rebuilds the determinant CI + multistate EN
downfold with the validated prototype routines, and compares spectra.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # tools/qmrsf_pathways_proto
import qmrsf_icpt2_ppp_proto as P                   # spinorb, gen_dets, build_H

NACT = 4


def main():
    f = open(os.path.join(HERE, "qmrsf_icpt2_full_live.dat"))
    norb, nPd = map(int, f.readline().split())
    h = np.array([[float(x) for x in f.readline().split()] for _ in range(norb)])
    eri = np.zeros((norb,)*4)
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                eri[p, q, r, :] = [float(x) for x in f.readline().split()]
    ecore = float(f.readline())
    eps = np.array([float(x) for x in f.readline().split()])
    eP_live = np.array([float(x) for x in f.readline().split()])
    en_live = np.array([float(x) for x in f.readline().split()])
    dy_live = np.array([float(x) for x in f.readline().split()])
    f.close()

    H1, g, _ = P.spinorb(h, eri)
    na = nb = 2
    dets = P.gen_dets(norb, na, nb)
    Hfull = P.build_H(dets, H1, g)
    Pidx = [i for i, d in enumerate(dets) if all((so % norb) < NACT for so in d)]
    Qidx = [i for i in range(len(dets)) if i not in set(Pidx)]

    HPP = Hfull[np.ix_(Pidx, Pidx)]
    eP, cP = np.linalg.eigh(HPP)
    HQP = Hfull[np.ix_(Qidx, Pidx)]
    Hqq = np.diag(Hfull)[Qidx]
    coup = HQP @ cP[:, :nPd]

    # Dyall sum of occupied-virtual MO energies per Q determinant (spatial = so%norb)
    sumv = np.array([sum(eps[so % norb] for so in dets[qi] if (so % norb) >= NACT) for qi in Qidx])

    def downfold(invd):
        Heff = np.diag(eP[:nPd]).astype(float)
        for k in range(nPd):
            for l in range(nPd):
                Heff[k, l] += 0.5 * np.sum(coup[:, k] * coup[:, l] * (invd[:, k] + invd[:, l]))
        return np.linalg.eigvalsh(0.5 * (Heff + Heff.T))

    def guard(d):
        return np.where(np.abs(d) < 1e-6, np.sign(d) * 1e-6 + 1e-30, d)

    inv_en = np.column_stack([1.0 / guard(eP[k] - Hqq) for k in range(nPd)])
    inv_dy = np.column_stack([1.0 / guard(-sumv) for _ in range(nPd)])
    edr_en = downfold(inv_en)
    edr_dy = downfold(inv_dy)

    dcas = np.max(np.abs(np.sort(eP[:nPd]) - np.sort(eP_live)))
    den = np.max(np.abs(np.sort(edr_en) - np.sort(en_live)))
    ddy = np.max(np.abs(np.sort(edr_dy) - np.sort(dy_live)))
    print("==== Gate: live QMRSF-icPT2 downfold vs NumPy prototype (H4/6-31G window) ====")
    print(f"  norb={norb}  ndet={len(dets)}  P={len(Pidx)}  Q={len(Qidx)}  nPdress={nPd}")
    print(f"  GATE bare CAS    max|live-oracle| = {dcas:.3e}  -> {'PASS' if dcas < 1e-7 else 'FAIL'}")
    print(f"  GATE icPT2 EN    max|live-oracle| = {den:.3e}  -> {'PASS' if den < 1e-7 else 'FAIL'}")
    print(f"  GATE icPT2 Dyall max|live-oracle| = {ddy:.3e}  -> {'PASS' if ddy < 1e-7 else 'FAIL'}")
    print(f"  ground (total): CAS={eP_live.min()+ecore:.8f}  EN={en_live.min()+ecore:.8f}"
          f"  Dyall={dy_live.min()+ecore:.8f}")
    ok = dcas < 1e-7 and den < 1e-7 and ddy < 1e-7
    print("  RESULT:", "PASS  (live EN + Dyall downfold reproduce the oracle)" if ok else "FAIL")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
