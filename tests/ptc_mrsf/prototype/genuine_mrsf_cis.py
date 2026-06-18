"""
Genuine MRSF-CIS: ROHF triplet reference + full-space spin-flip CIS.

This corrects the (2,2) active-space proxy used in the other prototypes. The real
method is:
  * reference  = ROHF high-spin (triplet) determinant -- NOT RHF, NOT CASSCF(2,2);
  * working space = single spin-flip excitations a^dag_{p beta} a_{q alpha} over
    the FULL orbital space (all occupied alpha -> all virtual beta), CIS-cost;
  * mixed reference (Ms = +1 and Ms = -1) restores spin completeness.

Demonstrated on H2 / 6-31G with the triplet as the ROHF reference (a 2-electron
case where full FCI is the exact benchmark):

  * SF-CIS (single reference) and MRSF-CIS span larger spaces than the dim-4
    (2,2) frontier: they include single excitations to the higher virtuals;
  * SF-CIS is spin-contaminated; the mixed reference makes MRSF-CIS spin-pure;
  * MRSF-CIS is closer to FCI than SF-CIS, but as a CIS-level method still misses
    dynamic correlation -- which the transcorrelation (built on the ROHF
    reference, over the full space) then recovers toward FCI.

Run:  python3 genuine_mrsf_cis.py
"""

import numpy as np
from scipy.linalg import expm
from pyscf import gto, scf, ao2mo, fci

from tc_finite_basis import ann, cre, ms0_determinants, build_fci_hamiltonian
from ptc_mrsf_cis import build_s2, spin_of
from nonsym_tda_eig import nonsym_tda_eig


def main():
    mol = gto.M(atom='H 0 0 0; H 0 0 1.2', basis='6-31g', spin=2, verbose=0)
    mfT = scf.ROHF(mol).run()                 # ROHF TRIPLET reference
    norb = mfT.mo_coeff.shape[1]
    C = mfT.mo_coeff
    h1 = C.T @ mfT.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    moe = mfT.mo_energy

    basis = ms0_determinants(norb, 1)          # full Ms=0 space (1 alpha, 1 beta)
    index = {d: i for i, d in enumerate(basis)}
    H = build_fci_hamiltonian(h1, eri, mol.energy_nuc(), norb, basis, index)
    S2 = build_s2(basis, norb, 1)
    wfci = np.linalg.eigvalsh(H)

    occ_a = [0, 1]                             # triplet singly-occupied (alpha) MOs
    allorb = list(range(norb))

    def didx(a_orbs, b_orbs):
        d = tuple(sorted([2 * p for p in a_orbs] + [2 * p + 1 for p in b_orbs]))
        return index.get(d)

    # SF-CIS: single spin-flips from the Ms=+1 triplet reference
    sf = set()
    for q in occ_a:
        rem = [o for o in occ_a if o != q]
        for p in allorb:
            i = didx(rem, [p])
            if i is not None:
                sf.add(i)
    sf = sorted(sf)
    # MRSF-CIS: also the Ms=-1 reference flips (mixed reference)
    mr = set(sf)
    for q in occ_a:
        rem = [o for o in occ_a if o != q]
        for p in allorb:
            i = didx([p], rem)
            if i is not None:
                mr.add(i)
    mr = sorted(mr)

    def states(sub_idx):
        sub = np.ix_(sub_idx, sub_idx)
        w, v = np.linalg.eigh(H[sub])
        ss = [spin_of(v[:, k], S2[sub]) for k in range(len(w))]
        return w, ss

    wsf, ssf = states(sf)
    wmr, smr = states(mr)

    print("=== genuine MRSF-CIS: ROHF triplet reference, full-space SF ===")
    print(f"ROHF triplet reference energy: {mfT.e_tot:.6f}\n")
    print(f"(2,2) frontier proxy would be dim 4)")
    print(f"SF-CIS   dim {len(sf):2d}   MRSF-CIS dim {len(mr):2d}   "
          f"full Ms=0 dim {len(basis)}\n")
    print(" state    SF-CIS  <S^2>      MRSF-CIS  <S^2>       FCI")
    for k in range(4):
        print(f"   {k}    {wsf[k]:9.5f} {ssf[k]:5.2f}    "
              f"{wmr[k]:9.5f} {smr[k]:5.2f}    {wfci[k]:9.5f}")

    # spin purity: SF-CIS contaminated, MRSF-CIS pure
    sf_contam = max(min(abs(s - 0), abs(s - 2)) for s in ssf[:4])
    mr_contam = max(min(abs(s - 0), abs(s - 2)) for s in smr[:4])
    print(f"\nmax spin contamination:  SF-CIS {sf_contam:.3f}   "
          f"MRSF-CIS {mr_contam:.3f}")

    # transcorrelation on the ROHF reference: active-pair -> external, downfolded
    eri_p = np.einsum('ikjl->ijkl', eri)  # not used; keep MO eri for amplitudes
    nso = 2 * norb
    occ_so = [2 * p for p in occ_a] + [2 * p + 1 for p in occ_a]  # active spin-orb
    ext = [P for P in range(nso) if P // 2 not in occ_a]
    esp = np.array([moe[P // 2] for P in range(nso)])

    def aphys(P, Q, R, S):
        p, sp = P // 2, P % 2
        q, sq = Q // 2, Q % 2
        r, sr = R // 2, R % 2
        s, ss = S // 2, S % 2
        v = 0.0
        if sp == sr and sq == ss:
            v += eri[p, r, q, s]
        if sp == ss and sq == sr:
            v -= eri[p, s, q, r]
        return v

    T2 = np.zeros((len(basis), len(basis)))
    for col, det in enumerate(basis):
        for ii, I in enumerate(occ_so):
            for J in occ_so[ii + 1:]:
                for ia, A in enumerate(ext):
                    for B in ext[ia + 1:]:
                        num = aphys(I, J, A, B)
                        if abs(num) < 1e-12:
                            continue
                        t = num / (esp[I] + esp[J] - esp[A] - esp[B])
                        g, d = ann(det, I)
                        if g == 0:
                            continue
                        g2, d = ann(d, J)
                        if g2 == 0:
                            continue
                        g *= g2
                        g2, d = cre(d, B)
                        if g2 == 0:
                            continue
                        g *= g2
                        g2, d = cre(d, A)
                        if g2 == 0:
                            continue
                        T2[index[d], col] += t * g * g2
    Hbar = expm(-T2) @ H @ expm(T2)
    sub = np.ix_(mr, mr)
    ee, vr, vl, info = nonsym_tda_eig(Hbar[sub], len(mr))
    ee = np.sort(ee)
    print(f"\npTC-MRSF-CIS (transcorrelated on ROHF ref): ground {ee[0]:.6f}")
    print(f"  MRSF-CIS {wmr[0]:.6f} -> pTC {ee[0]:.6f} -> FCI {wfci[0]:.6f}")
    rec = (ee[0] - wmr[0]) / (wfci[0] - wmr[0]) * 100
    print(f"  dynamic correlation recovered: {rec:.0f}%")

    # ---- validations ----
    assert len(sf) > 4 and len(mr) > 4, "full-space SF is larger than (2,2)"
    assert sf_contam > 0.01, "SF-CIS should be spin-contaminated"
    assert mr_contam < 1e-3, "MRSF-CIS (mixed ref) should be spin-pure"
    assert wmr[0] < wsf[0], "MRSF-CIS lower (better) than SF-CIS"
    assert wmr[0] > wfci[0] > -1e9
    assert ee[0] < wmr[0] and info["n_complex"] == 0, "pTC lowers toward FCI"
    print("\nVALIDATED: genuine MRSF-CIS is ROHF-referenced and full-space (not")
    print("(2,2)); the mixed reference removes SF-CIS spin contamination; and the")
    print("transcorrelation on the ROHF reference recovers dynamic correlation.")


if __name__ == "__main__":
    main()
