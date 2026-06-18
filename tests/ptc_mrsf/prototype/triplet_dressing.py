"""
Dressing the triplet: completing the dual-cusp picture of pTC-MRSF-CIS.

In ptc_mrsf_cis.py the closed-shell MP2 correlator dressed the singlet states
(including S0) but left the triplet exactly unchanged. The reason is structural:
the ground reference correlates only the doubly-occupied sigma_g^2 pair, so the
*active* (O1,O2) pair correlation that the excited states -- and especially the
triplet -- require is simply absent.

The fix is to correlate the active pair itself to the external space, for every
spin coupling. This dresses ALL MRSF states. Crucially, the triplet receives a
SMALLER correction than the singlets, because a parallel-spin (triplet) pair has
the weaker 1/4 coalescence cusp (a Fermi hole already keeps the electrons apart)
versus the 1/2 cusp of an antiparallel (singlet) pair. That singlet-vs-triplet
asymmetry is exactly the dual cusp condition Ten-no's pTC is built to satisfy,
and which MRSF's dual singlet/triplet targets require.

We demonstrate on stretched H2 / cc-pVDZ, benchmarking every state against FCI:
  * ground-reference correlator  -> singlets dressed, triplet untouched;
  * active-pair correlator        -> all states dressed, triplet less than singlets.

Run:  python3 triplet_dressing.py
"""

import numpy as np
from scipy.linalg import expm
from pyscf import gto, scf, ao2mo

from nonsym_tda_eig import nonsym_tda_eig
from tc_finite_basis import ann, cre, ms0_determinants, build_fci_hamiltonian
from ptc_mrsf_cis import build_s2, spin_of


def antisym_phys(eri, P, Q, R, S):
    """<PQ||RS> spin-orbital from spatial chemist eri (P = 2*p + spin)."""
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


def build_t2(eri, esp, occ, vir, basis, index):
    """Spin-orbital MP2-like cluster operator from occ -> vir."""
    T = np.zeros((len(basis), len(basis)))
    for col, det in enumerate(basis):
        for ii, I in enumerate(occ):
            for J in occ[ii + 1:]:
                for ia, A in enumerate(vir):
                    for B in vir[ia + 1:]:
                        num = antisym_phys(eri, I, J, A, B)
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
                        T[index[d], col] += t * g * g2
    return T


def downfold_states(H, T2, sub, n):
    Hbar = expm(-T2) @ H @ expm(T2)
    ee, vr, vl, info = nonsym_tda_eig(Hbar[sub], n)
    assert info["n_complex"] == 0
    return np.sort(ee)


def main():
    mol = gto.M(atom='H 0 0 0; H 0 0 1.6', basis='cc-pvdz', verbose=0)
    mf = scf.RHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    nocc = mol.nelectron // 2
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    esp = np.array([mf.mo_energy[P // 2] for P in range(2 * norb)])

    basis = ms0_determinants(norb, nocc)
    index = {d: i for i, d in enumerate(basis)}
    H = build_fci_hamiltonian(h1, eri, mol.energy_nuc(), norb, basis, index)
    S2 = build_s2(basis, norb, nocc)
    wfci = np.linalg.eigvalsh(H)

    act = [nocc - 1, nocc]
    compact = sorted(i for i, d in enumerate(basis)
                     if all((P // 2) in act for P in d))
    sub = np.ix_(compact, compact)
    wb, vb = np.linalg.eigh(H[sub])
    spins = [spin_of(vb[:, k], S2[sub]) for k in range(len(compact))]

    # (A) ground-reference correlator: occupied = ground occ spin-orbitals
    occ_g = [P for P in range(2 * norb) if P // 2 < nocc]
    vir_g = [P for P in range(2 * norb) if P // 2 >= nocc]
    ee_g = downfold_states(H, build_t2(eri, esp, occ_g, vir_g, basis, index),
                           sub, len(compact))

    # (B) active-pair correlator: occupied = the two active orbitals
    occ_a = [P for P in range(2 * norb) if P // 2 in act]
    vir_a = [P for P in range(2 * norb) if P // 2 > act[1]]
    ee_a = downfold_states(H, build_t2(eri, esp, occ_a, vir_a, basis, index),
                           sub, len(compact))

    print("=== Dressing the triplet (H2 / cc-pVDZ, vs FCI) ===\n")
    print(" state spin     bare     ground-ref   active-pair      FCI")
    for k in range(len(compact)):
        lbl = 'S' if spins[k] < 0.5 else 'T'
        print(f"   {k}    {lbl}   {wb[k]:9.5f}  {ee_g[k]:9.5f}   "
              f"{ee_a[k]:9.5f}   {wfci[k]:9.5f}")

    # triplet index
    kt = next(k for k in range(len(compact)) if spins[k] > 1.5)
    ks = 0  # ground singlet
    dt_ground = abs(ee_g[kt] - wb[kt])
    dt_active = abs(ee_a[kt] - wb[kt])
    ds_active = abs(ee_a[ks] - wb[ks])
    print(f"\ntriplet shift:  ground-ref {dt_ground:.2e}  (untouched)  ->  "
          f"active-pair {dt_active:.2e}  (dressed)")
    print(f"singlet S0 shift (active-pair): {ds_active:.2e}")
    print(f"triplet/singlet shift ratio = {dt_active/ds_active:.2f}  "
          "(< 1: triplet 1/4 cusp vs singlet 1/2 cusp)")

    # ---- validations ----
    assert dt_ground < 1e-6, "ground-ref correlator must leave the triplet untouched"
    assert dt_active > 1e-4, "active-pair correlator must dress the triplet"
    for k in range(len(compact)):                       # all states improve
        assert abs(ee_a[k] - wfci[k]) <= abs(wb[k] - wfci[k]) + 1e-9
    assert dt_active < ds_active, "triplet dressed less than singlet (1/4 vs 1/2 cusp)"
    print("\nVALIDATED: the active-pair correlator dresses every MRSF state,")
    print("including the triplet, and dresses the triplet less than the singlet")
    print("-- the singlet(1/2)/triplet(1/4) dual cusp that pTC is built to satisfy.")


if __name__ == "__main__":
    main()
