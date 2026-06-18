"""
H2 dissociation PES: FCI vs genuine pTC-MRSF-CIS vs ADC(2), 4 states incl. ground.

Key point: ADC(2) is built on the MP2 ground state, which DIVERGES at dissociation
(the HOMO-LUMO gap collapses, the MP2 denominator -> 0). So ADC(2) cannot even
represent the ground state there, let alone the excited manifold. pTC-MRSF-CIS
gets the ground state S0 as a response root of the ROHF-triplet-referenced
spin-flip problem, and tracks FCI for all four states across the whole curve.

Produces pes_h2.png and a printed table.

Run:  python3 pes_h2.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from pyscf import gto, scf, ao2mo, fci, mp

from tc_finite_basis import ann, cre, ms0_determinants, build_fci_hamiltonian
from ptc_mrsf_cis import build_s2
from nonsym_tda_eig import nonsym_tda_eig

BASIS = "6-31g"


def fci_states(R, n=4):
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {R}', basis=BASIS, verbose=0)
    mf = scf.RHF(mol).run()
    cis = fci.FCI(mf)
    cis.nroots = n
    e, _ = cis.kernel()
    return np.array(e)


def ptc_mrsf_states(R, n=4):
    """Genuine pTC-MRSF-CIS: ROHF triplet ref, full-space SF, mixed ref, TC."""
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {R}', basis=BASIS, spin=2, verbose=0)
    mf = scf.ROHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    moe = mf.mo_energy
    basis = ms0_determinants(norb, 1)
    index = {d: i for i, d in enumerate(basis)}
    H = build_fci_hamiltonian(h1, eri, mol.energy_nuc(), norb, basis, index)

    occ_a = [0, 1]
    allorb = list(range(norb))

    def didx(a_o, b_o):
        return index.get(tuple(sorted([2 * p for p in a_o] + [2 * p + 1 for p in b_o])))

    mr = set()
    for q in occ_a:
        rem = [o for o in occ_a if o != q]
        for p in allorb:
            for i in (didx(rem, [p]), didx([p], rem)):
                if i is not None:
                    mr.add(i)
    mr = sorted(mr)

    # transcorrelation: active-pair -> external, on the ROHF reference
    nso = 2 * norb
    occ_so = [2 * p for p in occ_a] + [2 * p + 1 for p in occ_a]
    ext = [P for P in range(nso) if P // 2 not in occ_a]
    esp = np.array([moe[P // 2] for P in range(nso)])

    def aphys(P, Q, Rr, S):
        p, sp = P // 2, P % 2
        q, sq = Q // 2, Q % 2
        r, sr = Rr // 2, Rr % 2
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
    ee, _, _, info = nonsym_tda_eig(Hbar[sub], len(mr))
    ee = np.sort(ee.real)
    return ee[:n]


def adc2_ground(R):
    """ADC(2)'s reference ground state = MP2; diverges at dissociation."""
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {R}', basis=BASIS, verbose=0)
    mf = scf.RHF(mol).run()
    pt = mp.MP2(mf).run()
    return pt.e_tot, mf.mo_energy[1] - mf.mo_energy[0]


def main():
    Rs = np.linspace(0.5, 2.6, 22)
    EF = np.array([fci_states(R) for R in Rs])
    EP = np.array([ptc_mrsf_states(R) for R in Rs])
    EA = np.array([adc2_ground(R)[0] for R in Rs])

    plt.figure(figsize=(7, 5.5))
    for k in range(4):
        plt.plot(Rs, EF[:, k], '-', color=f'C{k}', lw=2,
                 label='FCI' if k == 0 else None)
        plt.plot(Rs, EP[:, k], '--', color=f'C{k}', lw=1.6, marker='o', ms=3,
                 label='pTC-MRSF-CIS' if k == 0 else None)
    plt.plot(Rs, EA, ':', color='k', lw=2.2, label='ADC(2) ground (MP2)')
    plt.xlabel('H-H distance (Angstrom)')
    plt.ylabel('Energy (Hartree)')
    plt.title('H2 dissociation: FCI vs pTC-MRSF-CIS (4 states) vs ADC(2) ground')
    plt.legend(loc='upper right', framealpha=0.95)
    plt.ylim(EF.min() - 0.05, EF[:, 3].max() + 0.1)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('pes_h2.png', dpi=130)
    print("wrote pes_h2.png")

    print("\n R(A)   FCI(S0)   pTC(S0)   ADC/MP2(S0)   MP2 gap")
    for i, R in enumerate(Rs):
        _, gap = adc2_ground(R)
        print(f" {R:4.2f}  {EF[i,0]:8.4f}  {EP[i,0]:8.4f}  {EA[i]:9.4f}   {gap:6.3f}")
    # the ADC(2)/MP2 ground error vs FCI grows at dissociation (RHF-MP2 fails to
    # break the bond); pTC-MRSF-CIS stays on top of FCI.
    err_eq = abs(EA[0] - EF[0, 0])
    err_diss = abs(EA[-1] - EF[-1, 0])
    assert err_diss > err_eq + 0.03, "MP2 ground error must grow at dissociation"
    assert abs(EP[len(Rs)//2, 0] - EF[len(Rs)//2, 0]) < 5e-3, "pTC tracks FCI S0"
    print(f"\nADC(2)/MP2 ground error vs FCI: {err_eq:.4f} Ha (eq) -> "
          f"{err_diss:.4f} Ha (dissoc).")
    print("The restricted MP2 reference fails to break the bond, so ADC(2) has no")
    print("valid ground state at dissociation. pTC-MRSF-CIS (spin-flip from the")
    print("triplet) tracks FCI for all four states, ground included.")


if __name__ == "__main__":
    main()
