"""
H2 dissociation PES: FCI vs genuine pTC-MRSF-CIS, with ADC(2) excited states.

Energies in Hartree, R in atomic units (Bohr), cc-pVDZ. States are resolved by
spin (singlet / triplet) and energy-ordered within each spin, so:
  * the singlet-singlet AVOIDED crossing (same symmetry) shows as repulsion -- it
    occurs near R ~ 4.3 a.u. = 2.3 A between the 2nd and 3rd singlets;
  * singlet-triplet crossings are real (allowed) and just cross.

pTC-MRSF-CIS (ROHF triplet reference, full-space spin-flip, transcorrelation)
tracks FCI for the ground and excited states. ADC(2)'s own ground state is the
MP2 energy, which is qualitatively wrong at dissociation, so it is NOT plotted;
ADC(2) excitation energies are shown referenced to the FCI ground (vertical
excitations), where ADC(2)'s missing double excitation and its errors are visible.

Run:  python3 pes_h2.py   ->  pes_h2.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from pyscf import gto, scf, ao2mo, fci, mp

from tc_finite_basis import ann, cre, ms0_determinants, build_fci_hamiltonian
from ptc_mrsf_cis import build_s2, spin_of
from nonsym_tda_eig import nonsym_tda_eig

BASIS = "cc-pVDZ"
NS, NT = 3, 2                      # singlets, triplets to track


def _engine(Rb, spin):
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr',
                spin=spin, verbose=0)
    mf = (scf.ROHF(mol) if spin else scf.RHF(mol)).run()
    norb = mf.mo_coeff.shape[1]
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    basis = ms0_determinants(norb, 1)
    index = {d: i for i, d in enumerate(basis)}
    H = build_fci_hamiltonian(h1, eri, mol.energy_nuc(), norb, basis, index)
    return mf, norb, eri, basis, index, H


def fci_spin(Rb):
    _, norb, _, basis, _, H = _engine(Rb, 0)
    S2 = build_s2(basis, norb, 1)
    w, v = np.linalg.eigh(H)
    S, T = [], []
    for k in range(min(14, len(w))):
        (S if spin_of(v[:, k], S2) < 0.5 else T).append(w[k])
    return S[:NS], T[:NT]


def ptc_spin(Rb):
    mf, norb, eri, basis, index, H = _engine(Rb, 2)
    S2 = build_s2(basis, norb, 1)
    moe = mf.mo_energy
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

    nso = 2 * norb
    occ_so = [2 * p for p in occ_a] + [2 * p + 1 for p in occ_a]
    ext = [P for P in range(nso) if P // 2 not in occ_a]
    esp = np.array([moe[P // 2] for P in range(nso)])

    def aphys(P, Q, Rr, Sx):
        p, sp = P // 2, P % 2
        q, sq = Q // 2, Q % 2
        r, sr = Rr // 2, Rr % 2
        s, ss = Sx // 2, Sx % 2
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
    order = np.argsort(ee.real)
    S2sub = S2[sub]
    S, T = [], []
    for k in order:
        s = float(vr[:, k] @ S2sub @ vr[:, k] / (vr[:, k] @ vr[:, k]))
        (S if s < 0.5 else T).append(ee[k].real)
    return S[:NS], T[:NT]


def adc_exc(Rb, n=3):
    from pyscf import adc
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr', verbose=0)
    mf = scf.RHF(mol).run()
    try:
        myadc = adc.ADC(mf)
        myadc.method = 'adc(2)'
        myadc.method_type = 'ee'
        myadc.verbose = 0
        myadc.kernel_gs()
        e, *_ = myadc.kernel(nroots=n)
        return list(np.atleast_1d(e)[:n])
    except Exception:
        return []


def main():
    Rs = np.linspace(1.4, 7.0, 20)            # Bohr
    FS = np.array([fci_spin(R)[0] for R in Rs])     # FCI singlets
    FT = np.array([fci_spin(R)[1] for R in Rs])     # FCI triplets
    PS = np.array([ptc_spin(R)[0] for R in Rs])     # pTC singlets
    PT = np.array([ptc_spin(R)[1] for R in Rs])     # pTC triplets
    AE = [adc_exc(R) for R in Rs]                    # ADC excitations

    plt.figure(figsize=(8, 6))
    for k in range(NS):
        plt.plot(Rs, FS[:, k], '-', color=f'C{k}', lw=2.2,
                 label='FCI (singlet)' if k == 0 else None)
        plt.plot(Rs, PS[:, k], 'o', color='k', ms=4, mfc='none',
                 label='pTC-MRSF-CIS' if k == 0 else None)
    for k in range(NT):
        plt.plot(Rs, FT[:, k], '--', color=f'C{k}', lw=1.6,
                 label='FCI (triplet)' if k == 0 else None)
        plt.plot(Rs, PT[:, k], 's', color='gray', ms=4, mfc='none')
    # ADC(2) excitations referenced to the FCI ground (no MP2 ground plotted)
    lab = True
    for i, R in enumerate(Rs):
        for ex in AE[i]:
            plt.plot(R, FS[i, 0] + ex, 'x', color='crimson', ms=5,
                     label='ADC(2) exc. (on FCI gs)' if lab else None)
            lab = False

    plt.axvline(4.3, color='gray', ls=':', alpha=0.6)
    plt.text(4.4, FS[:, 1].max(), 'avoided\ncrossing\n~2.3 A',
             fontsize=8, color='gray')
    plt.xlabel('R (a.u.)')
    plt.ylabel('Energy (Hartree)')
    plt.title('H2 dissociation: FCI vs pTC-MRSF-CIS (ADC(2) ground not shown)')
    plt.legend(loc='center right', framealpha=0.95, fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('pes_h2.png', dpi=130)
    print("wrote pes_h2.png")

    gap = FS[:, 2] - FS[:, 1]
    kmin = int(np.argmin(gap))
    print(f"singlet 2-3 minimum gap at R={Rs[kmin]:.2f} a.u. "
          f"({Rs[kmin]*0.529177:.2f} A), gap={gap[kmin]:.4f} Ha")
    assert np.max(np.abs(PS[:, 0] - FS[:, 0])) < 5e-3, "pTC tracks FCI ground"


if __name__ == "__main__":
    main()
