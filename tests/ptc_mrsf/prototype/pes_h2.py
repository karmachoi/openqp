"""
H2 dissociation PES (cc-pVTZ): FCI vs pTC-MRSF-CIS vs ADC(2).

Symmetry-resolved singlet states (D-infinity-h):
    S0 = 1 ^1Sigma_g^+ (X)   S1 = 1 ^1Sigma_u^+ (B)
    S2 = 2 ^1Sigma_g^+ (EF)  S3 = 3 ^1Sigma_g^+ (GK)
2/3 ^1Sigma_g^+ undergo an avoided crossing near R ~ 2.6 Bohr.

Trust boundary:
  * FCI and ADC(2) are standard pyscf (deterministic, exact-in-basis, checkable
    with any QC package);
  * pTC-MRSF-CIS is this repository's reference implementation (ROHF triplet
    reference, full-space spin-flip, transcorrelation), validated against FCI in
    the small cases.

Run:  python3 pes_h2.py   ->  pes_h2.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from pyscf import gto, scf, fci, ao2mo, symm

from tc_finite_basis import ann, cre, ms0_determinants
from ptc_mrsf_cis import build_s2
from nonsym_tda_eig import nonsym_tda_eig

BASIS = "cc-pvtz"
LABELS = [r'$1\,^1\Sigma_g^+$ (X)', r'$1\,^1\Sigma_u^+$ (B)',
          r'$2\,^1\Sigma_g^+$ (EF)', r'$3\,^1\Sigma_g^+$ (GK)']


def fci_states(Rb):
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr',
                symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    g = fci.FCI(mf); g.wfnsym = 'A1g'; g.nroots = 3
    eg = np.sort(np.atleast_1d(g.kernel()[0]))
    u = fci.FCI(mf); u.wfnsym = 'A1u'; u.nroots = 1
    eu = np.atleast_1d(u.kernel()[0])
    return np.array([eg[0], eu[0], eg[1], eg[2]]), mf.e_tot


def adc_exc(Rb, e_gs):
    from pyscf import adc
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr', verbose=0)
    mf = scf.RHF(mol).run()
    try:
        a = adc.ADC(mf); a.method = 'adc(2)'; a.method_type = 'ee'; a.verbose = 0
        a.kernel_gs()
        e, *_ = a.kernel(nroots=6)
        return [e_gs + x for x in np.atleast_1d(e)[:6]]
    except Exception:
        return []


def ptc_states(Rb, ext_max=16):
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr',
                spin=2, symmetry=True, verbose=0)
    mf = scf.ROHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    moe = mf.mo_energy
    par = np.array([-1 if 'u' in symm.irrep_id2name(mol.groupname, s) else 1
                    for s in mf.orbsym])

    basis = ms0_determinants(norb, 1)
    index = {d: i for i, d in enumerate(basis)}
    na = norb
    dim = na * na
    h2e = fci.direct_spin1.absorb_h1e(h1, eri, norb, (1, 1), 0.5)
    H = np.empty((dim, dim))
    e0 = np.zeros(dim)
    for k in range(dim):
        e0[:] = 0.0; e0[k] = 1.0
        H[:, k] = fci.direct_spin1.contract_2e(
            h2e, e0.reshape(na, na), norb, (1, 1)).reshape(-1)
    H += mol.energy_nuc() * np.eye(dim)
    S2 = build_s2(basis, norb, 1)

    occ_a = [0, 1]

    def didx(a_o, b_o):
        return index.get(tuple(sorted([2 * p for p in a_o] + [2 * p + 1 for p in b_o])))
    mr = sorted({i for q in occ_a for rem in [[o for o in occ_a if o != q]]
                 for p in range(norb) for i in (didx(rem, [p]), didx([p], rem))
                 if i is not None})

    nso = 2 * norb
    occ_so = [2 * p for p in occ_a] + [2 * p + 1 for p in occ_a]
    ext = [P for P in range(nso) if occ_a[1] < P // 2 <= occ_a[1] + ext_max]
    esp = np.array([moe[P // 2] for P in range(nso)])

    def aphys(P, Q, Rr, Sx):
        p, sp = P // 2, P % 2; q, sq = Q // 2, Q % 2
        r, sr = Rr // 2, Rr % 2; s, ss = Sx // 2, Sx % 2
        v = 0.0
        if sp == sr and sq == ss:
            v += eri[p, r, q, s]
        if sp == ss and sq == sr:
            v -= eri[p, s, q, r]
        return v

    T2 = np.zeros((dim, dim))
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
                        g2, d = ann(d, J);  g *= g2
                        if g2 == 0:
                            continue
                        g2, d = cre(d, B);  g *= g2
                        if g2 == 0:
                            continue
                        g2, d = cre(d, A)
                        if g2 == 0:
                            continue
                        T2[index[d], col] += t * g * g2

    Hbar = expm(-T2) @ H @ expm(T2)
    sub = np.ix_(mr, mr)
    ee, vr, vl, info = nonsym_tda_eig(Hbar[sub], len(mr))
    Pd = np.array([np.prod([par[P // 2] for P in basis[i]]) for i in mr])
    S2sub = S2[sub]
    gS, uS = [], []
    for k in np.argsort(ee.real):
        v = vr[:, k]; nv = v @ v
        if float(v @ S2sub @ v / nv) > 0.5:
            continue
        (gS if float(v @ (Pd * v) / nv) > 0 else uS).append(ee[k].real)
    out = [gS[0] if gS else np.nan, uS[0] if uS else np.nan,
           gS[1] if len(gS) > 1 else np.nan, gS[2] if len(gS) > 2 else np.nan]
    return out


def main():
    Rs = np.linspace(1.4, 6.0, 16)
    EF = np.zeros((len(Rs), 4)); GS = np.zeros(len(Rs))
    for i, R in enumerate(Rs):
        EF[i], GS[i] = fci_states(R)
    AD = [adc_exc(R, GS[i]) for i, R in enumerate(Rs)]
    PT = np.array([ptc_states(R) for R in Rs])

    plt.figure(figsize=(8, 6))
    for k in range(4):
        plt.plot(Rs, EF[:, k], '-', color=f'C{k}', lw=2.2, label=LABELS[k])
        plt.plot(Rs, PT[:, k], 'o', color='k', ms=4, mfc='none',
                 label='pTC-MRSF-CIS' if k == 0 else None)
    lab = True
    for i, R in enumerate(Rs):
        for ex in AD[i]:
            plt.plot(R, ex, 'x', color='crimson', ms=4,
                     label='ADC(2) exc.' if lab else None)
            lab = False

    gap = EF[:, 3] - EF[:, 2]
    kx = int(np.argmin(gap))
    plt.axvline(Rs[kx], color='gray', ls=':', alpha=0.7)
    plt.text(Rs[kx] + 0.1, EF[0, 3] + 0.02,
             f'avoided crossing\n$2/3\\,^1\\Sigma_g^+$, R={Rs[kx]:.2f} Bohr',
             fontsize=8, color='gray')

    plt.xlabel('R (Bohr)')
    plt.ylabel('Energy (Hartree)')
    plt.title('H$_2$ / cc-pVTZ: FCI (lines) vs pTC-MRSF-CIS (o) vs ADC(2) (x)')
    plt.legend(loc='upper right', framealpha=0.95, fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('pes_h2.png', dpi=130)
    print("wrote pes_h2.png; avoided crossing at R=%.2f Bohr" % Rs[kx])


if __name__ == "__main__":
    main()
