"""
H2 dissociation, singlet states, cc-pVTZ. Energies relative to the S0 (X 1Sg+)
dissociation limit (-> 0). Three methods, same color per state:
    FCI  : solid     pTC-MRSF-CIS : dotted     ADC(2) : dashed

States (all singlets):
    S0 = 1 ^1Sigma_g^+ (X)            S1 = 1 ^1Sigma_u^+ (B)
    S2 = ^1Sigma_g^+ Rydberg (E/EF)   S3 = ^1Sigma_g^+ doubly-excited (sigma_u^2)
S3 is HIGH at short R (sigma_u^2 is a double excitation) and descends; it is a
multireference state that FCI and pTC-MRSF-CIS capture but ADC(2) (singles)
cannot. S3 is identified by its sigma_u^2 character, not by energy order.

FCI and ADC(2) are standard pyscf; pTC-MRSF-CIS is this repo's reference code.

Run:  python3 pes_h2.py  ->  pes_h2.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from pyscf import gto, scf, fci, ao2mo, symm

from tc_finite_basis import ann, cre, ms0_determinants
from ptc_mrsf_cis import build_s2

BASIS = "cc-pvtz"
COLORS = ['black', 'red', 'blue', 'green']     # match the reference figure
LABELS = [r'$1\,^1\Sigma_g^+$', r'$1\,^1\Sigma_u^+$',
          r'$2\,^1\Sigma_g^+$', r'$3\,^1\Sigma_g^+$']


def fci_states(Rb):
    """[S0, S1, S2, S3] singlets. S0=X(1Sg+), S1=B(1Su+ singlet),
    S2/S3 = lower/upper adiabats of the (Rydberg, sigma_u^2) 1Sg+ pair."""
    from pyscf.fci import spin_op
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr',
                symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    ne = mol.nelectron
    g = fci.FCI(mf); g.wfnsym = 'A1g'; g.nroots = 8
    eg, cg = g.kernel(); eg = np.atleast_1d(eg)
    sing = sorted((eg[k], cg[k][1, 1] ** 2) for k in range(len(eg))
                  if spin_op.spin_square(cg[k], norb, ne)[0] < 0.5)
    s0 = sing[0][0]
    exc = sing[1:]
    su2 = max(exc, key=lambda x: x[1])[0]                     # sigma_u^2 (doubly exc)
    ryd = next((e for e, w in exc if w < 0.3 and e != su2), np.nan)  # Rydberg
    s2, s3 = min(su2, ryd), max(su2, ryd)                     # adiabatic
    u = fci.FCI(mf); u.wfnsym = 'A1u'; u.nroots = 3
    eu, cu = u.kernel(); eu = np.atleast_1d(eu)
    s1 = min((eu[k] for k in range(len(eu))
              if spin_op.spin_square(cu[k], norb, ne)[0] < 0.5), default=np.nan)
    return np.array([s0, s1, s2, s3])


def adc_states(Rb, e_gs, fci_ref):
    from pyscf import adc
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr', verbose=0)
    mf = scf.RHF(mol).run()
    out = [np.nan] * 4              # ADC ground (MP2) is wrong: not plotted
    try:
        a = adc.ADC(mf); a.method = 'adc(2)'; a.method_type = 'ee'; a.verbose = 0
        a.kernel_gs()
        e, *_ = a.kernel(nroots=10)
        ex = np.array([e_gs + x for x in np.atleast_1d(e)])   # on the FCI ground
        for j in (1, 2, 3):                                   # S3 (sigma_u^2) is a
            d = np.abs(ex - fci_ref[j])                       # double -> ADC misses
            out[j] = ex[int(np.argmin(d))] if d.min() < 0.10 else np.nan
    except Exception:
        pass
    return out


def ptc_states(Rb, fci_ref, ext_max=16):
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
    na = norb; dim = na * na
    h2e = fci.direct_spin1.absorb_h1e(h1, eri, norb, (1, 1), 0.5)
    H = np.empty((dim, dim)); e0 = np.zeros(dim)
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
    from nonsym_tda_eig import nonsym_tda_eig
    Hbar = expm(-T2) @ H @ expm(T2)
    sub = np.ix_(mr, mr)
    ee, vr, vl, info = nonsym_tda_eig(Hbar[sub], len(mr))
    Pd = np.array([np.prod([par[P // 2] for P in basis[i]]) for i in mr])
    su2_det = didx([1], [1])                       # sigma_u^2 determinant
    su2_loc = mr.index(su2_det) if su2_det in mr else None
    S2sub = S2[sub]
    g_singlets = []
    for k in np.argsort(ee.real):
        v = vr[:, k]; nv = v @ v
        if float(v @ S2sub @ v / nv) > 0.5:
            continue
        if float(v @ (Pd * v) / nv) > 0:
            w_su2 = (v[su2_loc] ** 2 / nv) if su2_loc is not None else 0.0
            g_singlets.append((ee[k].real, w_su2))
    out = [np.nan] * 4
    if g_singlets:
        out[0] = g_singlets[0][0]
        exc = g_singlets[1:]
        if exc:
            su2 = max(exc, key=lambda x: x[1])[0]               # sigma_u^2 state
            ryd = next((e for e, w in exc if w < 0.3 and e != su2), np.nan)
            out[2], out[3] = min(su2, ryd), max(su2, ryd)       # adiabatic
    # 1Su+ : lowest u singlet
    u_singlets = sorted(ee[k].real for k in range(len(ee))
                        if float(vr[:, k] @ S2sub @ vr[:, k] /
                                 (vr[:, k] @ vr[:, k])) < 0.5 and
                        float(vr[:, k] @ (Pd * vr[:, k]) /
                              (vr[:, k] @ vr[:, k])) < 0)
    out[1] = u_singlets[0] if u_singlets else np.nan
    return np.array(out)


def main():
    Rs = np.unique(np.concatenate([np.linspace(0.7, 8.0, 20),
                                   np.linspace(2.0, 3.2, 13)]))
    EF = np.array([fci_states(R) for R in Rs])
    shift = EF[-1, 0]                                 # S0 dissociation limit
    PT = np.array([ptc_states(R, EF[i]) for i, R in enumerate(Rs)])
    AD = np.array([adc_states(R, EF[i, 0], EF[i]) for i, R in enumerate(Rs)])

    plt.figure(figsize=(8, 6))
    for k in range(4):
        c = COLORS[k]
        plt.plot(Rs, EF[:, k] - shift, '-', color=c, lw=2.4, label=LABELS[k])
        plt.plot(Rs, PT[:, k] - shift, ':', color=c, lw=2.4)
        plt.plot(Rs, AD[:, k] - shift, '--', color=c, lw=1.8)
    # legend: states (color) + methods (style)
    from matplotlib.lines import Line2D
    style = [Line2D([0], [0], color='k', ls='-', label='FCI'),
             Line2D([0], [0], color='k', ls=':', label='pTC-MRSF-CIS'),
             Line2D([0], [0], color='k', ls='--', label='ADC(2)')]
    leg1 = plt.legend(loc='upper right', framealpha=0.95, fontsize=9)
    plt.gca().add_artist(leg1)
    plt.legend(handles=style, loc='lower right', framealpha=0.95, fontsize=9)
    plt.xlabel('R (Bohr)')
    plt.ylabel('Energy relative to 2 H limit (Hartree)')
    plt.title('H$_2$ singlet states / cc-pVTZ')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('pes_h2.png', dpi=130)
    print("wrote pes_h2.png")


if __name__ == "__main__":
    main()
