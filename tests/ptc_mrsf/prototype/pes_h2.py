"""
H2 dissociation, singlet states, cc-pVTZ. Energies relative to the S0 (X 1Sg+)
dissociation limit (-> 0). Two methods, same color per state:
    FCI  : solid     pTC-MRSF-CIS : dotted

States (all singlets):
    S0 = 1 ^1Sigma_g^+ (X)            S1 = 1 ^1Sigma_u^+ (B)
    S2 = 2 ^1Sigma_g^+ (Rydberg E/EF) S3 = 3 ^1Sigma_g^+ (doubly-exc sigma_u^2)
S3 is HIGH at short R (sigma_u^2 is a double excitation) and descends; FCI and
pTC-MRSF-CIS both capture this multireference state.

FCI is standard pyscf; pTC-MRSF-CIS is this repo's reference code.

Run:  python3 pes_h2.py  ->  pes_h2.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from pyscf import gto, scf, fci, ao2mo, symm

from tc_finite_basis import ann, cre, ms0_determinants
from ptc_mrsf_cis import build_s2

BASIS = "cc-pvtz"
ONLY_FCI = False           # FCI(solid) + pTC-MRSF-CIS(dotted) + ADC(2)(dashed)
COLORS = ['black', 'red', 'blue', 'green']     # match the reference figure
LABELS = [r'$1\,^1\Sigma_g^+$', r'$1\,^1\Sigma_u^+$',
          r'$2\,^1\Sigma_g^+$', r'$3\,^1\Sigma_g^+$']


def fci_states(Rb):
    """[S0, S1, S2, S3] singlets, by ADIABATIC energy order within each spatial
    symmetry (eigenvalues are continuous in R, so the curves are smooth):
      S0 = 1 1Sg+ (X), S2 = 2 1Sg+, S3 = 3 1Sg+   (1st/2nd/3rd 1Sg+ singlets)
      S1 = 1 1Su+ (B)                              (lowest 1Su+ singlet).
    A SINGLET-restricted, symmetry-constrained solver (direct_spin0_symm with
    wfnsym) is essential: an unconstrained FCI mixes the near-degenerate A1g
    singlet/triplet roots (giving fractional <S^2> and discontinuous curves)."""
    mol = gto.M(atom=f'H 0 0 0; H 0 0 {Rb}', basis=BASIS, unit='Bohr',
                symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    norb = mf.mo_coeff.shape[1]
    ne = mol.nelectron
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C)
    cis = fci.direct_spin0_symm.FCI(mol)

    def solve(sym, nr):
        cis.wfnsym = sym; cis.nroots = nr
        e, _ = cis.kernel(h1, eri, norb, ne, ecore=mol.energy_nuc(),
                          orbsym=mf.orbsym)
        return np.sort(np.atleast_1d(e))
    g = solve('A1g', 3)        # 1,2,3 1Sg+  (singlet-pure, gerade)
    u = solve('A1u', 1)        # 1 1Su+ (B)  (singlet-pure, ungerade)
    return np.array([g[0], u[0], g[1], g[2]])


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
        # S3 (3 1Sg+) is sigma_u^2, a genuine DOUBLE excitation: it lives in the
        # 2p-2h block of ADC(2), described there only at zeroth order (no
        # correlation), so ADC(2) cannot represent it -> not plotted.
        for j in (1, 2):
            d = np.abs(ex - fci_ref[j])
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
    S2sub = S2[sub]
    # non-Hermitian Hbar: spin/parity must be the BIORTHONORMAL expectation
    # <vl|O|vr> (vl.T @ vr = I from nonsym_tda_eig), NOT <vr|O|vr>; otherwise the
    # 1Su+ triplet (-> 0) is misread as a singlet and stolen into the S1 curve.
    den = np.einsum('ik,ik->k', vl, vr)
    parity = (np.einsum('ik,ik->k', vl, Pd[:, None] * vr) / den).real
    # NOTE: the non-Hermitian transcorrelated H_bar spoils the <S^2> labels
    # (fractional/inverted values), so spin is NOT used here. Parity (gerade /
    # ungerade) IS clean. We track each physical state by parity + nearest FCI
    # energy -- the same scheme used for ADC. pTC keeps its OWN energies; this
    # only assigns which root is which state (and the lowest u root, the b 3Su+
    # triplet -> 0, is correctly NOT picked for the S1 = B 1Su+ singlet curve).
    g_states = sorted(ee[k].real for k in range(len(ee)) if parity[k] > 0)
    u_states = sorted(ee[k].real for k in range(len(ee)) if parity[k] < 0)

    def nearest(cands, ref, used, tol=0.2):
        best, bd = np.nan, tol
        for e in cands:
            if any(abs(e - x) < 1e-9 for x in used):
                continue
            if np.isfinite(ref) and abs(e - ref) < bd:
                bd, best = abs(e - ref), e
        return best

    out = [np.nan] * 4
    if g_states:
        out[0] = g_states[0]                              # S0 = lowest 1Sg+
    out[1] = nearest(u_states, fci_ref[1], [])            # S1 = B 1Su+ (singlet)
    out[2] = nearest(g_states[1:], fci_ref[2], [out[0]])  # S2 = 2 1Sg+
    out[3] = nearest(g_states[1:], fci_ref[3], [out[0], out[2]])  # S3 = 3 1Sg+
    return np.array(out)


CKPT = "pes_h2_data.npz"


def compute(Rs):
    """Compute (or resume) FCI/pTC/ADC for all R; checkpoint after each point.
    Re-run any time -- finished points are loaded from CKPT and skipped."""
    n = len(Rs)
    EF = np.full((n, 4), np.nan)
    PT = np.full((n, 4), np.nan)
    AD = np.full((n, 4), np.nan)
    if os.path.exists(CKPT):
        d = np.load(CKPT)
        oR = d["Rs"]                            # merge prior points by R value
        for i, R in enumerate(Rs):
            j = np.where(np.abs(oR - R) < 1e-9)[0]
            if len(j) and not np.isnan(d["EF"][j[0], 0]):
                EF[i], PT[i], AD[i] = d["EF"][j[0]], d["PT"][j[0]], d["AD"][j[0]]
        print(f"resumed from {CKPT}: "
              f"{int(np.sum(~np.isnan(EF[:, 0])))}/{n} points reused")
    for i, R in enumerate(Rs):
        if not np.isnan(EF[i, 0]):              # already computed
            continue
        EF[i] = fci_states(R)
        PT[i] = ptc_states(R, EF[i])
        AD[i] = adc_states(R, EF[i, 0], EF[i])
        np.savez(CKPT, Rs=Rs, EF=EF, PT=PT, AD=AD)   # checkpoint each point
        print(f"  done R={R:.3f} Bohr ({i + 1}/{n})", flush=True)
    return EF, PT, AD


def main():
    Rs = np.unique(np.concatenate([np.linspace(0.7, 8.0, 20),
                                   np.linspace(0.9, 2.0, 12),   # steep short-R region
                                   np.linspace(2.0, 3.2, 13)])) # avoided crossing
    EF, PT, AD = compute(Rs)
    shift = EF[-1, 0]                                 # S0 dissociation limit

    plt.figure(figsize=(8, 6))
    for k in range(4):
        c = COLORS[k]
        plt.plot(Rs, EF[:, k] - shift, '-', color=c, lw=2.4, label=LABELS[k])
        if not ONLY_FCI:
            plt.plot(Rs, PT[:, k] - shift, ':', color=c, lw=2.4)
    # legend: states (color) + methods (style)
    from matplotlib.lines import Line2D
    style = [Line2D([0], [0], color='k', ls='-', label='FCI'),
             Line2D([0], [0], color='k', ls=':', label='pTC-MRSF-CIS')]
    leg1 = plt.legend(loc='upper right', framealpha=0.95, fontsize=9)
    plt.gca().add_artist(leg1)
    if not ONLY_FCI:
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
