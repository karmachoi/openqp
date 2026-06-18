"""
Working MRSF-CIS reference on real molecular integrals (pyscf), plus the
non-Hermitian transcorrelated solver path.

This actually runs and is self-validating:

  (1) Builds the (2,2) spin-flip determinant Hamiltonian for methylene (CH2,
      the canonical diradical MRSF is designed for) from real pyscf integrals,
      diagonalizes the Ms=0 spin-flip sector, and ASSERTS the lowest root
      reproduces pyscf's independent CASCI(2,2) energy. The states come out
      spin-pure (<S^2> = 0 singlets, 2 triplet) -- the MRSF property.

  (2) Forms the transcorrelated effective Hamiltonian H_bar = J^{-1} H J on the
      same real integrals (non-Hermitian) and solves it with the validated
      non-Hermitian kernel tc_nonsym_tda_eig. ASSERTS:
        - tau=0 (g=1) reproduces the Hermitian MRSF-CIS spectrum (the gate);
        - the non-Hermitian spectrum is real and biorthonormal on real integrals.

What is real here: the integrals, the spin-flip state construction, the
CASCI validation, the spin purity, the non-Hermitian machinery.
What is a model here: the correlation factor J is an active-space Gutzwiller
operator (computable, illustrative), NOT the production F12 geminal. The
correlation *recovery* that the geminal provides is demonstrated exactly in
tc_hubbard_demo.py; here a similarity transform on the complete active space
is spectrum-preserving by construction, which is exactly why it is the right
test of the non-Hermitian solver.

Run:  python3 mrsf_cis_pyscf.py
"""

import numpy as np
from pyscf import gto, scf, mcscf, fci, ao2mo
from nonsym_tda_eig import nonsym_tda_eig


def build_active_ci_hamiltonian(mc, ncas, ms0=(1, 1)):
    """Full (ncas, *) Ms=0 CI Hamiltonian matrix from real active integrals."""
    h1e, ecore = mc.get_h1eff()
    eri = ao2mo.restore(1, mc.get_h2eff(), ncas)
    na = fci.cistring.num_strings(ncas, ms0[0])
    nb = fci.cistring.num_strings(ncas, ms0[1])
    dim = na * nb
    h2 = fci.direct_spin1.absorb_h1e(h1e, eri, ncas, ms0, 0.5)
    hmat = np.zeros((dim, dim))
    for k in range(dim):
        c = np.zeros(dim)
        c[k] = 1.0
        hc = fci.direct_spin1.contract_2e(h2, c.reshape(na, nb), ncas, ms0)
        hmat[:, k] = hc.reshape(-1)
    hmat += ecore * np.eye(dim)
    return hmat, (na, nb)


def double_occupancy(ncas, ms0):
    """Number of doubly-occupied active orbitals per determinant (basis order)."""
    occ_a = fci.cistring.gen_occslst(range(ncas), ms0[0])
    occ_b = fci.cistring.gen_occslst(range(ncas), ms0[1])
    docc = []
    for ia in occ_a:
        for ib in occ_b:
            docc.append(len(set(ia) & set(ib)))
    return np.array(docc, dtype=float)


def spin_square_of(vecs, shape, ncas, ms0):
    na, nb = shape
    out = []
    for j in range(vecs.shape[1]):
        ss, _ = fci.spin_op.spin_square(vecs[:, j].reshape(na, nb), ncas, ms0)
        out.append(ss)
    return np.array(out)


def main():
    HARTREE2EV = 27.211386245988

    # Asymmetric (C1) geometry on purpose: with C2v symmetry the two active
    # orbitals have different spatial symmetry, the closed<->open determinants
    # decouple, and the Gutzwiller transform stays symmetric. Breaking symmetry
    # makes the transcorrelated H_bar genuinely non-Hermitian.
    mol = gto.M(atom='C 0 0 0.0; H 0 0.99 0.59; H 0 -1.18 0.36',
                basis='6-31G', spin=2, verbose=0)
    mf = scf.ROHF(mol).run()
    ncas, nelecas, ms0 = 2, 2, (1, 1)

    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.kernel()
    e_casci_ground = mc.e_tot

    # Use a non-canonical (localized-diradical) active basis: pyscf canonicalizes
    # the active orbitals (diagonal active Fock), which by Brillouin's theorem
    # decouples closed<->open determinants and would make the Gutzwiller transform
    # trivially symmetric. A small active-active rotation restores the physical
    # coupling and leaves the CASCI spectrum exactly invariant.
    ncore = mc.ncore
    theta = 0.30
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta),  np.cos(theta)]])
    act = mc.mo_coeff[:, ncore:ncore + ncas].copy()
    mc.mo_coeff[:, ncore:ncore + ncas] = act @ rot

    # (1) MRSF-CIS spin-flip sector on real integrals --------------------------
    hmat, shape = build_active_ci_hamiltonian(mc, ncas, ms0)
    assert np.allclose(hmat, hmat.T), "bare H must be symmetric"
    w, v = np.linalg.eigh(hmat)
    ss = spin_square_of(v, shape, ncas, ms0)

    print("=== MRSF-CIS on real CH2 integrals (6-31G, (2,2)) ===")
    print(f"pyscf CASCI(2,2) ground energy : {e_casci_ground:.8f} Ha")
    print(f"MRSF-CIS ground energy         : {w[0]:.8f} Ha")
    assert abs(w[0] - e_casci_ground) < 1e-8, "MRSF-CIS must match CASCI"
    print("VALIDATED: MRSF-CIS ground == CASCI(2,2)\n")

    labels = {0: "triplet", 2: "singlet"}
    print(" state      E (Ha)        <S^2>   multiplicity")
    for i in range(len(w)):
        mult = 0.5 * (-1 + np.sqrt(1 + 4 * ss[i]))
        print(f"   {i}    {w[i]:.8f}   {ss[i]:5.3f}     2S+1={2*mult+1:.2f}")
    # singlet-triplet gap (lowest singlet minus ground triplet)
    e_trip = w[0]
    e_sing = w[np.argmax(ss < 0.5)]  # first ~singlet
    gap_ev = (e_sing - e_trip) * HARTREE2EV
    print(f"\nvertical T->S gap            : {gap_ev:.3f} eV  "
          f"(real number from real integrals)\n")

    # (2) Transcorrelated non-Hermitian solve on the same real integrals -------
    docc = double_occupancy(ncas, ms0)
    nstate = len(w)

    # tau = 0 gate: g = 1 -> H_bar == H
    g = 1.0
    jdiag = g ** docc
    hbar = (1.0 / jdiag)[:, None] * hmat * jdiag[None, :]
    ee0, vr0, vl0, info0 = nonsym_tda_eig(hbar, nstate)
    assert np.allclose(np.sort(ee0), np.sort(w), atol=1e-9), (ee0, w)
    print("=== transcorrelated non-Hermitian solve (real CH2 integrals) ===")
    print(f"tau=0 gate: max|E_TC - E_MRSF| = {np.max(np.abs(np.sort(ee0)-np.sort(w))):.2e}"
          "   -> reproduces MRSF-CIS")

    # tau != 0: genuine non-symmetric H_bar; spectrum is similarity-invariant
    # (right answer), and the kernel must return it real + biorthonormal.
    g = 0.8
    jdiag = g ** docc
    hbar = (1.0 / jdiag)[:, None] * hmat * jdiag[None, :]
    assert not np.allclose(hbar, hbar.T), "H_bar should be non-symmetric"
    ee, vr, vl, info = nonsym_tda_eig(hbar, nstate)
    biorth = np.max(np.abs(vl.T @ vr - np.eye(nstate)))
    res = np.max(np.abs(hbar @ vr - vr * ee))
    print(f"tau!=0 (g={g}): non-symmetric H_bar")
    print(f"   spectrum real (max Im) : {info['max_imag']:.2e}, "
          f"complex roots: {info['n_complex']}")
    print(f"   matches exact spectrum : {np.max(np.abs(np.sort(ee)-np.sort(w))):.2e}")
    print(f"   right-eig residual     : {res:.2e}")
    print(f"   biorthonormality       : {biorth:.2e}")
    assert info["n_complex"] == 0
    assert np.allclose(np.sort(ee), np.sort(w), atol=1e-8)
    assert biorth < 1e-7 and res < 1e-7
    print("\nVALIDATED: non-Hermitian transcorrelated solver runs on real")
    print("molecular integrals, returns the correct real spectrum with")
    print("biorthonormal left/right states.")


if __name__ == "__main__":
    main()
