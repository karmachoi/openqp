#!/usr/bin/env python3
"""
QMRSF Ŝ² spin-projection of the CAS(4,4) determinant-union backbone (pure NumPy, NO pyscf/scipy).

PROBLEM (production blocker):
  The QMRSF determinant-union backbone is built in the Ms=0 sector of CAS(4,4):
  4 active spatial orbitals, 4 electrons, 2 alpha + 2 beta -> C(4,2)*C(4,2) = 6*6 = 36 determinants.
  This Ms=0 space is spin-CONTAMINATED: it carries S=0 (singlet), S=1 (triplet) and S=2 (quintet)
  components all mixed together. The QMRSF target -- the lowest SINGLET -- is buried under the
  triplet/quintet contaminants. We must project the 36-det space into spin-pure blocks.

WHAT THIS PROTOTYPE DOES:
  1. Build the CAS(4,4) Ms=0 determinant space (36 dets) over a 4-orbital PPP active window.
  2. Build the Ŝ² operator MATRIX in that determinant basis from second-quantized
     Ŝ+ = sum_p a+_{p,alpha} a_{p,beta},  Ŝ- = sum_p a+_{p,beta} a_{p,alpha},  Ŝz,
     using   Ŝ² = Ŝ- Ŝ+ + Ŝz(Ŝz+1).
     Verify Ŝ² is Hermitian and [H, Ŝ²] = 0.
  3. Diagonalize Ŝ² -> the Ms=0 space splits into  S=0 (eig 0, dim 20),
     S=1 (eig 2, dim 15), S=2 (eig 6, dim 1).
  4. Project H into each spin block (P_S^T H P_S), diagonalize, and verify every
     block energy is a member of the full-FCI(=full-CAS) spectrum with the correct <S^2>.
     Show the lowest singlet (the QMRSF target) is cleanly extracted.
  5. Document the MATRIX-FREE recipe for the production Davidson sigma builder.

Counting check (Weyl/branching for 4 e- in 4 orbitals, Ms=0):
    full Ms=0 dim                = 36
    S=0 multiplicity count       = 20   (number of CSFs * (2S+1) folded into Ms=0 ... see note)
    S=1                          = 15
    S=2                          =  1
  The Ms=0 determinant count partitions by total spin as
    dim(Ms=0) = sum_S  N_CSF(S),   where in the Ms=0 row each spin S contributes
    N_CSF(S) determinant-space dimension = (# independent Ms=0 vectors of that S).
  For (4e,4o): N_CSF(S=0)=20, N_CSF(S=1)=15, N_CSF(S=2)=1  -> 20+15+1 = 36.  [verified below]

Run:  python3 qmrsf_s2_projection_proto.py
"""
import numpy as np

# Reuse the validated machinery from the icPT2 prototype (same directory).
from qmrsf_icpt2_ppp_proto import build_ppp, spinorb, gen_dets, melem, build_H


# ----------------------------------------------------------------------
# Spin operators in the determinant basis
# ----------------------------------------------------------------------
# Spin-orbital index convention (matches spinorb / gen_dets):
#   P in [0, n)     -> alpha electron in spatial orbital P
#   P in [n, 2n)    -> beta  electron in spatial orbital (P-n)
# A determinant is a sorted tuple of occupied spin-orbital indices.

def _apply_excitation(det, p_create, q_annih):
    """Apply a+_{p_create} a_{q_annih} to determinant `det` (sorted tuple of SO indices).
    Returns (new_det_sorted, sign) or (None, 0) if it annihilates the state.
    Jordan-Wigner / Slater-Condon sign bookkeeping consistent with melem()."""
    occ = list(det)
    if q_annih not in occ:                      # nothing to annihilate
        return None, 0
    sign = 1
    idx = occ.index(q_annih)                    # sign from a_{q}
    sign *= (-1) ** idx
    occ.pop(idx)
    if p_create in occ:                         # Pauli: orbital already occupied
        return None, 0
    idx = sum(1 for o in occ if o < p_create)   # sign from a+_{p}
    sign *= (-1) ** idx
    occ.insert(idx, p_create)
    return tuple(occ), sign


def build_Sz(dets, n):
    """Ŝz = 1/2 (n_alpha - n_beta), diagonal in the determinant basis."""
    N = len(dets)
    Sz = np.zeros((N, N))
    for j, d in enumerate(dets):
        na = sum(1 for so in d if so < n)
        nb = sum(1 for so in d if so >= n)
        Sz[j, j] = 0.5 * (na - nb)
    return Sz


def build_Sm_Sp(dets, n):
    """Build the matrix of the PRODUCT  Ŝ- Ŝ+  directly in the determinant basis.

    Ŝ- Ŝ+ = ( sum_q a+_{q,beta} a_{q,alpha} ) ( sum_p a+_{p,alpha} a_{p,beta} )
           = sum_{p,q} a+_{q,beta} a_{q,alpha} a+_{p,alpha} a_{p,beta}.

    CRUCIAL: this is a genuine TWO-body operator that maps Ms=0 -> Ms=0, but the
    INTERMEDIATE after Ŝ+ alone lives in the Ms=+1 sector. We therefore apply Ŝ+
    then Ŝ- as a composite, looking up only the FINAL determinant in the Ms=0
    basis. (Building Ŝ+ as a standalone matrix in the Ms=0 basis gives zero,
    because Ŝ+ leaves the Ms=0 space -- that was the trap.)
    Returns (N,N) ndarray O[i,j] = <det_i | Ŝ- Ŝ+ | det_j>."""
    N = len(dets)
    index = {d: i for i, d in enumerate(dets)}
    O = np.zeros((N, N))
    for j, d in enumerate(dets):
        for p in range(n):                       # Ŝ+ term: a+_{alpha p} a_{beta p}
            mid, s1 = _apply_excitation(d, p, p + n)
            if mid is None:
                continue
            for q in range(n):                   # Ŝ- term: a+_{beta q} a_{alpha q}
                fin, s2 = _apply_excitation(mid, q + n, q)
                if fin is None:
                    continue
                i = index.get(fin)               # final det IS back in Ms=0
                if i is None:
                    continue
                O[i, j] += s1 * s2
    return O


def build_Sp_Sm(dets, n):
    """Build the PRODUCT  Ŝ+ Ŝ-  directly (mirror of build_Sm_Sp), for the
    independent cross-check  Ŝ² = Ŝ+ Ŝ- + Ŝz(Ŝz-1).
    Ŝ+ Ŝ- = sum_{p,q} a+_{p,alpha} a_{p,beta} a+_{q,beta} a_{q,alpha}.
    Apply Ŝ- (to the Ms=-1 sector) then Ŝ+ (back to Ms=0)."""
    N = len(dets)
    index = {d: i for i, d in enumerate(dets)}
    O = np.zeros((N, N))
    for j, d in enumerate(dets):
        for q in range(n):                       # Ŝ- term: a+_{beta q} a_{alpha q}
            mid, s1 = _apply_excitation(d, q + n, q)
            if mid is None:
                continue
            for p in range(n):                   # Ŝ+ term: a+_{alpha p} a_{beta p}
                fin, s2 = _apply_excitation(mid, p, p + n)
                if fin is None:
                    continue
                i = index.get(fin)
                if i is None:
                    continue
                O[i, j] += s1 * s2
    return O


def build_S2(dets, n):
    """Ŝ² = Ŝ- Ŝ+ + Ŝz (Ŝz + 1)  (acts within the fixed Ms=0 determinant basis)."""
    SmSp = build_Sm_Sp(dets, n)
    Sz = build_Sz(dets, n)
    I = np.eye(len(dets))
    S2 = SmSp + Sz @ (Sz + I)
    return S2, SmSp, Sz


# ----------------------------------------------------------------------
# Spin-block projectors from the Ŝ² eigendecomposition
# ----------------------------------------------------------------------
def spin_blocks(S2, tol=1e-6):
    """Diagonalize Ŝ² and group eigenvectors by eigenvalue s(s+1).
    Returns list of (s2_value, S_quantum, projector_columns P_S (N x dim)) sorted by s2."""
    w, V = np.linalg.eigh(S2)
    # cluster eigenvalues
    blocks = []
    order = np.argsort(w)
    w = w[order]; V = V[:, order]
    i = 0
    N = len(w)
    groups = []
    while i < N:
        j = i + 1
        while j < N and abs(w[j] - w[i]) < tol:
            j += 1
        groups.append((np.mean(w[i:j]), V[:, i:j].copy()))
        i = j
    for s2val, P in groups:
        # s(s+1) = s2val -> s = (-1 + sqrt(1+4 s2val))/2
        s = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * max(s2val, 0.0)))
        blocks.append((s2val, s, P))
    return blocks


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------
def main():
    np.set_printoptions(precision=5, suppress=True)
    print("=" * 80)
    print("QMRSF  S^2 spin-projection of the CAS(4,4) Ms=0 determinant-union backbone")
    print("=" * 80)

    # ---- 1. Build a 4-orbital PPP active model -> CAS(4,4) Ms=0 space ----
    # A 4-center pi system IS a 4-orbital model: 4 active orbitals, 4 electrons, Ms=0.
    n = 4
    h_mo, eri_mo, eps = build_ppp(n)            # 4 orbitals, full space = the active space
    H1, g, _ = spinorb(h_mo, eri_mo)
    na = nb = 2                                  # 2 alpha + 2 beta  -> Ms = 0
    dets = gen_dets(n, na, nb)                   # all Ms=0 dets in the 4-orbital space
    N = len(dets)
    print(f"\n[1] CAS(4,4) Ms=0 space: n_active={n}, na={na}, nb={nb}, "
          f"determinants N = C(4,2)^2 = {N}")
    assert N == 36, f"expected 36 determinants, got {N}"

    # Build the determinant Hamiltonian (full CAS = full FCI in this active space)
    H = build_H(dets, H1, g)
    print(f"    H built: shape {H.shape},  Hermiticity max|H-H^T| = "
          f"{np.abs(H - H.T).max():.2e}")

    # ---- 2. Build S^2, verify Hermitian and [H,S^2]=0 ----
    S2, SmSp, Sz = build_S2(dets, n)
    herm = np.abs(S2 - S2.T).max()
    comm = H @ S2 - S2 @ H
    comm_res = np.abs(comm).max()
    print(f"\n[2] S^2 operator (S^2 = S- S+ + Sz(Sz+1)):")
    print(f"    S^2 Hermiticity   max|S^2 - S^2^T| = {herm:.2e}")
    print(f"    commutator        max|[H, S^2]|    = {comm_res:.2e}   ([H,S^2]=0 expected)")
    # cross-check the alternative identity S^2 = S+ S- + Sz(Sz-1) (built independently)
    SpSm = build_Sp_Sm(dets, n)
    I = np.eye(N)
    S2_alt = SpSm + Sz @ (Sz - I)
    print(f"    identity check    max|[S- S+ +Sz(Sz+1)] - [S+ S- +Sz(Sz-1)]| = "
          f"{np.abs(S2 - S2_alt).max():.2e}")

    # ---- 3. Diagonalize S^2, group by s(s+1) -> block dimensions ----
    blocks = spin_blocks(S2)
    print(f"\n[3] S^2 eigenvalue spectrum -> spin blocks of the Ms=0 space:")
    print(f"    {'s(s+1)':>8} {'S':>4} {'2S+1':>5} {'dim':>5}   label")
    label_map = {0: 'singlet', 1: 'triplet', 2: 'quintet'}
    total = 0
    expected = {0.0: 20, 2.0: 15, 6.0: 1}
    dim_by_s2 = {}
    for s2val, s, P in blocks:
        dim = P.shape[1]
        total += dim
        Sround = int(round(s))
        mult = 2 * Sround + 1
        lbl = label_map.get(Sround, f"S={Sround}")
        dim_by_s2[round(s2val, 3)] = dim
        print(f"    {s2val:>8.4f} {s:>4.1f} {mult:>5d} {dim:>5d}   {lbl}")
    print(f"    {'':>8} {'':>4} {'':>5} {total:>5d}   (sum)")
    # verify the 20/15/1 split
    ok_dims = (dim_by_s2.get(0.0) == 20 and dim_by_s2.get(2.0) == 15
               and dim_by_s2.get(6.0) == 1 and total == 36)
    print(f"    block-dimension gate: S=0 -> {dim_by_s2.get(0.0)} (want 20), "
          f"S=1 -> {dim_by_s2.get(2.0)} (want 15), S=2 -> {dim_by_s2.get(6.0)} (want 1)  "
          f"=> {'PASS' if ok_dims else 'FAIL'}")

    # ---- 4. Project H into each spin block, diagonalize, verify ----
    print(f"\n[4] Project H into each spin block:  H_S = P_S^T H P_S,  diagonalize")
    e_fci_full = np.sort(np.linalg.eigvalsh(H))             # full CAS/FCI spectrum (all spins)
    print(f"    full CAS/FCI spectrum ({N} states), lowest 8: "
          f"{np.round(e_fci_full[:8], 5)}")
    singlet_energies = None
    for s2val, s, P in blocks:
        Sround = int(round(s))
        lbl = label_map.get(Sround, f"S={Sround}")
        HS = P.T @ H @ P                                    # block Hamiltonian
        eS, cS = np.linalg.eigh(HS)
        # back-transform eigenvectors to full det basis to measure <S^2>
        full_vecs = P @ cS
        s2_exp = np.array([v @ (S2 @ v) for v in full_vecs.T])
        # confirm each block energy is in the full FCI spectrum
        in_fci = np.array([np.min(np.abs(e_fci_full - e)) for e in eS])
        max_match = in_fci.max()
        print(f"\n    --- {lbl}  (S={s:.1f}, s(s+1)={s2val:.3f}, dim={P.shape[1]}) ---")
        print(f"        energies (lowest {min(6, len(eS))}): {np.round(eS[:6], 5)}")
        print(f"        <S^2> per state (lowest 6): {np.round(s2_exp[:6], 5)}  "
              f"(should all = {s2val:.3f})")
        print(f"        max <S^2> deviation = {np.abs(s2_exp - s2val).max():.2e}")
        print(f"        max |E_block - nearest E_FCI| = {max_match:.2e}  "
              f"(0 => block energies are exact FCI states)")
        if Sround == 0:
            singlet_energies = eS

    # ---- the QMRSF target: lowest singlet, cleanly extracted ----
    e_singlet_lowest = singlet_energies.min()
    # what is the global ground state, and is it singlet or contaminant?
    e_global = e_fci_full[0]
    print(f"\n[5] QMRSF TARGET EXTRACTION:")
    print(f"    lowest SINGLET energy  (S=0 block min)  = {e_singlet_lowest:.6f}")
    print(f"    global CAS/FCI ground state             = {e_global:.6f}")
    # measure <S^2> of the raw (unprojected) global ground state
    w_all, V_all = np.linalg.eigh(H)
    gs = V_all[:, 0]
    gs_s2 = gs @ (S2 @ gs)
    gs_S = 0.5 * (-1 + np.sqrt(1 + 4 * max(gs_s2, 0)))
    print(f"    <S^2> of raw global ground state        = {gs_s2:.5f}  (S={gs_S:.2f})")
    # show the spin labels of the lowest raw states -> singlets ARE buried/interleaved
    s2_raw = np.array([V_all[:, k] @ (S2 @ V_all[:, k]) for k in range(min(8, N))])
    labs = []
    for v in s2_raw:
        Sk = int(round(0.5 * (-1 + np.sqrt(1 + 4 * max(v, 0)))))
        labs.append(label_map.get(Sk, f"S{Sk}")[0].upper())
    sing_ranks = [k for k in range(N) if abs(V_all[:, k] @ (S2 @ V_all[:, k])) < 1e-6][:3]
    print(f"    raw spectrum spin labels (lowest 8)     = {labs}")
    print(f"    -> singlets sit at global ranks {sing_ranks}: the 1st EXCITED singlet is")
    print(f"       BURIED beneath triplet states -- exactly the contamination the QMRSF")
    print(f"       backbone faces. The spin-pure S=0 block extracts ALL singlets cleanly,")
    print(f"       free of triplet/quintet contamination.")

    # ---- consistency: union of all spin-block spectra == full FCI spectrum ----
    all_block_e = []
    for s2val, s, P in blocks:
        HS = P.T @ H @ P
        all_block_e.extend(np.linalg.eigvalsh(HS).tolist())
    all_block_e = np.sort(np.array(all_block_e))
    union_match = np.abs(all_block_e - e_fci_full).max()
    print(f"\n[6] CONSISTENCY: sorted(union of all spin-block energies) vs full FCI spectrum")
    print(f"    max|E_union - E_FCI| = {union_match:.2e}  "
          f"(0 => projection is exact & complete, no states lost/duplicated)")

    # ---- final verdict ----
    print("\n" + "=" * 80)
    all_pass = (np.abs(H - H.T).max() < 1e-10 and herm < 1e-10 and comm_res < 1e-9
                and ok_dims and union_match < 1e-9)
    print(f"VERDICT: {'ALL GATES PASS' if all_pass else 'CHECK FAILURES ABOVE'}")
    print("=" * 80)

    # ---- 5/doc. Matrix-free recipe for the production Davidson ----
    print(__doc_matrix_free__)


__doc_matrix_free__ = r"""
------------------------------------------------------------------------------
MATRIX-FREE S^2 PROJECTION FOR THE PRODUCTION DAVIDSON (sigma-builder backbone)
------------------------------------------------------------------------------
The production QMRSF backbone never forms H explicitly; it provides a sigma
builder  sigma = H @ x  (and, analogously, can provide  s2x = S^2 @ x).
S^2 is a cheap ONE-BODY-like operator on determinants (it is the sum of the
single intra-spatial spin flips Ŝ+ Ŝ- + Ŝz^2 ...), so a matrix-free
"apply_S2(x)" routine is built exactly like apply_H but with the spin operators
above -- O(N_det) per spatial orbital, far cheaper than the sigma build.

Two production-ready strategies (both pure matrix-free, no explicit H or P_S):

  (A) SPIN-SHIFT (penalty) Davidson  -- recommended, simplest:
      Solve the shifted eigenproblem with sigma builder
          sigma_shift(x) = H@x + mu * (S^2 - s_t(s_t+1) I) @ x
      where s_t is the TARGET spin (0 for the QMRSF singlet) and mu is a large
      positive shift (mu ~ 10-100 * spectral width of H). Any vector with the
      wrong S^2 is pushed up by mu*(s(s+1)-s_t(s_t+1)) > 0, so the desired
      spin-pure roots fall to the bottom and Davidson converges to them. Needs
      ONLY apply_H and apply_S2 -- no projector, no extra storage. This is the
      "(S^2 - s(s+1)) penalty/shift" approach.

  (B) PROJECTED Davidson -- spin-adapt the Krylov subspace each iteration:
      Apply the spin projector to every trial/expansion vector before the sigma
      build, using the matrix-free spectral projector built from apply_S2:
          P_{s_t} = prod_{s' != s_t} (S^2 - s'(s'+1) I) / (s_t(s_t+1) - s'(s'+1))
      i.e.  x <- P_{s_t} x  via repeated apply_S2 (a degree-(#contaminants)
      polynomial in S^2; for CAS(4,4) only S in {0,1,2}, so a quadratic in S^2).
      Concretely, to kill triplet(2) and quintet(6) and keep singlet(0):
          x_pure = (S^2 - 2I)(S^2 - 6I) x / [(0-2)(0-6)]
      Because [H,S^2]=0, H maps the spin-pure subspace into itself, so the whole
      Davidson run stays in the singlet block; <S^2> of every Ritz vector is
      pinned at the target. Re-project after preconditioning to remove
      round-off leakage.

  Practical notes:
    * Build apply_S2 once from Ŝ+, Ŝ-, Ŝz on determinants (the build_spin_op
      logic above, lifted into the sigma-builder's determinant addressing).
    * Use (A) for robustness/speed; fall back to (B) if near-degenerate spin
      states of DIFFERENT multiplicity sit within mu-resolution.
    * Monitor convergence by <x|S^2|x> -> s_t(s_t+1); it is the spin-purity gate.
    * Cost: apply_S2 is ~ n_active single-flip passes over the CI vector, i.e.
      cheaper than one sigma build, so spin projection adds negligible overhead.
------------------------------------------------------------------------------
"""


if __name__ == "__main__":
    main()
