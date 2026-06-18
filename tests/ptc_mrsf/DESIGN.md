# pTC-MRSF-CIS: Projective Transcorrelated MRSF-CIS

Design notes and phased implementation plan for adding Ten-no's projective
transcorrelation (pTC) on top of MRSF-CIS in OpenQP.

## 1. Idea

MRSF-CIS gives the *right states* — spin-pure singlets/triplets, multireference
character at CIS cost, S0 as a response root, and the correct conical-intersection
topology — but at the bare HF/CIS level it has **no dynamic correlation** and does
**not** describe the electron–electron cusp, so excitation energies are off and
basis-set convergence is slow.

pTC (S. Ten-no, *J. Chem. Phys.* **159**, 171103 (2023)) injects exactly that
missing physics **into the Hamiltonian** rather than through an XC functional, via
a nonunitary similarity transform

```
H_bar = exp(-tau) H exp(tau)
```

whose effective Hamiltonian is non-Hermitian, terminates formally at four-body
interactions, is spin-contamination-free, and satisfies the singlet **and** triplet
first-order cusp conditions simultaneously.

The result, **pTC-MRSF-CIS**, keeps MRSF's states and adds correlation + fast
basis-set convergence, **without a functional and without the DFT
double-counting problem**.

## 2. Why this substrate

- **MRSF over SF-CIS.** The e–e cusp condition is spin-dependent (antiparallel
  slope 1/2, parallel 1/4). pTC installs a definite-spin cusp; a spin-contaminated
  SF-CIS state (single Ms=+1 reference) satisfies neither cusp cleanly. MRSF's
  mixed-reference spin purity matches what pTC is built to exploit.
- **ROHF MRSF-CIS as the pilot.** Chosen for the first implementation because it
  reuses OpenQP's validated MRSF infrastructure and gradients. Known limitation:
  the ROHF reference carries a coupling-coefficient/canonicalization ambiguity and
  triplet-biased orbitals. A future **SA-CAS(4,4)** substrate would remove that
  ambiguity and is the natural long-term target (a (4,4) CASCI is exactly the small
  FCI solver pTC is designed for), at the cost of MRSF's single-reference
  cheapness. Deferred — see §6.

## 3. Scaling

Final formal scaling is **O(N^5)** — MP2 class — set by the integral
transform / transcorrelated-Hamiltonian build, **not** by the response:

| Step | Order |
|---|---|
| ROHF reference (SCF) | O(N^4) |
| AO→MO transform | O(N^5) |
| pTC H_bar build (3-body normal-ordered + DF) | **O(N^5)** ← rate-limiting |
| MRSF-CIS non-Hermitian Davidson | O(N^4) / iteration |

Guaranteeing N^5 requires never materializing the higher-body terms:
explicit 3-body is O(N^6), explicit 4-body O(N^8). **Normal-ordering** the 3-/4-body
operators against the ROHF reference density (folding them into effective 1-/2-body
integrals) collapses the series back to N^5 — consistent with pTC(2) costing ~1/2 of
MP2-F12.

**Reducing scaling.** The N^5 is the MP2 floor for building a correlated
Hamiltonian; a dense canonical algorithm cannot beat it.
- **DF / RI / Cholesky** — keeps N^5, cuts the prefactor ~10x, and is what makes the
  3-body term tractable. **Locked in from Phase 2.**
- **Locality (PNO/DLPNO)** — the correlation factor is short-ranged, so a local
  treatment can reach ~O(N^3)→O(N) asymptotically. Deferred; large-system optimization.
- **THC / tensor decomposition** — can reach O(N^4). Deferred.

## 4. What changes vs. regular MRSF-CIS

| | regular MRSF-CIS | pTC-MRSF-CIS |
|---|---|---|
| Hamiltonian | bare H | similarity-transformed H_bar |
| Correlation | none (HF/CIS) | dynamic + e–e cusp in H_bar |
| Integrals | 1-/2-electron | + 3-body (normal-ordered, DF), modified 2-body |
| Matrix | real symmetric | non-Hermitian |
| Eigensolver | symmetric Davidson (`rpaeig`) | general Davidson, left+right vectors |
| Basis convergence | slow | near-CBS in small bases |

Unchanged: mixed-reference construction, spin completeness, S0-as-response-root,
conical-intersection topology.

## 5. Phased implementation

**Phase 1 — non-Hermitian reduced eigensolver (this commit).**
- `source/modules/tdhf_mrsf_ptc.F90 :: tc_nonsym_tda_eig` — DGEEV-based general
  reduced-space solve returning real eigenvalues with biorthonormal left/right
  Ritz vectors; flags complex (instability) roots.
- Validated two ways: a NumPy reference
  (`tests/ptc_mrsf/prototype/nonsym_tda_eig.py`) and a **compiled Fortran test**
  of the actual module (`tests/ptc_mrsf/prototype/tc_nonsym_eig_test.F90`, build:
  `gfortran source/precision.F90 source/modules/tdhf_mrsf_ptc.F90
  tests/ptc_mrsf/prototype/tc_nonsym_eig_test.F90 -llapack -lblas`).
  **Acceptance gate (met, both):** with tau=0 the non-symmetric solver reproduces
  the symmetric reference eigenpairs to machine precision (max|dE| ~ 1e-14), and
  the non-symmetric case yields a real spectrum with biorthonormal vectors
  (residuals ~1e-14).
- The Fortran module compiles cleanly against `precision` + LAPACK; not yet wired
  into the MRSF hot path (Phase 4).

Runnable, self-validating prototypes in `tests/ptc_mrsf/prototype/`:
- `nonsym_tda_eig.py` -- the non-Hermitian reduced eigensolver kernel + tests.
- `mrsf_cis_pyscf.py` -- MRSF-CIS on **real** methylene (CH2) integrals via pyscf:
  asserts the (2,2) spin-flip sector reproduces pyscf CASCI(2,2), reports spin-pure
  states (<S^2>) and the vertical T->S gap, and runs the non-Hermitian
  transcorrelated solve on the real integrals (tau=0 gate + real biorthonormal
  spectrum). The correlation factor used here is an active-space Gutzwiller model,
  not the production F12 geminal.
- `tc_hubbard_demo.py` -- exact Hubbard-dimer demonstration of *what
  transcorrelation buys*: a single mean-field determinant recovers 100% of the
  correlation energy through J, and is shown to be the exact right-eigenvector of
  the non-Hermitian H_bar (validated against the closed-form ground state).
- `tc_finite_basis.py` -- **Phase-2 mechanism on real molecular integrals**
  (H4/STO-3G). A genuine transcorrelation (the MP2 cluster operator T2 built from
  the integrals, not a model) forms H_bar = e^{-T2} H e^{T2}, which is downfolded
  into the compact HF+singles (CIS) space. Validated non-circularly: the
  from-scratch FCI engine reproduces pyscf FCI; the transcorrelated reference
  energy <HF|H_bar|HF> equals the pyscf MP2 energy exactly; the bare compact
  ground is just E_HF (Brillouin) while the transcorrelated compact ground
  recovers ~61% of the FCI correlation energy; tau=0 reproduces the bare result.
  This is exactly the role H_bar plays in pTC-MRSF-CIS (the production correlator
  is Ten-no's cusp-fixed geminal in place of T2; the downfolding machinery is
  identical).
- `cusp_convergence.py` -- **Phase-2b: the headline benefit, demonstrated
  exactly.** On Hooke's atom (real r12 cusp, exact via finite differences), an
  explicitly-correlated Gaussian basis {J g_k} (J with the correct coalescence
  cusp) converges to the exact energy far faster than the bare basis {g_k}: at
  n=4 the correlated basis is ~64x more accurate and already beats the bare basis
  at n=8. This is the basis-set-convergence acceleration pTC-MRSF-CIS inherits.
- `tc_nonsym_eig_test.F90` -- compiled Fortran validation of the actual
  `tc_nonsym_tda_eig` kernel (tau=0 gate vs LAPACK DSYEV; non-symmetric residual
  and biorthonormality), both at ~1e-14.
- `ptc_mrsf_cis.py` -- **the working method end-to-end on real molecular
  integrals (pyscf), no OpenQP integral engine needed.** On stretched H2/cc-pVDZ
  it builds bare MRSF-CIS (= CASCI(2,2): the ground state S0 as a response root
  plus the excited singlet/triplet manifold), folds external dynamic correlation
  in via H_bar = e^{-T2} H e^{T2}, downfolds onto the (2,2) space, and solves the
  non-Hermitian problem. S0 correlation recovered toward full FCI ~65%; tau=0
  reproduces bare MRSF-CIS; spectrum real. The closed-shell MP2 T2 dresses the
  singlet states (including S0) but not the triplet -- a concrete illustration of
  why pTC fixes the singlet AND triplet cusp conditions, which is what MRSF's dual
  targets require.

Integrals do not require a new Fortran geminal engine to obtain a *working*
method: `ptc_mrsf_cis.py` takes all 1e/2e integrals from pyscf and computes real
ground+excited states. The production OpenQP path would instead source bare
integrals from the Rys engine and add the geminal / normal-ordered 3-body
contributions; the pyscf route is the reference any Fortran implementation must
reproduce.

**Phase 2 — transcorrelated effective integrals. [PRIMITIVE LAYER DONE — validated]**
The native (pyscf-free) F12 geminal engine is implemented in
`source/modules/ptc_geminal.F90`: the closed-form Gaussian-geminal integral
`(ab|e^{-w r12^2}|cd)`, the intermediates X = geminal(2w), B = 8w^2·r2_geminal(2w),
V = <ab|e^{-w r12^2}/r12|cd> (mapped Gauss-Legendre; ->ERI as w->0), and the
STG-6G expansion of exp(-gamma r12). Validated standalone
(`tests/ptc_mrsf/prototype/ptc_geminal_test.F90`, build line in its header) against
pyscf-free oracles — all PASS: geminal(w->0)=overlap product (5e-9), finite-w vs
6D grid (~1e-13), r12^2·geminal vs grid (~1e-15), V(w=0)=Boys-F0 ERI (7e-16),
STG-6G vs exp(-r12) (~1e-4). Also compiles into liboqp (auto-globbed).
- **Remaining (the assembly):** `tc_build_eff_integrals` must loop these primitives
  over shell pairs of the real basis, contract with ROHF MOs, and DF/RI-factorize
  the correlation factor so the 3-body term (Phase 3) is never a dense O(N^6)
  tensor. Generalize the s-primitive kernels to higher angular momentum (or route
  through the existing Rys/Ishimura ERI engine for the 1/r12 parts).

**Phase 3 — normal-ordered higher-body terms.** `tc_normal_order_3body`:
fold 3-/4-body operators against the ROHF reference density into effective
1-/2-body integrals; keeps O(N^5).

**Phase 4 — integration. [INTEGRATION DONE + validated; H_bar molecular feed = remaining research]**
- **(a) tc control flag — DONE+validated.** `tc` (logical) + `tc_tau` (double) added to
  `tddft_parameters` (`source/types.F90`) and the C struct (`include/oqp.h`), plumbed
  through the Python input layer (`pyoqp/oqp/molecule/oqpdata.py`). `tc=True` under
  `[tdhf]` activates the pTC path end-to-end (Python -> cffi -> Fortran). The
  `OQP_PTC_MRSF` env seam is retained (OR'd in).
- **(b) route solver — DONE+validated.** `tdhf_mrsf_energy.F90` routes the reduced-space
  TDA diag (rpaeig call site) to `tc_nonsym_tda_eig` when enabled.
- **(d) tau=0 regression — PASS (bit-for-bit).** H2O/6-31G MRSF-s nstate=10: stock vs
  tc=True (and vs OQP_PTC_MRSF=1) at tau=0 give max|dE| = 0.0 eV. Fixtures
  `tests/ptc_mrsf/phase4_tau0_regression.inp`, `phase4_tcflag_tau0.inp`.
- **Mechanism validated end-to-end (exact oracle):** `tests/ptc_mrsf/prototype/tc_hubbard_test.F90`
  drives the production `tc_nonsym_tda_eig` on the Hubbard dimer: H_bar=J^{-1}HJ +
  non-Herm solve recover the EXACT spectrum (1.8e-15) and 100% of the correlation
  energy. Proves the transcorrelation -> non-Hermitian-solve chain.
- **Method works end-to-end on a molecule (native, pyscf-free).**
  `tests/ptc_mrsf/prototype/tc_mrsf_native_test.F90` runs the COMPLETE pTC-MRSF-CIS
  pipeline on H4/STO-3G built only from the native integral toolkit (ptc_geminal:
  S/T/V/ERI) + tc_nonsym_tda_eig: RHF, MP2, determinant FCI, H_bar=e^{-T2}He^{T2},
  compact downfold, non-Herm solve. Gated: H2 E_RHF vs textbook (4.5e-5),
  <HF|H_bar|HF>=E_MP2 (4e-16), 61.9% correlation recovered, tau=0 reproduces bare.
  This is a working pTC-MRSF-CIS producing real correlated numbers with the MP2-T2
  correlator.
- **(c) remaining production steps (engineering + one theory piece):**
  (i) swap the MP2-T2 correlator for the geminal (Ten-no pTC) effective integrals --
  the geminal engine (ptc_geminal: geminal/X/B/V + STG-6G + contraction) is built and
  validated; what's missing is the exact pTC effective-integral/normal-ordering FORMULA
  (the prototypes use the T2 proxy, so this is genuine theory work) and general
  angular-momentum geminal integrals (engine is s-only). (ii) inject H_bar into the
  LIVE OpenQP MRSF A.c (vs the standalone native demo) at
  `tdhf_mrsf_lib.F90 :: int2_mrsf_data_t_update` (cval/xval -> f3, ~L151-205), gated by
  the tau=0 bit-for-bit regression already in place.

## 6. Future: SA-CAS(4,4) substrate

Replace the single ROHF reference with a state-averaged CAS(4,4) to remove the
reference ambiguity entirely and present pTC with a genuine small-FCI solver
(cf. transcorrelated-multireference / TC-DMRG). Spin-pure by construction. Cost:
active-space selection + SA-CASSCF convergence, larger prefactor, new gradients.

## 7. Key code anchors (pilot)

- MRSF A·c build: `source/tdhf_mrsf_lib.F90 :: mrsfmntoia`, ERI feed
  `int2_mrsf_data_t_update` (f3 accumulation, ~L151–216).
- Symmetric assumption to replace: `source/tdhf_lib.F90 :: rpaeig` (TDA branch),
  `rparedms`, `rpavnorm`.
- MRSF solver loop / result storage: `source/modules/tdhf_mrsf_energy.F90`.
- Control struct: `source/types.F90 :: tddft_parameters` (C-bound; mirror in the
  Python ctypes layer when adding the `tc` flag).
