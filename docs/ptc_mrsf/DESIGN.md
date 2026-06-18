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
- Validated as a NumPy reference: `docs/ptc_mrsf/prototype/nonsym_tda_eig.py`.
  **Acceptance gate (met):** with tau=0 the non-symmetric solver reproduces the
  symmetric `rpaeig` TDA eigenpairs to machine precision (max|dE| ~ 7e-15), and the
  non-symmetric case yields a real spectrum with biorthonormal vectors
  (residuals ~1e-14).
- Not yet wired into the build/hot path.

**Phase 2 — transcorrelated effective integrals.** `tc_build_eff_integrals`:
DF/RI 1-/2-body effective integrals from the correlation factor + ROHF orbitals
(the N^5 step). Never materialize a dense 3-body tensor.

**Phase 3 — normal-ordered higher-body terms.** `tc_normal_order_3body`:
fold 3-/4-body operators against the ROHF reference density into effective
1-/2-body integrals; keeps O(N^5).

**Phase 4 — integration.** Add a `tc` control flag (off by default; plumb through
the Python input layer and the `tddft_parameters` C-bound struct), route the MRSF
solver loop in `tdhf_mrsf_energy.F90` to `tc_nonsym_tda_eig` when enabled, feed
H_bar through the existing AO-driven `f3` accumulation
(`tdhf_mrsf_lib.F90 :: int2_mrsf_data_t_update`). Regression: with the flag on but
tau=0, reproduce stock MRSF-CIS bit-for-bit.

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
