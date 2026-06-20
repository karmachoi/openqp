# QMRSF dual-pathway implementation design

Branch: `feat/qmrsf-dual-pathways` (worktree, off committed tip `2804a0ec` of `feat/mrsf-ensemble-reference`).
Companion theory write-up: Overleaf `6a3607a3...` (QMRSF: two pathways).
Created: 2026-06-20.

## 0. One-paragraph plan

The shared backbone already exists on `feat/mrsf-ensemble-reference`: the **determinant-basis,
weight-invariant collective response** `tdhf_mrsf_ensemble_sigma` in
`source/modules/tdhf_mrsf_energy.F90` (the "one consistent operator" in determinant basis;
spin-projection Ŝ² is the documented blocker, owned by the ensemble work). On top of that backbone
we add **two interchangeable dynamic-correlation layers**, selected by one input flag:

- **QMRSF-icPT2** (wavefunction picture): bare spin-pure DSF/RAS-SF(2) backbone + an
  internally-contracted external-Q second-order self-energy downfold
  `H_eff(E) = H_PP + Σ(E)`, `Σ(E) = H_PQ (E−H_QQ)^{-1} H_QP`, Dyall/Fink H₀, Hermitized
  (des Cloizeaux / NEVPT2-style partitioned denominator). O(N⁵) one-shot.
- **QMRSF-DK** (density-functional picture): KS orbitals + a dressed/frequency-dependent quadratic
  kernel `g_xc(ω)` on the closed-shell (0OS) double-spin-flip diagonal; correlation from the
  functional; **no** external-Q (avoids double counting). ~SCF + diagonalization cost.

Design rule for this branch (coordination): **add NEW files only; do not edit the shared
ensemble/spin-projection code** that the concurrent chat is changing. Wire-in happens behind a flag
once the backbone's Ŝ² projection lands.

## 1. Shared backbone (consumed, not modified here)

- Producer: `tdhf_mrsf_ensemble_sigma` (C-bound `tdhf_mrsf_ensemble_sigma`), flat-1D tagarray ABI
  (`mrsf_ens_perm`, `mrsf_ens_bvec` → `mrsf_ens_sigma`, `mrsf_ens_sigma_det`, `mrsf_ens_xm`).
- Driver: `pyoqp/oqp/library/single_point.py` runs the generalized Davidson `A c = ω S c` with the
  determinant configuration-overlap metric S.
- What both pathways need from it: the converged P-space eigenvectors {C_I^(state)}, the P-space
  effective Hamiltonian H_PP (or its action), and the active-space RDMs (1- and 2-RDM over O1..O4),
  which are cheap (4 electrons). These are the contraction inputs for icPT2 and the density inputs
  for the DK kernel.

## 2. Input/flag plumbing (new, additive)

Add to the tddft namelist (no change to existing defaults):
```
qmrsf_pathway = none | icpt2 | dk      ! default none (pure backbone, today's behavior)
qmrsf_0os_diag = backbone | seq | hve  ! 0OS diagonal source; default = backbone (bare, consistent)
qmrsf_icpt2_h0 = dyall | fink          ! zeroth-order H for the external-Q downfold
qmrsf_dk_gamma = <real>                ! dressed-kernel strength / freq parameter (DK)
```
Plumb through `pyoqp/oqp/.../single_point.py` → tagarray → new Fortran modules. The flag only routes;
it must leave `qmrsf_pathway=none` bit-identical to the current backbone.

## 3. Pathway I — QMRSF-icPT2 (new module `source/modules/tdhf_qmrsf_icpt2.F90`)

Stages (all post-backbone, operate on converged P vectors + RDMs):
1. Build the external-Q first-order interacting space (internally-contracted perturbers): the
   C→O, O→V, C→V and core/virtual double classes contracted against the P reference RDMs.
2. Form H₀ (Dyall: active two-body retained, inactive/virtual one-body) → orthonormalize perturbers
   (overlap metric S_Q), drop near-singular directions (intruder guard).
3. Assemble Σ via the partitioned resolvent on the orthonormalized perturbers; symmetrize
   (state-independent effective Hamiltonian) → Hermitian dressed H_PP.
4. Diagonalize the small dressed P matrix → corrected energies/vectors.

Validation: H₄ exact-vs-FCI (the existing external-Q pilot got this), then CBD / C₄H₆ dimer vs
XMS-CASPT2(4,4). Honest checks: first-order singles (orbital relaxation) magnitude on the
non-optimized ROHF reference; spin purity of the perturber space.

Novelty/benchmark to beat: MR-ADC(2) (Sokolov); contrast vs RAS-nSF-PT2 / RASCI(2) (determinant PT,
not internally contracted) and SF-ADC(2) (perturbative doubles; 0OS demoted to 0th order).

## 4. Pathway II — QMRSF-DK (new module `source/modules/tdhf_qmrsf_dk.F90`)

- Compute the dressed/quadratic kernel contribution `g_xc(ω)` for the six 0OS double-spin-flip
  diagonals on KS orbitals (frequency-dependent term restoring the double-excitation pole; an
  adiabatic kernel gives zero).
- Add it to the 0OS diagonal **only** (off-diagonals stay bare; spin-adaptation untouched). No
  external-Q. Keep one correlation source — assert no double counting with v_xc/f_xc already in the
  backbone diagonals.
- Well-definedness: the Maitra one-single/one-double dressing must be generalized to the coupled
  20-singlet block; document the chosen prescription and test it does not split degenerate multiplets.

Validation: show DK adds correlation the adiabatic MRSF kernel cannot (compare to backbone-only and
to icPT2 on the same geometries); PES/gradient continuity across an avoided crossing.

## 5. Shared validation gates (both pathways)
1. Spin purity: S=0↔S=1,2 off-blocks = 0 to 1e-10 after the existing spin-adaptation U.
2. Quintet self-consistency: dressed H recovers E₀ for the reference.
3. CAS-CI recovery: bare backbone == CAS(4,4) FCI in-space.
4. Continuity: smooth energy AND gradient across CBD square↔rectangular / C₄H₆-dimer scans.
5. Transferability: any DK/parameter set fixed once; leave-one-molecule-out, no refit.

## 6. Sequencing
1. [this branch] flag plumbing + two stub modules (compile, no behavior change at `none`).  ← next
2. QMRSF-icPT2 on H₄ → CBD/C₄H₆ vs XMS-CASPT2.
3. QMRSF-DK kernel derivation (see Overleaf Pathway II) → implement → calibrate against icPT2.
4. Optional explicit C→O/O→V/C→V channels (clean singles) once both cores validate.
