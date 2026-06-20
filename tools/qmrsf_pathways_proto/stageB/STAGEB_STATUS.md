# QMRSF Stage B — OpenQP live integration: status

Goal: wire the validated standalone Fortran (backbone + CSF + icPT2 + DK, all matching NumPy to
1e-14) into OpenQP with live quintet ROKS integrals, full build, real-molecule benchmarks.

## Done
- **Quintet ROHF reference runs in OpenQP** (the reference the method is built on).
  `h4_quintet_rohf.inp` (linear H4, STO-3G, `[scf] multiplicity=5 type=rohf`) converges:
  `Final ROHF energy -1.5198126991` (1 iter, high-spin single determinant). mult=5 confirmed.
  This validates the quintet (S=2, four singly-occupied orbitals) reference path end-to-end.

## Open issues / barriers
1. **Post-SCF BLAS error** `cblas_dgemm LDA=0 M=0` fires AFTER `SCF converged with diis`
   (a property/energy-wrapper step hits a zero-dimension block for the all-singly-occupied
   quintet). SCF result is fine; the run aborts before completing. Fix: guard the empty-block
   dgemm in the post-SCF path for high-spin (n_open = nbf or empty virtual/closed blocks).
2. **No off-the-shelf MO-integral transform.** OpenQP is integral-direct; MRSF uses `int2`
   (density -> J/K). There is no `ao2mo`/`moints`. The active-space one- and two-electron
   integrals h_act(4,4), (pq|rs) over the four active MOs of the quintet reference must be
   produced by a NEW routine. Two routes:
     (A) int2-based active-MO transform (Fortran): feed active-orbital outer-product densities
         D_pq = c_p c_q^T through the int2 digestor to assemble (pq|rs); reuse int2_compute_t as
         tdhf_mrsf_energy does. O(N^4), screened. The production route.
     (B) AO-integral dump + numpy transform (pyscf-free bridge): if OpenQP can emit AO 1e/2e
         integrals (conventional mode), do the AO->active-MO transform in numpy and write the
         .dat the Fortran already reads. Faster to validate; not production.

## Interface already in place
The validated standalone Fortran reads active integrals from a text file
(`qmrsf_cas_ref.dat`: nact, h_act(4,4), eri_act(4,4,4,4) chemist). Once route (A) or (B)
produces those integrals from the live quintet reference, the backbone + CSF + icPT2 + DK
plug in unchanged (all already validated to 1e-14 vs NumPy).

## Next
- Decide route (A) vs (B) for integral extraction (A = production, heavy; B = fast bridge).
- Fix the post-SCF empty-block dgemm for high-spin references.
- Then: H4/STO-3G CAS(4,4) = FCI gate; CBD square/rectangular; benzene vs XMS-CASPT2.
