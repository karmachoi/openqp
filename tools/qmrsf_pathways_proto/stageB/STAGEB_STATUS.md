# QMRSF Stage B — OpenQP live integration: status

Goal: wire the validated standalone Fortran (backbone + CSF + icPT2 + DK, all matching NumPy to
1e-14) into OpenQP with live quintet ROKS integrals, full build, real-molecule benchmarks.

## Done
- **Quintet ROHF reference runs in OpenQP** (the reference the method is built on).
  `h4_quintet_rohf.inp` (linear H4, STO-3G, `[scf] multiplicity=5 type=rohf`) converges:
  `Final ROHF energy -1.5198126991` (1 iter, high-spin single determinant). mult=5 confirmed.
  This validates the quintet (S=2, four singly-occupied orbitals) reference path end-to-end,
  **runs to exit 0** (full property/Mulliken output produced).
- **Post-SCF BLAS abort FIXED** (commit 1998a386). The TRAH stability safeguard built an
  orbital-rotation Hessian over an EMPTY occ-virt block for the fully-occupied quintet
  (nbf==n_occ), issuing `cblas_dgemm LDA=0` -> Accelerate hard-abort. Guard added in
  `single_point.py`: skip the safeguard when `nbf <= nocc_max` ("no virtual orbitals ->
  no orbital rotations; reference trivially stable"). Log confirms the skip + clean exit.
- **QMRSF plumbing made rebuild-tolerant** (same commit): `energy_func` dispatch + `qmrsf_*`
  setters degrade gracefully (getattr/try-except) against a liboqp predating the QMRSF symbols,
  so pyoqp runs pre-rebuild; an explicit QMRSF request still raises a clear "rebuild" error.
- **AO->MO active-space transform BUILT + VALIDATED end-to-end (the headline).**
  `source/modules/qmrsf_ao2mo.F90` (`qmrsf_active_integrals`) produces h_eff, (pq|rs), E_core
  over the four frontier active orbitals by REUSING the int2 Coulomb digestor (probe densities,
  one batched sweep, no AO tensor) -- OpenQP had no 4-index ao2mo before this.
  `source/modules/qmrsf_cas.F90` (`qmrsf_cas_solve`) ports the validated CAS(4,4) backbone.
  `tdhf_qmrsf_icpt2` now runs the real chain: active space (SOMOs) -> transform -> CAS -> diag.
  Standalone dispatch added: `[tdhf] type=qmrsf_icpt2` runs SCF (quintet ROHF) then the routine
  (no MRSF Davidson) via `Calculator.qmrsf_standalone` (single_point.py) + TDHF_TYPES.
  **GATE (route_a_oracle.py, pyscf-free closed-form H4/STO-3G oracle) ALL PASS:**
    C^T S C - I = 3.6e-9 ; h_act = 7.8e-9 ; (pq|rs) = 2.9e-9 ; CAS spectrum = 1.5e-8 ;
    E_core exact ; CAS(4,4)=FCI ground total -2.1026 (= oracle, < quintet ROHF -1.5198).
  Convention pinned empirically: probe density packed WITHOUT off-diagonal doubling
  (matches pack_matrix/dtrttp); raw int2_rhf consumer -> true J via fock_jk post-scaling
  (halve all, double diagonal); `scale_coulomb=1, scale_exchange=0` for pure J.
  Build: standalone ninja (gcc-15, USE_LIBINT=OFF, INT64=OFF), liboqp staged at a private
  OPENQP_ROOT=/tmp/qmrsf_root (worktree oqp.h has the QMRSF C decls) -- shared install untouched.

- **LIVE external-Q icPT2 downfold DONE + validated (the icPT2 payoff).** The full pathway
  now runs live: quintet ROHF -> int2-reuse AO->MO transform (frozen-core window of
  active+virtual MOs) -> brute-force determinant CI -> internally-contracted external-Q EN
  des-Cloizeaux multistate downfold. Modules: `qmrsf_icpt2_downfold.F90` (kernel),
  `qmrsf_icpt2_engine.F90` (det-CI perturber generation), driven by `tdhf_qmrsf_icpt2`.
  `[tdhf] type=qmrsf_icpt2` end-to-end. **Gates ALL PASS:**
    - downfold algebra (standalone vs NumPy multistate): EN 5.7e-14, Dyall 7.8e-14;
    - full pipeline (standalone, MO ints->det-CI->EN): 5.7e-14;
    - LIVE H4/6-31G (gate_icpt2_full.py vs NumPy on the dumped window integrals):
      bare CAS 4.0e-15, icPT2 EN 3.5e-13. ndet=784, P=36, Q=748, H_eff Hermitian to 0.
      H4/6-31G ground: CAS -2.15174 -> icPT2 -2.41298 (recovers 0.261 Ha dynamic corr).
    - H4/STO-3G regression: norb_w==nact, Q empty -> CAS=FCI -2.10261 (unchanged).
  Also fixed a second quintet BLAS abort: the POST-SCF TRAH stability check builds an empty
  occ-virt Hessian block when a spin channel is empty (nelec_B=0, M_s=+2). Guard extended in
  single_point.py: skip stability when nbf<=nocc OR min(nelec_A,nelec_B)==0. Quintet inputs
  set [scf] converger_type=diis (auto-manager else escalates to the empty-block TRAH).
  Brute-force det space is guarded (MAXDET=6e5); large windows fall back to CASCI with a note
  (the contracted/scalable perturber engine + Dyall denominators are the remaining work).

- **FIRST REAL-MOLECULE RESULTS: cyclobutadiene (CBD), rectangular, quintet ROHF.**
  The full pathway runs on a polyatomic with a real frozen core: quintet ROHF (28 e-,
  ncore=12) -> CAS(4,4) pi window -> contracted external-Q EN+Dyall downfold -> log + JSON.
    * CBD/STO-3G  (window 12, nvirt 8,  nQ=1024):  S0->S1 = 2.165 eV (EN) / 2.163 (Dyall);
      EN and Dyall agree (few virtuals, no intruder). ref quintet ROHF -151.4921.
    * CBD/6-31G   (window 32, nvirt 28, nQ=245,980): contracted engine handles ~2.5e5
      perturbers with NO FCI list (brute force ~1e9 dets, infeasible). Dyall S0->Sn =
      1.59 / 4.51 / 4.58 / 4.88 / 6.53 eV (physical). **EN catastrophically intruder-broken**
      (state-0 "correction" -7.6 Ha): with 2.5e5 Q states some are near-degenerate with P
      roots -> near-zero EN denominators. Textbook Epstein-Nesbet intruder catastrophe ->
      concrete evidence that **Dyall is the production denominator** (manuscript's choice).
  PRODUCTION NOTE: default to Dyall; EN needs a level shift / intruder regularization
  (proto has a level_shift hook; live engine uses 0). Inputs: cbd_quintet_icpt2{,_631g}.inp.

## Open issues / barriers
1. **Molden writer IndexError** (downstream, unrelated): `write_basis` -> `molden_bas[sh_at]`
   tuple index out of range for this system. Worked around with `save_molden=False`. Not on
   the QMRSF path; revisit only if molden output is wanted for QMRSF runs.
2. **No AO->MO integral transform exists (THE missing piece).** OpenQP is integral-direct;
   MRSF only ever does density -> J/K via `int2`. There is no `ao2mo`/`moints`. The active-space
   one- and two-electron integrals h_act, (pq|rs) over the active MOs of the quintet reference
   must be produced by a NEW Fortran routine.
   **DECISION (user, 2026-06-20): build it in Fortran; do NOT dump AO ints to python.**
   Route B (AO dump + numpy bridge) is DROPPED — it was only ever a throwaway validation hack
   and a python AO transform is not a production path.
   **Route A (adopted) — active-space transform that REUSES the validated `int2` digestor:**
     - 1e:  h_act(p,q) = sum_{mu,nu} C_{mu,p} Hcore_{mu,nu} C_{nu,q}  (dense, trivial).
     - 2e:  for each active pair (r,s): build AO density  D^{rs}_{lam,sig} = C_{lam,r} C_{sig,s},
            run it through int2 to get the Coulomb matrix  G^{rs}_{mu,nu} = sum_{lam,sig}(mu nu|lam sig) D^{rs}_{lam,sig},
            then  (pq|rs) = sum_{mu,nu} C_{mu,p} G^{rs}_{mu,nu} C_{nu,q}.
       => O(n_act^2) screened J-builds, NO AO 2e tensor ever stored. Same machinery MRSF uses.
       (n_act here = the 4 frontier MOs; the bulk virtuals/core are never transformed.)
   Once h_act, eri_act are assembled they feed the already-validated backbone/CSF/icPT2/DK
   (the standalone Fortran reads exactly this; in-process it is passed directly).

## Interface already in place
The validated standalone Fortran reads active integrals from a text file
(`qmrsf_cas_ref.dat`: nact, h_act(4,4), eri_act(4,4,4,4) chemist). Once route (A) or (B)
produces those integrals from the live quintet reference, the backbone + CSF + icPT2 + DK
plug in unchanged (all already validated to 1e-14 vs NumPy).

## Next
- Decide route (A) vs (B) for integral extraction (A = production, heavy; B = fast bridge).
- Fix the post-SCF empty-block dgemm for high-spin references.
- Then: H4/STO-3G CAS(4,4) = FCI gate; CBD square/rectangular; benzene vs XMS-CASPT2.
