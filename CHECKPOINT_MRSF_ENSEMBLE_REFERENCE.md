# MRSF Ensemble-Reference Checkpoint

Created: 2026-06-16 05:38:51 KST
Updated: 2026-06-16 10:12 KST

Repository: `/Users/cheolhochoi/Documents/openqp-private`
Branch: `feat/mrsf-ensemble-reference`
Private remote branch: `origin/feat/mrsf-ensemble-reference`
Current implementation head before this checkpoint update:
`2879f1f9 Add MRSF reference scan harness`
Upstream base synced first: `1aabd750 Allow macOS LP64 BLAS and log build metadata (#209)`

Local branch note: upstream tracking was intentionally removed. Push with an
explicit refspec:

```bash
git push origin HEAD:feat/mrsf-ensemble-reference
```

## Latest Checkpoint: MRSF Reference SCF Stability Fix (ROHF identity)

Timestamp: 2026-06-16 10:12 KST

Code commit: `18cd26f1 fix: apply SCF stability safeguard to MRSF (tdhf) reference`
Run: `tools/_mrsf_response_smoke/o2_dissociation_stabfix_pes_20260616_100733/`

### Problem (user-reported)

The standalone triplet ROHF (`method=hf`) and the ROHF reference used by regular
MRSF (`method=tdhf`) solve the same SCF with identical `[scf]` settings, so they
must give the same energy.  They did not: at `R=3.20 A` standalone gave
`-149.4337` but the MRSF reference gave `-149.1725`.

### Root cause

`single_point.py::_run_scf` Stage-3 stability safeguard (the TRAH pass that
detects a DIIS-converged but *unstable* open-shell solution and relaxes it to the
lowest one) was gated `... and self.method == 'hf'`.  So MRSF (`tdhf`) skipped it
and built on the unstable c-DIIS ROHF.  The standalone HF log shows the TRAH pass
finding a Hessian eigenvalue `-1.05` ("unstable - escaping") and dropping
`-149.1725 -> -149.4337`; the MRSF log never ran it.

### Fix

Gate changed to `self.method in ('hf', 'tdhf')`.  The existing
restore-on-no-improvement logic reverts energy-invariant re-canonicalization
(stable minimum untouched); only a genuinely lower solution is kept.

### Result

- MRSF reference ROHF now equals standalone ROHF at all seven O2 points
  (`o2_rohf_identity.pdf`).
- Regular MRSF S0 is now smooth and much closer to CASSCF(4,4)
  (`o2_mrsf_s0_stability_fix.pdf`); the earlier "jumpy regular MRSF" was this bug:

```text
R(A)   MRSF_S0 before     MRSF_S0 after      CASSCF(4,4)
2.10   -149.238385        -149.361422        -149.356294
2.60   -149.179000        -149.355516        -149.375155
3.20   -149.345655        -149.351202        -149.370645
```

- Side effect: the ensemble ROHF spike at `R=2.60` is gone (it now shifts to a
  milder bump at `R=2.10`); the ensemble mean field still has its own rough point.

Validation: focused suite + MRSF gradient/TRAH tests = 68 tests OK.

### Newly exposed (next work, separate from the ROHF bug)

With the references corrected, the *ensemble block path* shows root-selection
fragility independent of the redundant-response SI fix:

- At `R=1.40` the auto ensemble uses a single `[8,9]` reference (same SCF as
  regular MRSF) yet its block MRSF returns the `+2.32` eV singlet instead of the
  `+0.28` eV state regular MRSF finds -- the adaptive trial-vector/root targeting
  in the ensemble block invocation picks a different root.
- At `R=2.60` the ensemble S0 over-stabilizes again (`-149.555`, below CASSCF).

These are in the per-reference block MRSF invocation (trial vectors / Davidson
root selection), not in `collect_state_interaction_response` and not in the
reference SCF.

## Latest Checkpoint: Redundant-Response Fix (lowest root = S0)

Timestamp: 2026-06-16 10:05 KST

Code commit: `bba3a7f3 fix: keep ensemble MRSF lowest root on physical S0`
(bundles the earlier uncommitted response-pinning + O2 scan-target work, whose
diffs were entangled in the same files as this fix).

### Problem addressed

User principle: the lowest ensemble root must be S0, so ensemble MRSF S0 should
track ordinary MRSF S0.  Before this fix the equal-weight ensemble dived to
`-150.31 Eh` at `R=3.20 A` (vs CASSCF(4,4) `-149.37`).

Root cause (two compounding effects over a near-degenerate frontier; at
`R=3.20 A` the O2 π* gap was `~1e-14 Eh`, so the six auto open pairs
`[8,9],[7,9],[6,9],[7,8],[6,8],[6,7]` are redundant rotations of the *same*
manifold):

1. **Non-variational alternative-reference responses.** Only the frontier pair
   `[8,9]` (ref 1) gives the physical S0 (block state `-149.3368`, matching
   ordinary MRSF `-149.3457`).  Forcing the open shell into the other pairs
   produces block "state 1" energies far below it (e.g. `[6,9]` gave
   `-150.178`, below CASSCF).  Selecting the *lowest* root picked these
   artifacts.
2. **Variational collapse of the energy-only state interaction.** The made-up
   off-diagonal `H_ij = 0.5 S_ij (E_i + E_j)` over a near-singular candidate
   overlap metric (min eigenvalue `~1e-6`) drove the root further down, with
   expansion coefficients `~4.6` (a normalized state has `|c| <= 1`).

The 2x2 model makes the collapse explicit: for two states `a,b` with overlap
`s`, the generalized eigenproblem gives lowest root
`(a+b)/2 - |a-b|/(2*sqrt(1-s^2))`, which `-> -inf` as `s -> 1`.

### Fix

`pyoqp/oqp/utils/mrsf_reference.py :: collect_state_interaction_response`
now guards the lowest root with two redundancy-aware steps before the
generalized eigensolve, and reports diagnostics:

- **Redundancy elimination** -- candidates are processed in physical priority
  order (anchor/frontier reference first, then weight, then energy) and a
  candidate is dropped when its common-basis overlap with an already-kept state
  exceeds `mrsf_ref.overlap_threshold` (previously parsed but unused).
- **Variational floor** -- the anchor reference (the frontier ROHF open pair,
  passed from `single_point.py`) defines the physical S0; candidates below it
  (minus a small tolerance) are rejected as non-variational artifacts.
- Canonical-orthogonalization metric threshold raised `1e-8 -> 1e-6`.
- New returned `redundancy` block (kept/dropped counts, anchor pair, floor,
  full-metric min eigenvalue / condition number); logged via `file_utils.py`
  (`PyOQP MRSF SI kept states / dropped redundant / dropped below floor /
  metric min eig`) and parsed by `tools/mrsf_reference_scan.py`.

`single_point.py` resolves the anchor from `metadata['frontier']
['current_open_pair']` (fallback: first reference) and passes it plus
`overlap_threshold` into the call.

### Result (O2, 6-31G; `tools/_mrsf_response_smoke/o2_dissociation_fixed_pes_20260616_095747/`)

```text
R(A)   regROHF       ensROHF       regMRSF S0     ensMRSF S0     CASSCF(4,4)    (ens-reg)mEh
1.10   -149.511371   -149.511371   -149.469709    -149.469709    -149.549733       +0.0
1.21   -149.528851   -149.527835   -149.488068    -149.488068    -149.593794       +0.0
1.40   -149.514736   -149.479782   -149.460857    -149.460857    -149.599728       -0.0
1.70   -149.465003   -149.465003   -149.453141    -149.437830    -149.574059      +15.3
2.10   -149.444742   -149.382079   -149.238385    -149.375971    -149.356294     -137.6
2.60   -149.438409   -148.975645   -149.179000    -149.202434    -149.375155      -23.4
3.20   -149.433678   -149.435529   -149.345655    -149.336780    -149.370645       +8.9
```

- The `-150.31` collapse at `R=3.20` is gone; ensemble S0 now equals the anchor
  S0 (`+0.0987` excitation -> `-149.3368`), within `~9 mEh` of ordinary MRSF and
  closer to CASSCF at `R=2.10`.
- At the multi-reference points the log confirms `kept=1`, `dropped below
  floor=5`, `max abs offdiag = 0.0`, anchor coefficient `1.0`.

### Plots requested

- `o2_s0_before_after_fix.png/pdf` -- before/after S0 vs regular MRSF vs CASSCF.
- `o2_rohf_vs_ensemble_rohf.png/pdf` -- regular triplet ROHF vs ensemble ROHF,
  plus the ROHF reference ordinary MRSF actually uses.

Diagnostic finding from the ROHF plot: standalone triplet ROHF is the smoothest
reference; ordinary MRSF's own ROHF reference gets stuck high at stretched
bonds; the ensemble ROHF recovers the good solution at `R=3.20` but still has a
spurious spike at `R=2.60` (`-148.976`).  So the *response* collapse is fixed,
but a **mean-field-level** continuity problem at `R=2.60` remains for later.

### Validation

```bash
python3 -m unittest \
  tests.test_mrsf_reference_scan_tool tests.test_mrsf_reference_metadata \
  tests.test_mrsf_reference_scf tests.test_umrsf_energy_regression \
  tests.test_single_point_scf_fallback
# Ran 60 tests OK  (incl. 2 new SI guard tests)
```

Installed package was refreshed by copying the three edited Python files into
`site-packages/oqp` (Python-only change; no Fortran rebuild needed).

### Next

- Stabilize the ensemble *mean field* (the `R=2.60` ROHF spike) -- the remaining
  discontinuity is now at SCF/root-selection, not in the response.
- The variational floor is anchored on the frontier pair; revisit if a genuine
  avoided crossing makes another reference physically lower (overlap tracking).

## Previous Checkpoint: O2/CASSCF Comparison and Response Pinning

Timestamp: 2026-06-16 09:45:48 KST

Code repository state:

- Checkout: `/Users/cheolhochoi/Documents/openqp-private`
- Branch: `feat/mrsf-ensemble-reference`
- Last pushed code commit: `699c877c feat: expose ensemble MRSF off-diagonal coupling`
- Current working tree has uncommitted code changes and generated smoke outputs.
- Overleaf manuscript was updated separately and pushed to Overleaf `main` at
  `f0a1cc2 docs: add O2 MRSF CASSCF comparison`.

Important local code changes since `699c877c`:

- `pyoqp/oqp/utils/mrsf_reference.py`
  - Added `freeze_mrsf_reference_config`.
  - Purpose: once ensemble SCF has applied an auto-selected reference set, the
    response metadata is pinned to that same applied set instead of recomputing
    references from the final orbitals and accidentally collapsing back to one
    reference.

- `pyoqp/oqp/library/single_point.py`
  - `_update_mrsf_reference_metadata` now freezes response references to the
    SCF-applied open pairs when prior SCF metadata exists.
  - Logs `response_reference_source = scf_applied_open_pairs` and a warning
    that response references are pinned to the ensemble-SCF reference set.

- `tools/mrsf_reference_scan.py`
  - Added `o2_dissociation` scan target.
  - Added `equal_block` variant using existing
    `mrsf_ref.coupling=block_diagonal` as a control.

- Tests updated:
  - `tests/test_mrsf_reference_metadata.py`
  - `tests/test_mrsf_reference_scf.py`
  - `tests/test_mrsf_reference_scan_tool.py`

Focused validation after these local changes:

```bash
python3 -m unittest \
  tests/test_mrsf_reference_scan_tool.py \
  tests/test_mrsf_reference_metadata.py \
  tests/test_mrsf_reference_scf.py
```

Result:

```text
Ran 43 tests
OK
```

The installed package was refreshed with `pip install .` using Homebrew
GCC/GFortran compilers.  Do not use clang for this branch, and do not remove or
disable the external build cache.

Working install command used:

```bash
CC=/opt/homebrew/bin/gcc-15 \
CXX=/opt/homebrew/bin/g++-15 \
FC=/opt/homebrew/bin/gfortran-15 \
CMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
CMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
CMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 \
pip install .
```

O2 validation artifacts:

- Auto ensemble comparison:
  `tools/_mrsf_response_smoke/o2_dissociation_auto_pes_20260616_093241/`
- Manual six-reference block/off-diagonal control:
  `tools/_mrsf_response_smoke/o2_dissociation_manual6_control_pes_20260616_092700/`
- CASSCF comparison:
  `tools/_mrsf_response_smoke/o2_mrsf_casscf44_compare_20260616_093614/`

Key comparison report:

```text
tools/_mrsf_response_smoke/o2_mrsf_casscf44_compare_20260616_093614/comparison_report.md
```

Main O2 result table:

```text
R (A)  regular MRSF Eh    equal auto Eh      equal forced-6 Eh   CASSCF(4,4) Eh
1.10   -149.4697085100   -149.4697085100   -149.6460348300   -149.5497331264
1.21   -149.4880682400   -149.4880682400   -149.6534740900   -149.5937940140
1.40   -149.4428861800   -149.4428861800   -149.5852219300   -149.5997279097
1.70   -149.3208086900   -149.3228609200   -149.5841396700   -149.5740585699
2.10   -149.2383850700   -149.3759711200   -149.3535165800   -149.3562944200
2.60   -149.1790001800   -149.4033668200   -149.3787234400   -149.3751553008
3.20   -149.3456546600   -150.3142956400   -150.2170868500   -149.3706446669
```

Interpretation:

- CASSCF(4,4) singlet, computed with PySCF 2.13.0 in 6-31G with active
  orbitals `[6,7,8,9]`, converged at all seven O2 bond lengths with
  near-zero `<S^2>`.
- CASSCF(4,4) does not show the huge long-bond collapse.
- Regular MRSF is also much milder at stretched O2.
- Equal-weight ensemble MRSF, both auto and forced-six, dives strongly at
  stretched O2.  At `R=3.20 A`, auto equal gives `-150.31429564 Eh` while
  CASSCF(4,4) gives `-149.37064467 Eh`.
- Therefore the current equal-weight ensemble prototype has an ensemble
  SCF/root-selection problem for O2 dissociation.  This is not a physical
  CASSCF-like cusp.
- The response pinning fix is still useful: it prevents auto mode from using a
  six-reference ensemble SCF but then silently recomputing one response block
  after SCF.  However, pinning alone does not solve the physical/root-selection
  issue.

Overleaf update:

- Repo: `/Users/cheolhochoi/Documents/Manuscripts/overleaf-6a30598d`
- Remote: `https://git@git.overleaf.com/6a30598df69f710f964ccab2`
- Commit pushed: `f0a1cc2 docs: add O2 MRSF CASSCF comparison`
- Added figures:
  - `figures/o2_mrsf_casscf44_main_compare.pdf`
  - `figures/o2_mrsf_casscf44_equal_compare.pdf`
- Local validation:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Result: successful 11-page build.

Overleaf branch policy note:

- A branch push was attempted first:

```bash
git push origin HEAD:refs/heads/codex/o2-mrsf-casscf44-overleaf
```

- Overleaf rejected it with `wrong branch` and `Please use the main branch`.
- The validated commit was then pushed to `main`.

Current next technical step:

1. Stabilize the ensemble reference/root-selection criterion before doing
   gradients, RT-MRSF, EKT, NAC, or property work.
2. The O2 result suggests lowest-root selection in the current ensemble
   state-interaction space is unsafe.
3. Candidate fixes to evaluate next:
   - select/track the physically relevant root by overlap with ordinary MRSF or
     CASSCF-like reference character,
   - regularize or truncate nearly linearly dependent response components,
   - revise the ensemble SCF target so the stretched O2 reference does not
     collapse into the over-stabilized root,
   - compare spin-adapted CASSCF state characters against MRSF block
     components before changing the response equation.

Do not start gradients, RT-MRSF, EKT, UMRSF, NAC, or Davidson work yet.

## Scope

This branch prototypes an ensemble-reference mean-field layer for ambiguous
ROHF references in MRSF-TDDFT. The target problem is discontinuous PES behavior
when more than one triplet ROHF open-shell configuration is plausible, such as
near-degenerate frontier orbitals, benzene-like degeneracy, or weak dimers.

The implemented part is ensemble-reference SCF. The coupled ensemble MRSF
linear-response solver is deliberately not implemented yet and is guarded by a
clear `NotImplementedError`.

This branch does not implement EKT changes and does not change Davidson.

## Commits

- `0d97003e Add ensemble-reference MRSF SCF prototype`
- `ba8f333f Add auto MRSF reference pair selection`
- `6a025a38 Add gap-softmax MRSF reference weights`
- `2879f1f9 Add MRSF reference scan harness`

## Implemented Input

New section:

```ini
[mrsf_ref]
mode=off | diagnostic | state_average
open_pairs=auto | 5:6;4:7
weights=equal | 0.5,0.5
weight_temperature=0.05
max_refs=2
gap_threshold=0.01
overlap_threshold=0.85
strict=False
```

Manual mode:

```ini
[mrsf_ref]
mode=state_average
open_pairs=5:6;4:7
weights=0.5,0.5
```

Auto mode:

```ini
[mrsf_ref]
mode=state_average
open_pairs=auto
weights=equal
max_refs=2
```

Gap-softmax weights:

```ini
[mrsf_ref]
mode=state_average
open_pairs=auto
weights=gap_softmax
weight_temperature=0.05
max_refs=2
```

For a triplet ROHF reference with `nelec_A = nelec_B + 2`, auto mode selects a
frontier-window ensemble. It prioritizes:

1. the current ROHF open pair, normally `[nelec_B + 1, nelec_B + 2]`
2. the balanced frontier alternative, normally `[nelec_B, nelec_A + 1]`
3. additional frontier-window pairs up to `max_refs`

For the H2O triplet 6-31G smoke test, auto mode selected:

```text
[[5, 6], [4, 7]]
```

## Main Files Changed

- `pyoqp/oqp/utils/mrsf_reference.py`
  - Parses `[mrsf_ref]`.
  - Builds manual and auto candidate open-pair lists.
  - Builds weighted triplet ROHF ensemble metadata.
  - Produces dense alpha/beta occupation vectors.
  - Resolves `weights=equal`, explicit numeric weights, and
    `weights=gap_softmax`.

- `pyoqp/oqp/library/single_point.py`
  - Stages `OQP::mrsf_ref_occ_a` and `OQP::mrsf_ref_occ_b` before SCF.
  - Logs reference metadata.
  - Records `scf.applied_open_pairs`, `scf.applied_weights`, and
    `scf.applied_weight_model` because softmax diagnostics can be recomputed
    after SCF from converged MOs.
  - Blocks unfinished coupled response for `mode=state_average`.

- `source/guess.F90`
  - Allows explicit spin occupation vectors when building densities.

- `source/scf.F90`
  - Reads and validates ensemble-reference occupation tags.
  - Rebuilds ROHF densities from fixed fractional occupations.
  - Uses an ensemble-compatible ROHF orbital-update Fock path.

- `source/tagarray_driver.F90`
  - Adds `OQP::mrsf_ref_occ_a` and `OQP::mrsf_ref_occ_b`.

- `pyoqp/oqp/utils/input_checker.py`
  - Validates `[mrsf_ref]`.
  - Allows both manual and auto `state_average`.
  - Rejects pFON with `state_average`.

- `pyoqp/oqp/utils/file_utils.py`
  - Logs pair selection mode and ensemble metadata.

- `tools/mrsf_reference_scan.py`
  - Generates H2O triplet OH-stretch scan inputs.
  - Runs `rohf`, `equal`, and `gap_softmax` variants.
  - Parses SCF energy, convergence/escalation status, selected open pairs,
    reference weights, and SCF-applied weights.
  - Writes ignored scratch summaries under
    `tools/_mrsf_reference_scan_scratch/`.

- `tests/test_mrsf_reference_metadata.py`
- `tests/test_mrsf_reference_scf.py`
- `tests/test_mrsf_reference_scan_tool.py`

## Build Command

Use `pip install .`. On this macOS checkout, CMake selected AppleClang by
default and failed OpenMP discovery, so the working install command was:

```bash
CC=/opt/homebrew/bin/gcc-15 \
CXX=/opt/homebrew/bin/g++-15 \
FC=/opt/homebrew/bin/gfortran-15 \
CMAKE_ARGS="-DLINALG_LIB_INT64=OFF -DENABLE_OPENTRAH=OFF" \
pip install .
```

Do not force `OQP_REUSE_EXTERNALS=OFF`; that discards the external cache and
makes the build much slower.

## Validation Passed

Focused tests:

```bash
python3 -m unittest \
  tests.test_mrsf_reference_scan_tool \
  tests.test_mrsf_reference_metadata \
  tests.test_mrsf_reference_scf \
  tests.test_umrsf_energy_regression \
  tests.test_symmetry_parser_checker \
  tests.test_single_point_scf_fallback
```

Result:

```text
Ran 47 tests
OK
```

Final `pip install .` succeeded and installed `OpenQP 0.1.0`.

Installed-package auto smoke:

```bash
openqp --nompi --omp 1 tmp_mrsf_ref_auto_smoke.inp
```

Expected behavior was observed:

- input check passed with one warning about unfinished coupled response
- pair selection: `auto/frontier_window`
- selected open pairs: `[[5, 6], [4, 7]]`
- ensemble occupations were applied
- SCF converged in 11 iterations
- final ROHF energy: `-75.4080504337`
- run then stopped at expected `NotImplementedError` before coupled response

Temporary smoke input/log files were removed after validation.

Installed-package gap-softmax smoke:

```bash
openqp --nompi --omp 1 tmp_mrsf_ref_gap_softmax_smoke.inp
```

Expected behavior was observed:

- input check passed with one warning about unfinished coupled response
- pair selection: `auto/frontier_window`
- selected open pairs: `[[5, 6], [4, 7]]`
- weight model: `gap_softmax`
- `weight_temperature=0.05`
- SCF-applied weights from the pre-SCF proxy:
  `[0.9873932056032152, 0.012606794396784957]`
- post-SCF diagnostic weights from converged MO proxies:
  `[0.8866801587843114, 0.11331984121568874]`
- SCF converged in 23 iterations
- final ROHF energy: `-75.6144482009`
- run then stopped at expected `NotImplementedError` before coupled response

H2O triplet OH-stretch scan:

```bash
python3 tools/mrsf_reference_scan.py --points 0.99,1.00,1.01
```

Result summary path from the validation run:

```text
tools/_mrsf_reference_scan_scratch/h2o_triplet_20260616_055653/summary.csv
```

Observed scan table:

```text
scale  variant       status                   energy          iter  converged  escalated  applied_weights
0.99   rohf          ok                       -75.7146014183  1     yes        no         -
0.99   equal         expected_response_guard  -75.4028526940  11    yes        no         [0.5, 0.5]
0.99   gap_softmax   expected_response_guard  -75.5425319983  25    yes        yes        [0.9880432004455875, 0.011956799554412643]
1.00   rohf          ok                       -75.7192307219  1     yes        no         -
1.00   equal         expected_response_guard  -75.4080504337  11    yes        no         [0.5, 0.5]
1.00   gap_softmax   expected_response_guard  -75.6144482009  18    yes        yes        [0.9873932056032152, 0.012606794396784957]
1.01   rohf          ok                       -75.7235134936  1     yes        no         -
1.01   equal         expected_response_guard  -75.4129042826  11    yes        no         [0.5, 0.5]
1.01   gap_softmax   expected_response_guard  -75.6189269677  19    yes        yes        [0.9867152385230168, 0.013284761476983297]
```

Immediate interpretation:

- Auto pair selection remained stable as `[[5, 6], [4, 7]]`.
- Equal-weight ensemble SCF was smooth across these three points and did not
  require SCF escalation.
- Gap-softmax weights changed smoothly, but those runs escalated to SOSCF and
  showed a large energy step from scale 0.99 to 1.00. Treat gap-softmax as a
  useful diagnostic but not yet the preferred continuity model.

## Known Limitations

- `mode=state_average` currently supports ensemble-reference SCF only.
- Coupled block MRSF response over multiple references is not implemented.
- `weights=equal`, explicit numeric weights, and `weights=gap_softmax` are
  implemented. Gap-softmax weights are fixed before SCF from the available
  orbital-energy proxy; post-SCF diagnostics may report different proxy weights
  from converged MOs, so check `scf.applied_weights` for what was actually used.
- `gap_threshold` currently contributes diagnostics/warnings; it is not yet a
  hard filter for auto candidate selection.
- MO overlap tracking across geometry is not implemented yet.
- The scan harness currently implements only the H2O triplet OH-stretch sanity
  target. Ethylene torsion and benzene-like degeneracy scans are next.
- Full `python3 -m unittest discover -s tests` is not a useful pass/fail signal
  in this checkout; it entered broader native tests and crashed with a bus error
  in `append_shell`. The focused suite and installed smoke above are the current
  validation basis for this branch.
- The sandbox prints an OpenMPI TCP bind warning even with local runs. Use
  `openqp --nompi` for these small smoke tests.

## Next Steps

1. Add the first real PES continuity test:
   - ethylene torsion
   - benzene or benzene-like degenerate pi system
   - weak dimer frontier-degeneracy case

2. Add MO-overlap tracking:
   - preserve diabatic-like open-pair labels along geometry
   - use overlaps to avoid arbitrary pair relabeling

3. Add more weight optimization models:
   - later: variational or entropy-regularized weight optimization
   - decide whether gap-softmax should remain fixed from initial MOs or be made
     self-consistent inside SCF

4. Only after SCF continuity is demonstrated, implement coupled ensemble MRSF
   response:
   - build one spin-flip block per selected ROHF reference
   - add coupling between reference blocks
   - remove the current `NotImplementedError` guard when the solver is real

## Minimal Restart Commands

```bash
cd /Users/cheolhochoi/Documents/openqp-private
git switch feat/mrsf-ensemble-reference
git status --short --branch
git log --oneline --decorate -5
```

Then rebuild if needed with the `pip install .` command above and run the
focused test command above.
