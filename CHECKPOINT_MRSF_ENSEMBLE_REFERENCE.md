# MRSF Ensemble-Reference Checkpoint

Created: 2026-06-16 05:38:51 KST
Updated: 2026-06-16 05:48:05 KST

Repository: `/Users/cheolhochoi/Documents/openqp-private`
Branch: `feat/mrsf-ensemble-reference`
Private remote branch: `origin/feat/mrsf-ensemble-reference`
Current head: `6a025a38 Add gap-softmax MRSF reference weights`
Upstream base synced first: `1aabd750 Allow macOS LP64 BLAS and log build metadata (#209)`

Local branch note: upstream tracking was intentionally removed. Push with an
explicit refspec:

```bash
git push origin HEAD:feat/mrsf-ensemble-reference
```

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

- `tests/test_mrsf_reference_metadata.py`
- `tests/test_mrsf_reference_scf.py`

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
  tests.test_mrsf_reference_metadata \
  tests.test_mrsf_reference_scf \
  tests.test_umrsf_energy_regression \
  tests.test_symmetry_parser_checker \
  tests.test_single_point_scf_fallback
```

Result:

```text
Ran 44 tests
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
- Full `python3 -m unittest discover -s tests` is not a useful pass/fail signal
  in this checkout; it entered broader native tests and crashed with a bus error
  in `append_shell`. The focused suite and installed smoke above are the current
  validation basis for this branch.
- The sandbox prints an OpenMPI TCP bind warning even with local runs. Use
  `openqp --nompi` for these small smoke tests.

## Next Steps

1. Add PES continuity tests:
   - small H2O triplet scan first
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
