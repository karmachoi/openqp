# ISPHER basis option plan

## Semantics

OpenQP now recognizes a GAMESS-compatible `[input]` keyword named `ispher`:

- `ispher=-1`: Cartesian basis functions/SALC, matching the historical OpenQP basis path and the GAMESS default Cartesian-compatible variational space.
- `ispher=0`: accepted compatibility mode. GAMESS builds SALCs from spherical harmonics while retaining Cartesian contaminants; OpenQP currently has no separate SALC-only population-analysis distinction, so this is documented and logged as a Cartesian-equivalent compatibility mode.
- `ispher=1`: requested true pure/spherical-harmonic basis functions. This would reduce shell sizes such as d 6->5, f 10->7, and g 15->9. The current native OpenQP integral/basis backend is Cartesian, so this mode is rejected explicitly instead of silently pretending to run a pure basis.

The default remains `ispher=-1`, preserving current Cartesian-compatible behavior.

## Current basis architecture

- `pyoqp/oqp/molecule/oqpdata.py` owns the input schema and parser conversion. Unknown options are rejected before runtime.
- `pyoqp/oqp/library/set_basis.py` reads Basis Set Exchange electron shells and appends shell angular momentum/exponents/coefficients to the native OpenQP basis API.
- `source/types.F90` exposes `control_parameters` through the C ABI.
- `source/basis_api.F90::map_shell2basis_set` maps staged shells into `basis_set` arrays using Cartesian shell sizes (`NUM_CART_BF`).
- PySCF bridge/export code must not claim pure OpenQP basis behavior until the native basis representation carries pure shell metadata and all AO-order-dependent consumers are validated.

## Implementation path

1. Add parser/schema support for `[input] ispher` with validation limited to `-1`, `0`, and `1`.
2. Preserve Cartesian shell-size helpers and add explicit pure shell-size helpers for diagnostics/tests.
3. Log `ispher=0` as a Cartesian-equivalent compatibility mode with the SALC limitation stated.
4. Wire `control.ispher` into the native control structure and reject `ispher=1` in both Python setup and native basis mapping until true pure/spherical support exists.
5. Add a documented `ispher=1` input under `docs/dev/examples/` as documentation of the currently blocked pure-basis request; do not place it under the default `examples/` tree because `openqp --run_tests all` must not collect a known unsupported runtime path.

## Validation commands

- `python -m pytest -q tests/test_ispher_basis_option.py`
- `git diff --check`
- If runtime is rebuilt: run a small `openqp` smoke with `ispher=-1` and `ispher=0`, and verify `ispher=1` fails with the explicit pure/spherical NotImplemented/abort message.

## Claim boundaries / limitations

- OpenQP does **not** yet implement true `ispher=1` pure/spherical AO contractions or AO-order-dependent native integral support.
- OpenQP does **not** claim full GAMESS `ispher=0` SALC/population-analysis semantics; it accepts the keyword as Cartesian-equivalent compatibility mode and logs the limitation.
- No default behavior changes are intended; existing inputs without `ispher` continue through the Cartesian-compatible path.
