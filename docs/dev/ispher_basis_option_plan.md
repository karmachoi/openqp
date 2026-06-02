# ISPHER basis option plan

## Semantics

OpenQP now recognizes a GAMESS-compatible `[input]` keyword named `ispher`:

- `ispher=-1`: Cartesian basis functions/SALC, matching the historical OpenQP basis path and the GAMESS default Cartesian-compatible variational space.
- `ispher=0`: accepted compatibility mode. GAMESS builds SALCs from spherical harmonics while retaining Cartesian contaminants; OpenQP currently has no separate SALC-only population-analysis distinction, so this is documented and logged as a Cartesian-equivalent compatibility mode.
- `ispher=1`: request the GAMESS-style pure/spherical variational/SALC space that removes Cartesian contaminants (commonly summarized as d 6->5, f 10->7, and g 15->9 in the active AO space). OpenQP now accepts the keyword and maps native shell sizes through pure/spherical bookkeeping, but full AO-order-dependent integral/SCF validation is still required before broad production claims.

The default remains `ispher=-1`, preserving current Cartesian-compatible behavior. `control.ispher == 1` now reaches `source/basis_api.F90`, where the native shell-size mapping uses `NUM_PURE_BF` for `basis%nbf` and per-shell `basis%naos`.

## Current basis architecture

- `pyoqp/oqp/molecule/oqpdata.py` owns the input schema and parser conversion. Unknown options are rejected before runtime.
- `pyoqp/oqp/library/set_basis.py` reads Basis Set Exchange electron shells and appends shell angular momentum/exponents/coefficients to the native OpenQP basis API.
- `source/types.F90` exposes `control_parameters` through the C ABI.
- `source/basis_api.F90::map_shell2basis_set` maps staged shells into `basis_set` arrays using Cartesian shell sizes (`NUM_CART_BF`) for `ispher=-1/0` and pure shell sizes (`NUM_PURE_BF`) for `ispher=1`.
- PySCF bridge/export code must not claim fully validated pure OpenQP production behavior until AO-order-dependent consumers are validated against oracles.

## Implementation path

1. Add parser/schema support for `[input] ispher` with validation limited to `-1`, `0`, and `1`.
2. Preserve Cartesian shell-size helpers and add explicit pure shell-size helpers for diagnostics/tests.
3. Log `ispher=0` as a Cartesian-equivalent compatibility mode with the SALC limitation stated.
4. Wire `control.ispher` into the native control structure and use native shell-size mapping for `ispher=1` (`NUM_PURE_BF`) without the previous Python/native abort guards.
5. Keep the documented `ispher=1` input under `docs/dev/examples/` until runtime oracle validation is complete; do not place it under the default `examples/` tree before it has a stable reference result.

## Validation commands

- `python3 -m unittest tests/test_ispher_basis_option.py -v`
- `git diff --check`
- If runtime is rebuilt: run small `openqp` smoke tests with `ispher=-1`, `ispher=0`, and `ispher=1`, then compare `ispher=1` AO dimensions/energies against a trusted pure/spherical oracle before promoting examples to the default test tree.

## Claim boundaries / limitations

- OpenQP has initial `ispher=1` parser, Python setup, and native shell-size mapping support; full pure/spherical AO-order/integral/SCF validation remains incomplete.
- OpenQP does **not** claim full GAMESS `ispher=0` SALC/population-analysis semantics; it accepts the keyword as Cartesian-equivalent compatibility mode and logs the limitation.
- No default behavior changes are intended; existing inputs without `ispher` continue through the Cartesian-compatible path.
