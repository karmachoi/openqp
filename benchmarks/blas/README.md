# OpenQP BLAS backend benchmark

Decides the native-BLAS platform policy from data instead of assumption — in
particular the open question **"on Intel Mac, is Accelerate or MKL faster?"**

## Why
OpenQP defaults to ILP64 BLAS. Apple Accelerate only exposes an LP64 BLAS to
Fortran, so `FindBLAS` silently falls back to the bundled **NetLib reference
BLAS** (unoptimized) on macOS and on any Linux without an ILP64 BLAS installed —
including the shipped PyPI wheels. This harness rebuilds against each available
backend and times a BLAS-bound MRSF-TDDFT single point so the policy PR can pick
the fastest native BLAS per platform and enforce it (no silent NetLib).

## Run it (on each machine)
```bash
cd <openqp-repo>
./benchmarks/blas/run_blas_bench.sh          # uses all physical cores
# or: OMP_NUM_THREADS=8 BENCH_REPEAT=3 ./benchmarks/blas/run_blas_bench.sh
```
Each backend does a full `pip install . --force-reinstall` (~1–3 min). Backends
whose libraries are missing are skipped. Prereqs: the **GCC toolchain**
(`brew install gcc@15` on macOS; `gfortran` on Linux). For MKL set `MKLROOT` (or
install oneAPI); for ILP64 OpenBLAS install `libopenblas64-dev` (Linux) so
`pkg-config --exists openblas64` passes.

## Machines to cover
| host | HW | backends the harness will try |
|------|----|-------------------------------|
| (this Apple Silicon mac) | macOS arm64 | netlib, **accelerate**, (openblas64) |
| **zeus** | macOS x86_64 (Intel) | netlib, **accelerate**, **mkl**, (openblas64) — *this answers the open question* |
| **chc4** | Linux x86_64 | netlib, **mkl**, **openblas64** |

## Backends
- `netlib_ILP64` — reference baseline (what a naive build ships). Slowest.
- `accelerate_LP64` — Apple Accelerate. Built with the **decoupling**:
  `LINALG_LIB_INT64=OFF` + `-fdefault-integer-8` (8-byte internals so the
  `int64_t` Python ABI stays valid) + `ENABLE_OPENTRAH=OFF`. The policy PR must
  make this automatic and re-enable OpenTrustRegion with correct integer widths.
- `mkl_ILP64` — Intel MKL (ILP64, native to the int64_t layer, no decoupling).
- `openblas64_ILP64` — ILP64 OpenBLAS.

## Output
Prints a table and writes `results_<host>_<os>-<arch>.csv` (wall_s, cpu_s,
speedup vs netlib, linked BLAS). **Attach all machines' CSVs to the PR** and set
the policy from the winners; on Intel Mac pick Accelerate only if it beats MKL
enough to justify losing MKL's ILP64 simplicity.
