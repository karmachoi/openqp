# DFT XC-grid "compute-once / reuse" optimizations

Two independent, **default-OFF**, env-gated optimizations of the SCF DFT cost.
Both leave the converged SCF energy/gradient unchanged (Opt 1 bit-exact; Opt 2
within the SCF convergence tolerance). OpenQP already reuses the J/K (HF) part of
the Fock matrix incrementally; these target the XC build, which did not reuse
anything across SCF iterations.

| Env var | Default | Effect |
|---|---|---|
| `OQP_XC_PHI_CACHE=1` | off | **Opt 1** — cache the collocation matrix Φ across SCF iterations |
| `OQP_XC_INCDFT=1` | off | **Opt 2** — incremental XC reuse (experimental; *measured non-improvement*, see below) |
| `OQP_XC_TIMING=1` | off | per-iteration timers: `[SCFTIME]` (J/K vs XC wall split) and `[XCTIME]` (geom-Φ vs density-driven) |

All measurements: Apple-silicon macOS, gfortran-15, Release, 4 OpenMP threads,
6-31G(d), Becke grid. The box was heavily loaded during runs, so **per-iteration
`[SCFTIME]` ratios and SCF iteration counts are the robust metrics**; total-wall
numbers are min-of-3 but still noisy.

---

## Opt 1 — Collocation-Φ cache (SHIPPED, exact)

`Φ[μ,i] = φ_μ(r_i)` (basis values + grid derivatives) depends only on
geometry + basis + grid + integration threshold, **not** on the density, yet the
native grid loop (`run_xc` → `compAOs`/`pruneAOs`) recomputed it on every Fock
build. Opt 1 stores the post-pruning Φ block (plus weights and significant-AO
metadata) per grid slice on the first build of a geometry and replays it on
subsequent SCF iterations.

* New module `source/dftlib/dft_gridint_phi_cache.F90`; hooks in
  `run_xc` (`dft_gridint.F90`); opt-in set only by the repeated SCF Fock build
  (`dmatd_blk`), so one-shot consumers (gradients, response) never cache.
* Validity keyed by a geometry hash + derivative order + DFT threshold → any
  change rebuilds transparently. The replayed Φ is bit-for-bit identical to the
  recomputed Φ.

**Exactness** (converged energy, `OQP_XC_PHI_CACHE` 1 vs 0):

| system / functional | ΔE (Ha) | iters off → on |
|---|---|---|
| benzene / B3LYPV5 | 0.0e+00 | 9 → 9 |
| benzene / M06-2X (meta-GGA) | 0.0e+00 | 9 → 9 |
| C20H42 / B3LYPV5 | 1e-10 | 11 → 11 |
| C20H42 / PBE | 0.0e+00 | 12 → 12 |
| C20H42 / PBE, fine grid (96×302, unpruned) | 0.0e+00 | 12 → 12 |

(The 1e-10 is floating-point summation-order noise from dynamic OpenMP
scheduling, present in the baseline too — not a cache error.)

**Speedup.** The win scales with the XC fraction of the Fock build, which is set
by grid density and functional. The cache removes ~85–90 % of the collocation
(`compAOs`) cost; what that buys in total wall depends on how big XC is:

| case | XC build / iter (off → replay) | XC frac of Fock | total wall (off → on) |
|---|---|---|---|
| C20H42 B3LYPV5, SG1 grid | 0.263 → 0.185 s (−30 %) | ~13 % | 20.5 → 19.8 s (~3 %) |
| C20H42 PBE, **fine** grid | 1.26 → 0.41 s (−67 %) | ~41–59 % | 34.7 → 25.7 s (**−26 %**) |

The geom-Φ phase per iteration (aggregate thread-seconds, `[XCTIME]`) drops on
replay, confirming `compAOs` is eliminated, e.g. C20H42 PBE SG1: 0.45 → 0.038 s.

**Cost / when it helps.**
* Memory: the cache is `Σ_slices numAOs_pruned × numPts × numAOVecs × 8 B`
  (benzene SG1 ≈ 94 MB; C20H42 SG1 ≈ 1.7 GB; fine grids larger). Logged as
  `cacheMB` under `OQP_XC_TIMING`. Off by default precisely because of this
  memory/recompute trade.
* For **hybrid** functionals on small (SG1) grids the 2e J/K build dominates the
  Fock matrix (~87 %), so XC — and therefore this optimization — is a small share
  of total SCF time. The win is large for **dense/fine grids, meta-GGA, and pure
  functionals**, and it grows further when the sibling 2e-integral screening work
  (PR #238) shrinks the J/K cost (XC becomes a bigger fraction). The features
  compose.

---

## Opt 2 — Incremental XC (IncDFT): experimental, **measured non-improvement**

Implemented as whole-matrix XC reuse: store `V_xc[D_ref]`/`E_xc[D_ref]` from the
last full build and reuse it on iterations inside a late-SCF DIIS-error window,
with a periodic forced full rebuild and a return to full builds near convergence
so the fixed point stays exact. Gated `OQP_XC_INCDFT=1`
(`source/dftlib/dft_incdft.F90`; plumbed `scf.F90` → `calc_fock` → `calc_jk_xc`).

**Result: it does not help — it hurts convergence and is a net loss.** The
converged energy stays within the 1e-7 gate (the forced full builds near
convergence preserve correctness), but the iteration count blows up:

| system | baseline iters | IncDFT iters | net wall |
|---|---|---|---|
| benzene / PBE (default window) | 9 | 17 | slower |
| C20H42 / PBE, SG1 (default window) | 12 | 31 | 25.4 → **47.6 s** |

A window sweep on benzene/PBE shows the trade is fundamental, not a tuning miss:

```
start  stop   refresh   iters  reuses  ΔE
3e-2   1e-4   5          17      6      -5.9e-9
5e-3   5e-4   2          13      2       3.8e-9
1e-3   1e-4   1          10      1       4.0e-10
1e-3   5e-4   2           9      0       0          (window so narrow nothing reused)
```

**Each reused XC build costs ≈ one extra SCF iteration.** Because XC is nonlinear
and the reused matrix is *global*, freezing it produces a Fock inconsistent with
the still-evolving density, which DIIS must work off. And one extra SCF iteration
(which includes a J/K build — always ≥ the XC build here) costs more than the XC
build it skipped. So in this codebase whole-matrix XC reuse can only lose.

**The correct approach is per-batch screening** (Q-Chem-style incremental DFT):
keep `V_xc[D]` *accurate* by recomputing only the grid batches where `ΔP` is
significant and reusing the rest, via reference subtraction
`V = V_ref + Σ_active (V_s[D] − V_s[D_ref])`. Because the assembled Fock stays
correct to the screening threshold, convergence is **not** penalized, and the
saving is the screened-out batch fraction (grows late in SCF as `ΔP` sparsifies).
This needs the engine to evaluate a reference density per active batch (cheap
given the Opt-1 Φ cache) and per-slice `ΔP` screening — a larger change left as
future work. The current `OQP_XC_INCDFT` path is retained, off by default and
clearly experimental, as the gating/reference-store scaffolding for that work and
to make the negative result reproducible (`[SCFTIME] ... xc_reused=T/F`).

## Reproduce

```sh
# build (note LP64 BLAS flag, required on this Mac):
cmake -G Ninja -B build -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 \
  -DCMAKE_Fortran_COMPILER=gfortran-15 -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON -DENABLE_OPENMP=ON -DENABLE_PYTHON=OFF -DUSE_LIBINT=OFF \
  -DLINALG_LIB=auto -DLINALG_LIB_INT64=OFF
ninja -C build oqp

# run with timing, toggling the optimizations:
OQP_XC_PHI_CACHE=1 OQP_XC_TIMING=1  pyoqp input.inp     # Opt 1
OQP_XC_INCDFT=1    OQP_XC_TIMING=1  pyoqp input.inp     # Opt 2 (experimental)
```
