#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# OpenQP BLAS backend benchmark.
#
# Rebuilds OpenQP against each BLAS backend AVAILABLE on this machine and times
# a BLAS-bound MRSF-TDDFT single point, so the native-BLAS platform policy can
# be quantified from data. Policy: native BLAS per platform -- macOS (arm64 AND
# x86_64/Intel) -> Apple Accelerate; Linux x86_64 -> MKL; Linux aarch64 -> OpenBLAS.
# NetLib reference BLAS is the (slow) baseline everywhere for the speedup number.
#
# Run on each machine (Apple Silicon, zeus=Intel Mac, chc4=x86 Linux):
#     cd <openqp-repo> && ./benchmarks/blas/run_blas_bench.sh
#
# Env knobs: OMP_NUM_THREADS (default = physical cores), BENCH_REPEAT (default 2),
#            OQP_REPO (default = git toplevel of this script).
#
# NOTE: each backend does a full --force-reinstall rebuild (~1-3 min). Backends
# whose libraries are absent are skipped automatically. Accelerate is built with
# the LP64+8-byte-internal decoupling (LINALG_LIB_INT64=OFF + -fdefault-integer-8,
# OpenTrustRegion temporarily OFF) that the policy PR must make automatic.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="${OQP_REPO:-$(git -C "$HERE" rev-parse --show-toplevel 2>/dev/null || echo "$HERE/../..")}"
INPUT="$HERE/bench_mrsf.inp"
REPEAT="${BENCH_REPEAT:-2}"
OS="$(uname -s)"; ARCH="$(uname -m)"; HOSTN="$(uname -n | cut -d. -f1)"
CORES="$( (getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) )"
NTHREADS="${OMP_NUM_THREADS:-$CORES}"

pick() { for c in "$@"; do command -v "$c" >/dev/null 2>&1 && { echo "$c"; return; }; done; }
FC="$(pick gfortran-15 gfortran-14 gfortran)"
CC="$(pick gcc-15 gcc-14 gcc)"
CXX="$(pick g++-15 g++-14 g++)"
[ -z "$FC" ] && { echo "ERROR: no gfortran (need the GCC toolchain: brew install gcc@15 / apt install gfortran)"; exit 1; }

COMMON="-C cmake.define.USE_LIBINT=OFF -C cmake.define.ENABLE_OPENMP=ON -C cmake.define.ENABLE_MPI=OFF \
        -C cmake.define.CMAKE_C_COMPILER=$CC -C cmake.define.CMAKE_CXX_COMPILER=$CXX \
        -C cmake.define.CMAKE_Fortran_COMPILER=$FC"

have_mkl()  { [ -n "${MKLROOT:-}" ] || ldconfig -p 2>/dev/null | grep -qi 'libmkl_rt' || ls /opt/intel/oneapi/mkl/latest 2>/dev/null >/dev/null; }
have_ob64() { pkg-config --exists openblas64 2>/dev/null; }

# POLICY: native BLAS per platform -- macOS (arm64 AND x86_64/Intel) uses Apple
# Accelerate; Linux uses MKL (x86_64) / OpenBLAS (aarch64). MKL is therefore a
# Linux-only backend here. On Linux, if a oneAPI MKL is installed but its env was
# never sourced (MKLROOT unset), source it now: have_mkl() can see the install
# dir, but FindBLAS needs MKLROOT / the runtime env to actually locate and link
# MKL. On module-based clusters `module load imkl/...` sets MKLROOT, so this is
# skipped. NEVER source MKL on macOS -- it would make LINALG_LIB=auto pick the
# Apple-deprecated MKL (which also needs libiomp5 at runtime) over Accelerate.
if [ "$OS" != Darwin ] && [ -z "${MKLROOT:-}" ]; then
  for _mklenv in /opt/intel/oneapi/mkl/latest/env/vars.sh \
                 "${HOME}/intel/oneapi/mkl/latest/env/vars.sh"; do
    if [ -f "$_mklenv" ]; then
      set +u; . "$_mklenv" intel64 >/dev/null 2>&1 || true; set -u
      [ -n "${MKLROOT:-}" ] && break
    fi
  done
  if [ -z "${MKLROOT:-}" ]; then
    for _sv in /opt/intel/oneapi/setvars.sh "${HOME}/intel/oneapi/setvars.sh"; do
      if [ -f "$_sv" ]; then
        set +u; . "$_sv" >/dev/null 2>&1 || true; set -u
        [ -n "${MKLROOT:-}" ] && break
      fi
    done
  fi
  [ -n "${MKLROOT:-}" ] && echo "run_blas_bench: sourced oneAPI MKL env (MKLROOT=$MKLROOT)"
fi

# name|extra cmake defines  (reference NetLib is the baseline everywhere)
CONFIGS=( "netlib_ILP64|-C cmake.define.LINALG_LIB=netlib -C cmake.define.LINALG_LIB_INT64=ON" )
if [ "$OS" = Darwin ]; then
  CONFIGS+=( "accelerate_LP64|-C cmake.define.LINALG_LIB=auto -C cmake.define.LINALG_LIB_INT64=OFF -C cmake.define.ENABLE_OPENTRAH=OFF -C cmake.define.CMAKE_Fortran_FLAGS=-fdefault-integer-8" )
fi
# MKL is a Linux-only backend (macOS uses native Accelerate, selected via auto above).
[ "$OS" != Darwin ] && have_mkl && CONFIGS+=( "mkl_ILP64|-C cmake.define.LINALG_LIB=Intel10_64ilp -C cmake.define.LINALG_LIB_INT64=ON" )
have_ob64 && CONFIGS+=( "openblas64_ILP64|-C cmake.define.LINALG_LIB=OpenBLAS -C cmake.define.LINALG_LIB_INT64=ON" )

CSV="$HERE/results_${HOSTN}_${OS}-${ARCH}.csv"
echo "backend,build,linked_blas,wall_s,scf_energy_hartree,speedup_vs_netlib" > "$CSV"
echo "=================================================================="
echo " OpenQP BLAS benchmark   host=$HOSTN  $OS/$ARCH  cores=$CORES  threads=$NTHREADS"
echo " toolchain: $FC / $CC / $CXX   repeat=$REPEAT"
echo " backends: ${CONFIGS[*]%%|*}"
echo "=================================================================="
printf "%-20s %-8s %-14s %10s %16s %9s\n" backend build linked wall_s scf_energy vs_netlib

base_wall=""
for cfg in "${CONFIGS[@]}"; do
  name="${cfg%%|*}"; defs="${cfg#*|}"
  ( cd "$REPO" && eval pip install . --force-reinstall --no-deps $COMMON $defs ) >"/tmp/blasbench_build_$name.log" 2>&1
  if [ $? -ne 0 ]; then
    printf "%-20s %-8s\n" "$name" "FAIL"
    echo "$name,FAIL,,,," >> "$CSV"; continue
  fi
  # Resolve the built liboqp path. NOTE: `import oqp` prints a benign
  # "Failed to import mpi4py" line to stdout when mpi4py is absent; keep ONLY the
  # library path or it corrupts the linkage probe below and mislabels every
  # backend "static/netlib".
  LIB="$(python3 -c 'import oqp,glob,os;d=os.path.dirname(oqp.__file__);g=glob.glob(d+"/lib/liboqp*");print(g[0] if g else "")' 2>/dev/null | grep -E 'liboqp' | tail -1)"
  if [ "$OS" = Darwin ]; then
    linked="$(otool -L "$LIB" 2>/dev/null | grep -oiE 'Accelerate|mkl|openblas' | head -1)"
  else
    linked="$(ldd "$LIB" 2>/dev/null | grep -oiE 'mkl|openblas|lapack' | head -1)"
  fi
  # ILP64 MKL / bundled NetLib may not surface via ldd/otool; fall back to the
  # BLAS libs CMake recorded in the build log, else the configured backend name.
  [ -z "$linked" ] && linked="$(grep -hoiE 'libmkl|libopenblas|Accelerate|Reference BLAS|netlib' "/tmp/blasbench_build_$name.log" 2>/dev/null | grep -oiE 'mkl|openblas|Accelerate|netlib' | head -1)"
  [ -z "$linked" ] && linked="${name%%_*}"
  # OpenQP writes <input-basename>.log in the INPUT FILE'S directory, NOT the CWD.
  # Copy the input into an isolated per-backend dir and run THAT copy, so the log
  # (and scratch) stay isolated and land where we read them. Accept a timing ONLY
  # on a clean "PyOQP terminated" with a converged SCF energy -- never trust the
  # wall clock alone (a crash/abort must not be recorded as a result).
  RUNDIR="/tmp/blasbench_run_${name}"; rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"
  cp "$INPUT" "$RUNDIR/"; runin="$RUNDIR/$(basename "$INPUT")"
  runlog="$RUNDIR/$(basename "${INPUT%.inp}").log"
  best=""; energy=""; run_ok=1
  for i in $(seq 1 "$REPEAT"); do
    t0="$(python3 -c 'import time;print(time.time())')"
    ( cd "$RUNDIR" && OMP_NUM_THREADS="$NTHREADS" openqp "$runin" ) >"$RUNDIR/stdout_$i.log" 2>&1
    rc=$?
    t1="$(python3 -c 'import time;print(time.time())')"
    if [ "$rc" -ne 0 ] || ! grep -qi 'PyOQP terminated' "$runlog" 2>/dev/null \
                       || ! grep -qiE 'Final .* energy is' "$runlog" 2>/dev/null; then
      run_ok=0; break
    fi
    energy="$(grep -iE 'Final .* energy is' "$runlog" | tail -1 | grep -oE '\-?[0-9]+\.[0-9]+' | head -1)"
    w="$(python3 -c "print(round($t1-$t0,2))")"
    best="$(python3 -c "b='$best';print(min(float(b),$w) if b else $w)")"
  done
  if [ "$run_ok" -ne 1 ]; then
    cp -f "$runlog" "/tmp/blasbench_runfail_${name}.log" 2>/dev/null
    printf "%-20s %-8s %-14s %10s\n" "$name" "RUNFAIL" "$linked" "log:/tmp/blasbench_runfail_${name}.log"
    echo "$name,RUN_FAIL,$linked,,," >> "$CSV"; continue
  fi
  [ "$name" = "netlib_ILP64" ] && base_wall="$best"
  sp="$(python3 -c "bw='$base_wall';print(round(float(bw)/$best,2) if bw else 1.0)" 2>/dev/null)"
  printf "%-20s %-8s %-14s %10s %16s %9s\n" "$name" "ok" "$linked" "$best" "$energy" "${sp}x"
  echo "$name,ok,$linked,$best,$energy,${sp}" >> "$CSV"
done
echo "------------------------------------------------------------------"
echo "wall_s = best-of-$REPEAT total wallclock; speedup vs netlib baseline."
echo "scf_energy_hartree must agree across backends (numerical-correctness check)."
echo "results CSV -> $CSV"
echo "Attach the CSVs from all machines to the BLAS-policy PR."
