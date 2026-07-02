#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# OpenQP BLAS backend benchmark.
#
# Rebuilds OpenQP against each BLAS backend AVAILABLE on this machine and times
# a BLAS-bound MRSF-TDDFT single point, so the native-BLAS platform policy can
# be decided from data (esp. the open "Intel Mac: Accelerate vs MKL?" question).
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

# name|extra cmake defines  (reference NetLib is the baseline everywhere)
CONFIGS=( "netlib_ILP64|-C cmake.define.LINALG_LIB=netlib -C cmake.define.LINALG_LIB_INT64=ON" )
if [ "$OS" = Darwin ]; then
  CONFIGS+=( "accelerate_LP64|-C cmake.define.LINALG_LIB=auto -C cmake.define.LINALG_LIB_INT64=OFF -C cmake.define.ENABLE_OPENTRAH=OFF -C cmake.define.CMAKE_Fortran_FLAGS=-fdefault-integer-8" )
fi
have_mkl  && CONFIGS+=( "mkl_ILP64|-C cmake.define.LINALG_LIB=Intel10_64ilp -C cmake.define.LINALG_LIB_INT64=ON" )
have_ob64 && CONFIGS+=( "openblas64_ILP64|-C cmake.define.LINALG_LIB=OpenBLAS -C cmake.define.LINALG_LIB_INT64=ON" )

CSV="$HERE/results_${HOSTN}_${OS}-${ARCH}.csv"
echo "backend,build,linked_blas,wall_s,cpu_s,speedup_vs_netlib" > "$CSV"
echo "=================================================================="
echo " OpenQP BLAS benchmark   host=$HOSTN  $OS/$ARCH  cores=$CORES  threads=$NTHREADS"
echo " toolchain: $FC / $CC / $CXX   repeat=$REPEAT"
echo " backends: ${CONFIGS[*]%%|*}"
echo "=================================================================="
printf "%-20s %-8s %-14s %10s %10s %9s\n" backend build linked wall_s cpu_s vs_netlib

base_wall=""
for cfg in "${CONFIGS[@]}"; do
  name="${cfg%%|*}"; defs="${cfg#*|}"
  ( cd "$REPO" && eval pip install . --force-reinstall --no-deps $COMMON $defs ) >"/tmp/blasbench_build_$name.log" 2>&1
  if [ $? -ne 0 ]; then
    printf "%-20s %-8s\n" "$name" "FAIL"
    echo "$name,FAIL,,,," >> "$CSV"; continue
  fi
  LIB="$(python3 -c 'import oqp,glob,os;d=os.path.dirname(oqp.__file__);g=glob.glob(d+"/lib/liboqp*");print(g[0] if g else "")' 2>/dev/null)"
  if [ "$OS" = Darwin ]; then
    linked="$(otool -L "$LIB" 2>/dev/null | grep -oiE 'Accelerate|mkl|openblas' | head -1)"
  else
    linked="$(ldd "$LIB" 2>/dev/null | grep -oiE 'mkl|openblas|lapack' | head -1)"
  fi
  [ -z "$linked" ] && linked="static/netlib"
  best=""
  for i in $(seq 1 "$REPEAT"); do
    t0="$(python3 -c 'import time;print(time.time())')"
    OMP_NUM_THREADS="$NTHREADS" openqp "$INPUT" >"/tmp/blasbench_run_$name.log" 2>&1
    t1="$(python3 -c 'import time;print(time.time())')"
    w="$(python3 -c "print(round($t1-$t0,2))")"
    best="$(python3 -c "b='$best';print(min(float(b),$w) if b else $w)")"
  done
  cpu="$(grep -i 'Total CPU time' "/tmp/blasbench_run_$name.log" | tail -1 | grep -oE '[0-9]+\.[0-9]+' | head -1)"
  [ "$name" = "netlib_ILP64" ] && base_wall="$best"
  sp="$(python3 -c "bw='$base_wall';print(round(float(bw)/$best,2) if bw else 1.0)" 2>/dev/null)"
  printf "%-20s %-8s %-14s %10s %10s %9s\n" "$name" "ok" "$linked" "$best" "${cpu:-?}" "${sp}x"
  echo "$name,ok,$linked,$best,${cpu:-},${sp}" >> "$CSV"
done
echo "------------------------------------------------------------------"
echo "wall_s = best-of-$REPEAT total wallclock; speedup vs netlib baseline."
echo "results CSV -> $CSV"
echo "Attach the CSVs from all machines to the BLAS-policy PR."
