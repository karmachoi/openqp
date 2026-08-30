#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage: $0 RUNTIME_ROOT INPUT_ROOT RELATIVE_CALC RESULT_ROOT JOB_ID" >&2
  exit 64
fi

runtime_root=$1
input_root=$2
relative_calc=$3
result_root=$4
job_id=$5

scratch=$(mktemp -d "/tmp/oqp-uracil-${job_id}.XXXXXX")
result_dir="${result_root}/${relative_calc}"
mkdir -p "$result_dir"
cp -R "${input_root}/${relative_calc}/." "$scratch/"

export OMP_NUM_THREADS=7
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export PYTHONNOUSERSITE=1
export VECLIB_MAXIMUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
start_epoch=$(date +%s)
{
  echo "job_id=${job_id}"
  echo "host=$(hostname)"
  echo "pid=$$"
  echo "start_utc=${start_utc}"
  echo "runtime_root=${runtime_root}"
  echo "input=${relative_calc}/uracil.oqp"
  echo "input_sha256=$(shasum -a 256 "$scratch/uracil.oqp" | awk '{print $1}')"
  echo "omp_num_threads=${OMP_NUM_THREADS}"
  echo "veclib_maximum_threads=${VECLIB_MAXIMUM_THREADS}"
} > "$scratch/run-manifest.txt"

set +e
(
  cd "$scratch"
  "${runtime_root}/python/bin/openqp" --nompi --omp 7 uracil.oqp \
    > stdout.log 2> stderr.log
)
exit_code=$?
set -e

end_epoch=$(date +%s)
{
  echo "exit_code=${exit_code}"
  echo "end_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "wall_seconds=$((end_epoch - start_epoch))"
} >> "$scratch/run-manifest.txt"

cp -R "$scratch/." "$result_dir/"
if [[ $exit_code -eq 0 ]]; then
  touch "$result_dir/COMPLETE"
else
  touch "$result_dir/FAILED"
fi
exit "$exit_code"
