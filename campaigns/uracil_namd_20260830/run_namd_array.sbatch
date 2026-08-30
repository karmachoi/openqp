#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -euo pipefail

campaign_root=${OQP_URACIL_CAMPAIGN_ROOT:?set OQP_URACIL_CAMPAIGN_ROOT}
phase=${OQP_URACIL_PHASE:?set OQP_URACIL_PHASE}
array_index=${SLURM_ARRAY_TASK_ID:?run as a Slurm array}
runtime_root="$campaign_root/build/attempt-002/venv"
source_root="$campaign_root/source/openqp-ec04414dc2226c67cd1b59355060faa2d5c0cd4b"
jobs_file="$campaign_root/inputs/$phase-jobs.tsv"
gcc_root=/opt/soft/install/GCCcore/12.3.0
mkl_root=/opt/soft/install/imkl/2023.1.0/mkl/2023.1.0

export MODULEPATH=/etc/scl/modulefiles:/etc/modulefiles:/usr/share/modulefiles:/opt/nvidia/hpc_sdk/modulefiles:/opt/soft/modulefiles/all
set +u
source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load GCC/12.3.0 Python/3.11.3-GCCcore-12.3.0 imkl/2023.1.0
set -u
export LD_LIBRARY_PATH="$gcc_root/lib64:$mkl_root/lib/intel64:${LD_LIBRARY_PATH:-}"

line=$(sed -n "$((array_index + 1))p" "$jobs_file")
test -n "$line"
IFS=$'\t' read -r listed_index identity method relative_root <<< "$line"
test "$listed_index" = "$array_index"
input_root="$campaign_root/inputs/$relative_root"
result_root="$campaign_root/results/$phase/$method/$(basename "$relative_root")"
test -f "$input_root/uracil.oqp"
test -f "$input_root/geometry.xyz"
test -f "$input_root/velocity.au"
test -f "$input_root/run.env"
test -x "$runtime_root/bin/openqp"
test -f "$campaign_root/build/attempt-002/VERIFIED"
test ! -e "$result_root"

scratch=$(mktemp -d "/tmp/oqp-uracil-${SLURM_JOB_ID}-${array_index}.XXXXXX")
result_stage="${result_root}.partial-${SLURM_JOB_ID}"
mkdir -p "$result_stage"

finalize() {
  code=$?
  trap - EXIT INT TERM HUP
  {
    printf 'exit_code=%s\n' "$code"
    date -Ins
  } >> "$scratch/run-manifest.txt"
  find "$scratch" -type f ! -name SHA256SUMS -print0 \
    | sort -z | xargs -0 sha256sum > "$scratch/SHA256SUMS"
  cp -a "$scratch/." "$result_stage/"
  if [ "$code" -eq 0 ]; then
    mv "$result_stage" "$result_root"
  else
    mv "$result_stage" "${result_root}.failed-${SLURM_JOB_ID}"
  fi
  rm -rf "$scratch"
  exit "$code"
}
trap finalize EXIT INT TERM HUP

cp -a "$input_root/." "$scratch/"
cd "$scratch"
set -a
source run.env
set +a

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export MKL_INTERFACE_LAYER=ILP64
export MKL_THREADING_LAYER=GNU
export OPENQP_ROOT=$(
  "$runtime_root/bin/python" -c \
    'import importlib.util; print(next(iter(importlib.util.find_spec("oqp").submodule_search_locations)))'
)

{
  printf 'identity=%s\nmethod=%s\nphase=%s\n' "$identity" "$method" "$phase"
  printf 'job_id=%s\narray_index=%s\n' "$SLURM_JOB_ID" "$array_index"
  printf 'host=%s\npartition=%s\n' "$(hostname)" "${SLURM_JOB_PARTITION:-unknown}"
  printf 'source_commit=%s\n' "$(git -C "$source_root" rev-parse HEAD)"
  printf 'runtime_python=%s\nopenqp_root=%s\n' "$runtime_root/bin/python" "$OPENQP_ROOT"
  "$runtime_root/bin/python" --version
  sha256sum uracil.oqp geometry.xyz velocity.au initial-condition.json run.env \
    "$runtime_root/bin/openqp" "$OPENQP_ROOT/lib/liboqp.so"
  date -Ins
} > run-manifest.txt

"$runtime_root/bin/openqp" --nompi --omp 1 uracil.oqp \
  > stdout.log 2> stderr.log
