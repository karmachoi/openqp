#!/bin/bash
#SBATCH --job-name=oqp-nac-ec04414-verify
#SBATCH --partition=chc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830/build/logs/slurm-%j.out

set -euo pipefail

root=/bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830
attempt="$root/build/attempt-002"
venv="$attempt/venv"
source_root="$root/source/openqp-ec04414dc2226c67cd1b59355060faa2d5c0cd4b"
expected_wheel_sha=6ca50f554ee868aed591dee979d93ce443ef361c642e77a638cee47f16e8fab7
wheel=$(find "$attempt/wheelhouse" -maxdepth 1 -type f -name 'openqp-*.whl' -print -quit)

export MODULEPATH=/etc/scl/modulefiles:/etc/modulefiles:/usr/share/modulefiles:/opt/nvidia/hpc_sdk/modulefiles:/opt/soft/modulefiles/all
set +u
source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load GCC/12.3.0 Python/3.11.3-GCCcore-12.3.0 imkl/2023.1.0
set -u
export MKL_INTERFACE_LAYER=ILP64
export MKL_THREADING_LAYER=GNU
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OPENQP_ROOT=$("$venv/bin/python" -c \
  'import importlib.util; print(next(iter(importlib.util.find_spec("oqp").submodule_search_locations)))')

test "$(sha256sum "$wheel" | awk '{print $1}')" = "$expected_wheel_sha"
native_lib=$(find "$OPENQP_ROOT" -type f -name 'liboqp*.so' -print -quit)
test -n "$native_lib"
ldd "$native_lib" > "$attempt/ldd.txt"
! grep -Eq 'not found|libmkl_(intel|gf)_lp64' "$attempt/ldd.txt"
grep -Eq 'libmkl_(intel|gf)_ilp64|libmkl_rt' "$attempt/ldd.txt"

cd "$source_root"
"$venv/bin/python" -m pytest -q -rs \
  tests/test_namd_analytic_directional.py \
  tests/test_namd_baeck_an.py \
  tests/test_namd_diagnostics.py \
  tests/test_namd_rng_binding.py \
  tests/test_mrsf_nac_fortran_driver.py \
  tests/test_mrsf_nac_zvector_predictor.py \
  tests/test_mrsf_nac_zvector_residual_floor.py \
  tests/test_cphf_rohf_minres_contract.py \
  > "$attempt/regression-final.log" 2>&1
grep -Fq '58 passed, 1 skipped' "$attempt/regression-final.log"
grep -Fq 'SKIPPED [1] tests/test_namd_baeck_an.py:1638: compiled OpenQP runtime is not available' \
  "$attempt/regression-final.log"
test "$(grep -c '^SKIPPED ' "$attempt/regression-final.log")" -eq 1

{
  printf 'status=verified\nsource_commit=ec04414dc2226c67cd1b59355060faa2d5c0cd4b\n'
  printf 'wheel_sha256=%s\n' "$expected_wheel_sha"
  printf 'regression=58 passed; 1 intentional OpenMM-only force-field fingerprint skip\n'
  printf 'openmm=not installed; gas-phase Uracil NAMD does not use it\n'
  printf 'parent_build_job=5430249\nparent_finalizer_job=5430250\n'
  sha256sum "$wheel" "$native_lib" "$attempt/ldd.txt" "$attempt/regression-final.log"
  date -Ins
} > "$attempt/VERIFIED"

cat "$attempt/VERIFIED"
