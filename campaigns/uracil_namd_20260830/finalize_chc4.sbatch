#!/bin/bash
#SBATCH --job-name=oqp-nac-deb378d6-finalize
#SBATCH --partition=chc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=/bighome/cheolho.choi/openqp-uracil-namd-deb378d6-20260830/build/logs/slurm-%j.out

set -euo pipefail

campaign_root=/bighome/cheolho.choi/openqp-uracil-namd-deb378d6-20260830
source_commit=deb378d6e062495fddda0b4391f4d7821b7f6c1b
source_root="$campaign_root/source/openqp-$source_commit"
parent_root="$campaign_root/build/attempt-002"
attempt_root="$campaign_root/build/attempt-003"
expected_wheel_sha=221d9fa3d172575476be736e676258fea492888c5a71521dcb865f05cebe6bc4
gcc_root=/opt/soft/install/GCCcore/12.3.0
mkl_root=/opt/soft/install/imkl/2023.1.0/mkl/2023.1.0
python_root=/opt/soft/install/Python/3.11.3-GCCcore-12.3.0

test ! -e "$attempt_root"
mkdir -p "$attempt_root/wheelhouse"

export MODULEPATH=/etc/scl/modulefiles:/etc/modulefiles:/usr/share/modulefiles:/opt/nvidia/hpc_sdk/modulefiles:/opt/soft/modulefiles/all
set +u
source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load GCC/12.3.0 Python/3.11.3-GCCcore-12.3.0 imkl/2023.1.0
set -u

export PATH="$gcc_root/bin:$python_root/bin:$PATH"
export LD_LIBRARY_PATH="$gcc_root/lib64:$python_root/lib:$mkl_root/lib/intel64:${LD_LIBRARY_PATH:-}"
export MKLROOT="$mkl_root"
export MKL_INTERFACE_LAYER=ILP64
export MKL_THREADING_LAYER=GNU
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

wheel=$(find "$parent_root/wheelhouse" -maxdepth 1 -type f -name 'openqp-*.whl' -print -quit)
test -n "$wheel"
test "$(sha256sum "$wheel" | awk '{print $1}')" = "$expected_wheel_sha"
cp "$wheel" "$attempt_root/wheelhouse/"

"$python_root/bin/python3.11" -m venv "$attempt_root/venv"
"$attempt_root/venv/bin/python" -m pip install --upgrade pip
"$attempt_root/venv/bin/python" -m pip install 'numpy<2.2' scipy cffi basis_set_exchange geometric pytest
"$attempt_root/venv/bin/python" -m pip install --no-deps "$attempt_root/wheelhouse/$(basename "$wheel")"

package_root=$("$attempt_root/venv/bin/python" -c \
  'import importlib.util; print(next(iter(importlib.util.find_spec("oqp").submodule_search_locations)))')
native_lib=$(find "$package_root" -type f -name 'liboqp*.so' -print -quit)
test -n "$native_lib"
ldd "$native_lib" > "$attempt_root/ldd.txt"
if grep -Eq 'not found|libmkl_(intel|gf)_lp64' "$attempt_root/ldd.txt"; then
  echo 'FAIL: unresolved or LP64-linked library' >&2
  cat "$attempt_root/ldd.txt" >&2
  exit 1
fi
grep -Eq 'libmkl_(intel|gf)_ilp64|libmkl_rt' "$attempt_root/ldd.txt"

export OPENQP_ROOT="$package_root"
cd "$source_root"
"$attempt_root/venv/bin/python" -m pytest -q \
  tests/test_namd_analytic_directional.py \
  tests/test_namd_baeck_an.py \
  tests/test_namd_diagnostics.py \
  tests/test_namd_rng_binding.py \
  tests/test_mrsf_nac_fortran_driver.py \
  tests/test_mrsf_nac_zvector_predictor.py \
  tests/test_mrsf_nac_zvector_residual_floor.py \
  tests/test_cphf_rohf_minres_contract.py \
  > "$attempt_root/regression.log" 2>&1
if grep -Eq '[0-9]+ skipped' "$attempt_root/regression.log"; then
  echo 'FAIL: required regression was skipped' >&2
  cat "$attempt_root/regression.log" >&2
  exit 1
fi

{
  date -Ins
  hostname
  uname -a
  "$gcc_root/bin/gcc" --version | head -n 1
  "$gcc_root/bin/gfortran" --version | head -n 1
  "$attempt_root/venv/bin/python" -V
  printf 'source_commit=%s\n' "$(git -C "$source_root" rev-parse HEAD)"
  printf 'parent_attempt=attempt-002\nparent_job=5430247\n'
  sha256sum "$attempt_root/wheelhouse/$(basename "$wheel")" "$attempt_root/ldd.txt" "$attempt_root/regression.log"
} > "$attempt_root/VERIFIED"

echo "BUILD_VERIFIED $attempt_root"
