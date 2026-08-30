#!/bin/bash
#SBATCH --job-name=oqp-nac-ec04414-skipdiag
#SBATCH --partition=chc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830/build/logs/slurm-%j.out

set -euo pipefail

root=/bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830
venv="$root/build/attempt-002/venv"
source_root="$root/source/openqp-ec04414dc2226c67cd1b59355060faa2d5c0cd4b"
gcc_root=/opt/soft/install/GCCcore/12.3.0
mkl_root=/opt/soft/install/imkl/2023.1.0/mkl/2023.1.0

export LD_LIBRARY_PATH="$gcc_root/lib64:$mkl_root/lib/intel64:${LD_LIBRARY_PATH:-}"
export MKL_INTERFACE_LAYER=ILP64
export MKL_THREADING_LAYER=GNU
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OPENQP_ROOT=$("$venv/bin/python" -c \
  'import importlib.util; print(next(iter(importlib.util.find_spec("oqp").submodule_search_locations)))')

cd "$source_root"
"$venv/bin/python" -m pytest -q -rs \
  tests/test_namd_analytic_directional.py \
  tests/test_namd_baeck_an.py \
  tests/test_namd_diagnostics.py \
  tests/test_namd_rng_binding.py \
  tests/test_mrsf_nac_fortran_driver.py \
  tests/test_mrsf_nac_zvector_predictor.py \
  tests/test_mrsf_nac_zvector_residual_floor.py \
  tests/test_cphf_rohf_minres_contract.py
