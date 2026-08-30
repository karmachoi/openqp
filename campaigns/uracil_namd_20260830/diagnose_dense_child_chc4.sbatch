#!/bin/bash
#SBATCH --job-name=oqp-nac-ec04414-densechild
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

scratch=$(mktemp -d /tmp/oqp-dense-child-${SLURM_JOB_ID}.XXXXXX)
trap 'rm -rf "$scratch"' EXIT
"$venv/bin/python" - "$source_root/tests/test_namd_baeck_an.py" "$scratch/child.py" <<'PY'
import ast, pathlib, sys
tree = ast.parse(pathlib.Path(sys.argv[1]).read_text())
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name == 'test_dense_trajectory_is_appendable_and_restart_manifest_is_runnable':
        for stmt in node.body:
            if isinstance(stmt, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'script' for t in stmt.targets):
                pathlib.Path(sys.argv[2]).write_text(ast.literal_eval(stmt.value))
                raise SystemExit(0)
raise SystemExit('dense child script not found')
PY
export OQP_NAMD_TEST_ROOT="$scratch/runtime"
mkdir "$OQP_NAMD_TEST_ROOT"
cd "$source_root"
export PYTHONPATH="$source_root/pyoqp${PYTHONPATH:+:$PYTHONPATH}"
"$venv/bin/python" "$scratch/child.py"
