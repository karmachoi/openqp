#!/bin/bash
# QMRSF-icPT2 consolidated regression suite.
# Runs every icPT2 validation gate (live OpenQP runs + their pyscf-free oracles, plus the
# standalone Fortran algebra tests) and prints one PASS/FAIL summary. Run from this dir
# (tools/qmrsf_pathways_proto/stageB) after building liboqp in the worktree.
#
#   ./run_icpt2_tests.sh
#
# Requires: a built liboqp staged at $OPENQP_ROOT (default /tmp/qmrsf_root) and gfortran-15.
set -u
ROOT=/Users/cheolhochoi/Documents/openqp-private-qmrsf-pathways
export OPENQP_ROOT=${OPENQP_ROOT:-/tmp/qmrsf_root}
export PYTHONPATH=$ROOT/pyoqp
PYOQP="python3 $ROOT/pyoqp/oqp/pyoqp.py"
HERE=$(pwd); FORT=$ROOT/tools/qmrsf_pathways_proto/fortran
pass=0; fail=0
ok(){ echo "  [PASS] $1"; pass=$((pass+1)); }
no(){ echo "  [FAIL] $1"; fail=$((fail+1)); }

echo "================ QMRSF-icPT2 regression suite ================"

# --- live runs + pyscf-free oracle gates (each gate exits 0 on PASS) ---
echo "-- live: H4/STO-3G transform gate (CAS(4,4)=FCI) --"
rm -f qmrsf_icpt2_live.dat qmrsf_cact_live.dat
$PYOQP h4_quintet_icpt2.inp >/tmp/t_h4sto.out 2>&1
python3 route_a_oracle.py >/tmp/t_h4sto_gate.out 2>&1 && ok "H4/STO-3G transform vs closed-form oracle" || no "H4/STO-3G transform (see /tmp/t_h4sto_gate.out)"

echo "-- live: H4/6-31G external-Q downfold gate (EN+Dyall) --"
rm -f qmrsf_icpt2_full_live.dat
$PYOQP h4_quintet_icpt2_631g.inp >/tmp/t_h4631.out 2>&1
python3 gate_icpt2_full.py >/tmp/t_h4631_gate.out 2>&1 && ok "H4/6-31G downfold vs NumPy oracle" || no "H4/6-31G downfold (see /tmp/t_h4631_gate.out)"

echo "-- live: H6/STO-3G frozen-core gate (ncore=1) --"
rm -f qmrsf_icpt2_full_live.dat qmrsf_cfull_live.dat
$PYOQP h6_quintet_icpt2.inp >/tmp/t_h6.out 2>&1
python3 gate_frozencore.py >/tmp/t_h6_gate.out 2>&1 && ok "H6/STO-3G frozen-core vs closed-form oracle" || no "H6/STO-3G frozen-core (see /tmp/t_h6_gate.out)"

# --- standalone Fortran algebra tests (print 'RESULT: PASS') ---
echo "-- standalone Fortran algebra --"
( cd "$FORT" && gfortran-15 -O2 qmrsf_backbone_core.f90 -o /tmp/bb -framework Accelerate >/dev/null 2>&1 && /tmp/bb 2>&1 | grep -q "RESULT: PASS" ) && ok "backbone CAS(4,4) det-CI" || no "backbone det-CI"
( cd "$FORT" && gfortran-15 -O2 qmrsf_icpt2_downfold.f90 -o /tmp/df -framework Accelerate >/dev/null 2>&1 && /tmp/df 2>&1 | grep -q "RESULT: PASS" ) && ok "downfold kernel (EN+Dyall)" || no "downfold kernel"
( cd "$FORT" && gfortran-15 -O2 qmrsf_icpt2_full.f90 -o /tmp/fl -framework Accelerate >/dev/null 2>&1 && /tmp/fl 2>&1 | grep -q "RESULT: PASS" ) && ok "full pipeline (det-CI -> downfold)" || no "full pipeline"
( cd "$FORT" && cp -f "$HERE/qmrsf_icpt2_full_live.dat" . 2>/dev/null; gfortran-15 -O2 qmrsf_icpt2_contracted.f90 -o /tmp/ct -framework Accelerate >/dev/null 2>&1 && /tmp/ct 2>&1 | grep -q "RESULT: PASS" ) && ok "contracted engine (no FCI list)" || no "contracted engine"

echo "============================================================="
echo "  SUMMARY: $pass passed, $fail failed"
[ "$fail" -eq 0 ] && { echo "  RESULT: ALL PASS"; exit 0; } || { echo "  RESULT: FAILURES"; exit 1; }
