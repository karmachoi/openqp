#!/usr/bin/env python3
"""Dump MO integrals + reference for the standalone Fortran full-pipeline icPT2 test.

Writes fortran/icpt2_full_ref.dat: the spatial MO integrals (h_mo, chemist eri_mo)
plus (norb, na, nb, ncore, nact) and the reference EN dressed energies. The Fortran
qmrsf_icpt2_full.f90 then builds the spin-orbital tensors, enumerates the full
determinant space, assembles H by Slater-Condon, partitions P (CAS window) / Q,
and runs the multistate EN downfold -- reproducing this reference from integrals.
"""
import os
import numpy as np
import qmrsf_icpt2_multistate as M

N, NELEC, THOP, DELTA, NPD = 6, 6, 1.0, 0.2, 4
NCORE = (N - 4) // 2          # = 1
NACT = 4
NA = NB = NELEC // 2          # = 3

case = M.build_case(N, NELEC, thop=THOP, delta=DELTA)
h_mo, eri_mo = case['h_mo'], case['eri_mo']
edEN, _, eP, _ = M.icpt2_multistate(case, nP=NPD, denom='EN')

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "fortran", "icpt2_full_ref.dat")
with open(out, "w") as f:
    f.write(f"{N} {NA} {NB} {NCORE} {NACT} {NPD}\n")
    for p in range(N):
        f.write(" ".join(f"{h_mo[p, q]:.16e}" for q in range(N)) + "\n")
    for p in range(N):
        for q in range(N):
            for r in range(N):
                f.write(" ".join(f"{eri_mo[p, q, r, s]:.16e}" for s in range(N)) + "\n")
    f.write(" ".join(f"{x:.16e}" for x in edEN) + "\n")

print(f"wrote {out}")
print(f"  norb={N} na={NA} nb={NB} ncore={NCORE} nact={NACT} nPdress={NPD}")
print(f"  CAS eP[:4]  = {eP}")
print(f"  ref EN      = {edEN}")
