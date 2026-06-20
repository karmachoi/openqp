#!/usr/bin/env python3
"""Dump a reference case for the standalone Fortran QMRSF-icPT2 downfold test.

Writes fortran/icpt2_downfold_ref.dat from the validated NumPy multistate
prototype (qmrsf_icpt2_multistate.py): the raw matrices the Fortran downfold
must consume (H_PP, H_QP, the EN Q-diagonal, and the per-root Dyall H0) plus
the reference dressed energies for EN and Dyall denominators. The Fortran
independently diagonalizes H_PP, contracts the couplings, builds the
des-Cloizeaux symmetric effective Hamiltonian, and must reproduce both
reference spectra.
"""
import os
import numpy as np
import qmrsf_icpt2_multistate as M

N, NELEC, THOP, DELTA, NPD = 6, 6, 1.0, 0.2, 4

case = M.build_case(N, NELEC, thop=THOP, delta=DELTA)
Hfull, Pidx, Qidx = case['Hfull'], case['Pidx'], case['Qidx']
HPP = Hfull[np.ix_(Pidx, Pidx)]
HQP = Hfull[np.ix_(Qidx, Pidx)]
Hqq_en = np.diag(Hfull)[Qidx]
eP, _ = np.linalg.eigh(HPP)
H0dy = np.column_stack([M.dyall_denoms(case, eP[k]) for k in range(NPD)])

edEN, _, _, _ = M.icpt2_multistate(case, nP=NPD, denom='EN')
edDy, _, _, _ = M.icpt2_multistate(case, nP=NPD, denom='Dyall')

dimP, dimQ = len(Pidx), len(Qidx)
here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "fortran", "icpt2_downfold_ref.dat")


def row(vals):
    return " ".join(f"{v:.16e}" for v in vals) + "\n"


with open(out, "w") as f:
    f.write(f"{dimP} {dimQ} {NPD}\n")
    for i in range(dimP):
        f.write(row(HPP[i, :]))
    for i in range(dimQ):
        f.write(row(HQP[i, :]))
    f.write(row(Hqq_en))
    for i in range(dimQ):
        f.write(row(H0dy[i, :]))
    f.write(row(edEN))
    f.write(row(edDy))

print(f"wrote {out}")
print(f"  dimP={dimP} dimQ={dimQ} nPdress={NPD}")
print(f"  ref dressed EN = {edEN}")
print(f"  ref dressed Dy = {edDy}")
