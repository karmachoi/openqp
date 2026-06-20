#!/usr/bin/env python3
"""Dump the P/Q blocks + reference icPT2 result for the standalone Fortran downfold test.
Uses the validated NumPy proto on the 6-orbital PPP model with a CAS(4,4) active window
(core MO1 + active MO2-5 + virtual MO6), so there is a genuine external Q space to downfold."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qmrsf_icpt2_ppp_proto import run_case, icpt2

_, _, _, dets, Hfull, Pidx, part = run_case(6, 6, 1.0)
ECAS, Eic, sig = icpt2(Hfull, dets, Pidx, nroots=1)[0]      # state-specific EN downfold (ground)
P = list(Pidx); Q = [i for i in range(len(dets)) if i not in set(P)]
HPP = Hfull[np.ix_(P, P)]; HPQ = Hfull[np.ix_(P, Q)]; Hqq = np.diag(Hfull)[Q]

out = os.path.join(os.path.dirname(__file__), "qmrsf_icpt2_ref.dat")
with open(out, "w") as f:
    f.write("%d %d\n" % (len(P), len(Q)))
    for i in range(len(P)): f.write(" ".join("%.16e" % HPP[i, j] for j in range(len(P))) + "\n")
    for i in range(len(P)): f.write(" ".join("%.16e" % HPQ[i, j] for j in range(len(Q))) + "\n")
    f.write(" ".join("%.16e" % x for x in Hqq) + "\n")
    f.write("%.16e %.16e %.16e\n" % (ECAS, Eic, sig))
print("dumped nP=%d nQ=%d  E_CAS=%.8f  icPT2=%.8f  sigma=%.8f" % (len(P), len(Q), ECAS, Eic, sig))
