#!/usr/bin/env python3
"""Dump CAS(4,4) active integrals + reference Ms=0 spectrum for the Fortran backbone test.
Uses the validated NumPy proto on a 4-orbital / 4-electron PPP model (the full space IS CAS(4,4),
so no frozen-core fold is needed -> a clean apples-to-apples reference for the Fortran)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qmrsf_icpt2_ppp_proto import build_ppp, spinorb, gen_dets, build_H

h, eri, eps = build_ppp(4, thop=1.0)          # 4 MOs, h(4,4), eri(4,4,4,4) chemist (pq|rs)
H1, g, _ = spinorb(h, eri)
dets = gen_dets(4, 2, 2)                       # all 36 Ms=0 determinants (2a2b in 4 orbitals)
H = build_H(dets, H1, g)
ev = np.linalg.eigvalsh(H)

out = os.path.join(os.path.dirname(__file__), "qmrsf_cas_ref.dat")
with open(out, "w") as f:
    f.write("4\n")
    for p in range(4):
        f.write(" ".join("%.16e" % h[p, q] for q in range(4)) + "\n")
    for p in range(4):
        for q in range(4):
            for r in range(4):
                f.write(" ".join("%.16e" % eri[p, q, r, s] for s in range(4)) + "\n")
    f.write("%d\n" % len(ev))
    f.write(" ".join("%.16e" % x for x in ev) + "\n")
print("dumped %s : ndet=%d  FCI ground=%.10f  (max ev=%.6f)" % (out, len(ev), ev[0], ev[-1]))
