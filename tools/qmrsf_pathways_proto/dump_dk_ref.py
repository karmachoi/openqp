#!/usr/bin/env python3
"""Dump a reference case for the standalone Fortran QMRSF-DK dressed-kernel test.

Writes fortran/dk_ref.dat from the validated NumPy live prototype
(qmrsf_dk_live_proto.py): the raw model the Fortran dressed-kernel pole search
must consume

    A0       (Ns x Ns)  adiabatic single-spin-flip response block (symmetric),
    omega_d  (Nd,)      bare 0OS closed-shell double-spin-flip energies,
    V        (Ns x Nd)  single<->0OS couplings,

plus the reference EXACT spectrum (eigvalsh of the explicit augmented matrix
[[A0, V],[V^T, diag(omega_d)]]) and the reference DRESSED spectrum from the
NumPy secular pole search. The Fortran independently assembles the dressed
kernel g_xc(omega), forms the pole-cancelled secular function, root-finds all
Ns+Nd roots, and must reproduce both spectra to < 1e-9.

Mirrors dump_icpt2_downfold_ref.py (the established pattern for the standalone
Fortran drafts under fortran/).
"""
import os
import numpy as np
import qmrsf_dk_live_proto as DK

A0, omega_d, V = DK.build_live_model()
Ns, Nd = A0.shape[0], omega_d.size

exact = DK.exact_spectrum(A0, omega_d, V)
dressed, _ = DK.dressed_roots(A0, omega_d, V)

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "fortran", "dk_ref.dat")


def row(vals):
    return " ".join(f"{v:.16e}" for v in np.atleast_1d(vals)) + "\n"


with open(out, "w") as f:
    f.write(f"{Ns} {Nd}\n")
    for i in range(Ns):
        f.write(row(A0[i, :]))            # A0, row-major
    f.write(row(omega_d))                 # bare 0OS double energies
    for i in range(Ns):
        f.write(row(V[i, :]))             # couplings, row-major
    f.write(row(exact))                   # reference exact spectrum (Ns+Nd)
    f.write(row(dressed))                 # reference dressed spectrum (Ns+Nd)

print(f"wrote {out}")
print(f"  Ns={Ns} Nd={Nd}  (Ns+Nd={Ns + Nd} roots)")
print(f"  ref exact   = {np.array2string(exact, precision=6)}")
print(f"  ref dressed = {np.array2string(dressed, precision=6)}")
print(f"  max|dressed-exact| (NumPy) = {np.abs(dressed - exact).max():.3e}")
