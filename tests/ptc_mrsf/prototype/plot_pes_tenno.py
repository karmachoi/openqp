#!/usr/bin/env python3
"""H2 dissociation with the genuine Ten-no transcorrelation through MRSF-CIS.
MRSF-CIS = ROHF high-spin (triplet) reference + mixed-reference spin-flip singles
(no CASSCF(2,2) -- the open-shell frontier of the ROHF reference plus spin-flips
into all virtuals). Left: S0 PES -- bare MRSF-CIS sits above FCI (missing dynamic
correlation); the genuine TC pulls it onto the in-basis FCI across the whole
curve. Right: T1-S0 and S1-S0 excitation energies. The FCI S0 curve is validated
against pyscf to the printed digits at every R. (cc-pVDZ, Slater geminal gamma=0.7.)"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("pes_tenno.dat")
R = d[:, 0]
fci  = d[:, 1:4]   # S0 T1 S1
bare = d[:, 4:7]
tc   = d[:, 7:10]

fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))

# ---- (a) S0 potential energy curve ----
ax[0].plot(R, fci[:, 0],  color="#222222", lw=1.1, label="FCI (in-basis)")
ax[0].plot(R, bare[:, 0], color="#c62828", lw=1.1, ls="--", label="bare MRSF-CIS (ROHF ref)")
ax[0].plot(R, tc[:, 0],   color="#1565c0", lw=1.1, ls="-.", label="TC-MRSF-CIS (Ten-no)")
ax[0].set_xlabel("H--H distance  R  (bohr)")
ax[0].set_ylabel(r"$S_0$ energy  (Hartree)")
ax[0].set_title("(a) Ground-state PES")
ax[0].legend(loc="lower right", fontsize=8, frameon=False)

# ---- (b) excitation energies ----
ev = 27.211386
ax[1].plot(R, (fci[:, 1]-fci[:, 0])*ev,  color="#222222", lw=1.1, label=r"FCI  $T_1$")
ax[1].plot(R, (tc[:, 1]-tc[:, 0])*ev,    color="#1565c0", lw=1.1, ls="-.", label=r"TC  $T_1$")
ax[1].plot(R, (bare[:, 1]-bare[:, 0])*ev,color="#c62828", lw=1.1, ls="--", label=r"bare $T_1$")
ax[1].plot(R, (fci[:, 2]-fci[:, 0])*ev,  color="#222222", lw=1.1, alpha=0.55, label=r"FCI  $S_1$")
ax[1].plot(R, (tc[:, 2]-tc[:, 0])*ev,    color="#1565c0", lw=1.1, ls="-.", alpha=0.55, label=r"TC  $S_1$")
ax[1].set_xlabel("H--H distance  R  (bohr)")
ax[1].set_ylabel("excitation energy  (eV)")
ax[1].set_title("(b) $T_1$--$S_0$ and $S_1$--$S_0$")
ax[1].legend(loc="upper right", fontsize=7.5, frameon=False, ncol=2)

for a in ax:
    a.tick_params(labelsize=8)
    a.axhline(0 if a is ax[1] else -1.0, color="0.8", lw=0.6, zorder=0)

fig.tight_layout()
fig.savefig("pes_tenno.png", dpi=300)
print("wrote pes_tenno.png")
