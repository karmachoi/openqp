"""Plot the genuine-geminal pTC-MRSF-CIS H2 curves from pes_h2_geminal.dat
(tc_h2_pes_geminal.F90). Left: bare MRSF-CIS (ROHF, dashed) / pTC-MRSF-CIS
(genuine cusp-fixed geminal, dotted) / FCI (solid) for S0,T1,S1. Right: the
dynamic-correlation recovered by the geminal for the singlet S0/S1 AND the
triplet T1 -- the key point, the geminal dresses both spin channels (the
spin-resolved cusp), unlike a closed-shell MP2 correlator.

matplotlib/numpy only, no pyscf.  Run: python3 plot_pes_geminal.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

d = np.loadtxt("pes_h2_geminal.dat")
R = d[:, 0]
fci, bare, ptc = d[:, 1:4], d[:, 4:7], d[:, 7:10]
states = [(0, 'black', r'$S_0\ (^1\Sigma_g^+)$'),
          (1, 'red',   r'$T_1\ (^3\Sigma_u^+)$'),
          (2, 'blue',  r'$S_1$ (open singlet)')]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for k, c, lab in states:
    ax1.plot(R, fci[:, k],  '-',  color=c, lw=1.2, label=lab)
    ax1.plot(R, ptc[:, k],  ':',  color=c, lw=1.4)
    ax1.plot(R, bare[:, k], '--', color=c, lw=1.0, alpha=0.9)
style = [Line2D([0], [0], color='k', ls='-',  label='FCI (6-311G)'),
         Line2D([0], [0], color='k', ls=':',  label='pTC-MRSF-CIS (geminal)'),
         Line2D([0], [0], color='k', ls='--', label='bare MRSF-CIS (ROHF)')]
l1 = ax1.legend(loc='upper right', fontsize=9, framealpha=0.95)
ax1.add_artist(l1)
ax1.legend(handles=style, loc='lower right', fontsize=9, framealpha=0.95)
ax1.set_xlabel('R (Bohr)'); ax1.set_ylabel('Energy (Hartree)')
ax1.set_title('H$_2$ dissociation (ROHF ref, genuine geminal pTC-MRSF-CIS)')
ax1.grid(alpha=0.3)

# absolute geminal correction (mHa) per state = how much the geminal lowers
# each state. Both the singlets AND the triplet get a nonzero correction --
# the spin-resolved cusp at work (a closed-shell MP2 correlator gives the
# triplet exactly zero).
for k, c, lab in states:
    corr = (bare[:, k] - ptc[:, k]) * 1000.0
    ax2.plot(R, corr, '-o', color=c, ms=2.5, lw=1.2, label=lab)
ax2.axhline(0, color='gray', lw=0.8)
ax2.set_xlabel('R (Bohr)')
ax2.set_ylabel('pTC correction  $E_{bare}-E_{pTC}$  (mHa)')
ax2.set_title('Both spin channels dressed (singlet AND triplet)')
ax2.legend(fontsize=9, framealpha=0.95)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("pes_h2_geminal.png", dpi=300)
print("wrote pes_h2_geminal.png")
