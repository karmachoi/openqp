"""Plot the NATIVE (pyscf-free) H2 dissociation curves from pes_h2_native.dat
produced by tc_h2_pes_native.F90. Three methods per state -- bare MRSF-CIS
(dashed), pTC-MRSF-CIS (dotted), FCI (solid) -- so the dynamic-correlation
recovery is visible (unlike a FCI-vs-pTC-only plot where the two coincide).

Reads only the native data file; matplotlib/numpy only, no pyscf.
Run:  python3 plot_pes_native.py   ->  pes_h2_native.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

d = np.loadtxt("pes_h2_native.dat")
R = d[:, 0]
# columns: R | S0 T1 S1 (fci) | S0 T1 S1 (bare) | S0 T1 S1 (ptc)
fci  = d[:, 1:4]
bare = d[:, 4:7]
ptc  = d[:, 7:10]

states = [(0, 'black', r'$S_0\ (^1\Sigma_g^+)$'),
          (1, 'red',   r'$T_1\ (^3\Sigma_u^+)$'),
          (2, 'blue',  r'$S_1$ (open-shell singlet)')]

plt.figure(figsize=(8, 6))
for k, c, lab in states:
    plt.plot(R, fci[:, k],  '-',  color=c, lw=2.4, label=lab)
    plt.plot(R, ptc[:, k],  ':',  color=c, lw=2.4)
    plt.plot(R, bare[:, k], '--', color=c, lw=1.6, alpha=0.9)

style = [Line2D([0], [0], color='k', ls='-',  label='FCI (6-311G)'),
         Line2D([0], [0], color='k', ls=':',  label='pTC-MRSF-CIS'),
         Line2D([0], [0], color='k', ls='--', label='bare MRSF-CIS (2,2)')]
leg1 = plt.legend(loc='upper right', framealpha=0.95, fontsize=9)
plt.gca().add_artist(leg1)
plt.legend(handles=style, loc='lower right', framealpha=0.95, fontsize=9)

plt.xlabel('R (Bohr)')
plt.ylabel('Energy (Hartree)')
plt.title(r'H$_2$ dissociation, native pTC-MRSF-CIS / 6-311G, pyscf-free')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("pes_h2_native.png", dpi=300)
print("wrote pes_h2_native.png")

# quick numeric summary: S0 dynamic-correlation recovery vs R
print("\n S0 recovery (pTC closes bare->FCI gap):")
for i in (0, 2, 5, 10, len(R) - 1):
    gap = fci[i, 0] - bare[i, 0]
    rec = (ptc[i, 0] - bare[i, 0]) / gap * 100 if abs(gap) > 1e-9 else 0.0
    print(f"  R={R[i]:.2f}  bare={bare[i,0]:.5f}  pTC={ptc[i,0]:.5f}  "
          f"FCI={fci[i,0]:.5f}  recovered={rec:4.0f}%")
