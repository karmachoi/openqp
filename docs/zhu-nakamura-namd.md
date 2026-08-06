# Zhu–Nakamura global-switching NAMD

OpenQP provides the Yu–Zhu multidimensional Zhu–Nakamura (ZN) global-switching
algorithm for localized, same-spin avoided crossings. Select it with
`hop_method=zhu_nakamura`; the default remains `hop_method=fssh`.

```text
mrsf(nstate=2)/bhhlyp/sto-3g
namd(T0,nstep=20,dt=0.1,velocity="velocity.dat",hop_method=zhu_nakamura)
geom="formaldehyde.xyz"
```

For each adjacent state pair, OpenQP detects a local minimum of the adiabatic
energy gap from three consecutive nuclear points. The endpoint adiabatic
gradients are cross-interpolated to estimate the two diabatic forces. The
native Fortran kernel then evaluates the Zhu–Nakamura effective coupling
`a2`, effective collision energy `b2`, and global switching probability. The
`thrshe` value is the largest centre-point gap, in Hartree, that may be treated
as a ZN crossing.

The probability follows the two branches of the published avoided-crossing
formula: the inner term is `b2*b2 + 1` when the mass-weighted diabatic-force
product is positive and `abs(b2*b2 - 1)` when it is negative (`b2` denotes the
paper's \(b^2\)). OpenQP also applies the published limits `a2 <= 0.001`
(adiabatic, zero switching probability) and `a2 >= 1000` (diabatic, unit
switching probability).

A centre-point event can only be decided after the following point has been
calculated. If a hop is accepted, OpenQP discards that trial point, restores
the centre geometry, rescales only the momentum components along the
self-consistent force-difference direction, and repeats velocity Verlet and
the electronic calculation on the new active surface. Thus the first possible
event is centre step 1 after the step-2 electronic result is available;
`first_hop_step=1` is appropriate.

ZN uses all-state gradients at every nuclear point and is therefore more
expensive per point than FSSH with only the active-state gradient. Wavefunction
overlaps are still evaluated for MO/root/phase tracking and NACME diagnostics,
but FSSH coefficient flux does not decide a ZN hop. Electronic coefficients in
the trajectory are diagnostic in this mode.

The packed trajectory stores `hop_probabilities`, `zn_event_step`, `zn_a2`,
and `zn_b2`. A restart checkpoint also stores the last two structures,
adiabatic energies, and all-state gradients so the first post-restart decision
uses exactly the same three-point history.

Current scope: fixed-step, same-spin gas-phase/ODP dynamics with localized
avoided crossings. SOC crossing/parallel formulas, nonadiabatic-tunnelling
surfaces, and QM/MM rollback are intentionally rejected rather than silently
approximated.

Method reference: L. Yu, C. Xu, Y. Lei, C. Zhu, and Z. Wen, *Phys. Chem. Chem.
Phys.* **16**, 25883–25895 (2014), DOI: 10.1039/C4CP03498H.
