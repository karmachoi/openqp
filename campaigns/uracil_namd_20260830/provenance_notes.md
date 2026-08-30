# Uracil initial-condition provenance

The archived NX file contains 10,000 Wigner initial conditions at 300 K.  The
coordinates are in bohr and the velocities are in atomic units.  Input
generation converts only the coordinates to angstrom because OpenQP's text
geometry parser expects angstrom; the velocity values are copied without a
unit conversion.

The legacy population directory contains 100 trajectories that have 10,001
saved state labels (0--5000 fs at 0.5-fs spacing): 13 start on S1 and 87 start
on S2.  The published lifetime tables instead use a denominator of 88, and the
S2-only fitting table uses a denominator of 82.  Thus the reported ensemble
contains six initial-S1 and 82 initial-S2 trajectories.  No surviving file in
the archive lists those 88 trajectory identifiers directly.

An exact binary reconstruction of the 82-member S2 subset was attempted from
all 87 complete S2 state histories and every 0--200-fs value in
`s2_population_200.tbl`.  No subset satisfies all tabulated counts, including
after testing adjacent half-femtosecond index offsets.  The likely causes are
that the aggregate table and the surviving per-trajectory tables were produced
from different archive revisions or restart concatenations.  This mismatch is
retained as provenance; the 88-member list is not guessed.

The new matched comparison therefore draws 50 conditions reproducibly from
the 100 unambiguously complete trajectories: six S1 and 44 S2 conditions,
using NumPy PCG64 seed 20260830 without replacement.  Every NAMD method uses
the same conditions, initial velocities, OpenQP counter-RNG seed and stream.
The 500-fs duration resolves the early S2-to-S1 and fast S1-to-S0 channels, but
does not determine the reported 1.9-ps slow lifetime or the 5-ps product yield.

The four matched methods are:

- `baseline_npi_iso`: overlap/NPI time-derivative coupling and isotropic
  energy-conserving velocity rescaling;
- `full_analytic`: analytic time-derivative coupling and analytic-NAC-directed
  velocity rescaling at every nuclear step;
- `rlzt10`: the same analytic formulation, with the transported linear
  Z-vector replacing the exact response solve after two exact warm-up steps
  and an exact refresh every ten nuclear steps.
- `ht_nac`: OpenQP's own NPI/FSSH propagation first selects an uncommitted hop
  candidate; only then does the same OpenQP NAMD driver evaluate the resident
  exact analytic NAC and use its direction for energy-conserving velocity
  rescaling.  It performs no second propagation or random draw and uses no
  external dynamics or QM/MM engine.

All four inputs run with the same OpenQP commit so that the comparison changes
only the named NAMD option.  HT-NAC changes when the established analytic NAC
is evaluated, not the MRSF reference, electronic response equations, or the
OpenQP nuclear/electronic propagation algorithm.

The native inputs retain the legacy Uracil MRSF setup: a triplet ROHF reference
(`SCF MULT=3`), singlet MRSF target fold (`TDDFT MULT=1`), BH&HLYP/6-31G*, and
0.5-fs nuclear steps.  The legacy S2 generator used `IROOT=3, NSTATE=4`; its S1
counterpart used `IROOT=2, NSTATE=3`.  The matched comparison evaluates four
singlet roots for every trajectory and selects active root 3 (S2) or 2 (S1),
so the electronic state space remains identical across methods.

Routine trajectories use matched SCF, TD-Davidson, and Z-vector thresholds of
`1e-8`, the least tight thresholds admitted by the current analytic NAC driver.
This is a cost-conscious NAMD setting, not an assertion of universal NAC
convergence.  One representative initial condition is evaluated separately at
`1e-10`; the production setting is accepted only if the change in gap, force,
full h vector, and hopping probability is negligible relative to the 0.5-fs
time-discretization effect.
