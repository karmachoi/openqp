# Real OpenQP MRSF-CIS reference for H2 (after the zero-closed-shell fix)

MRSF-CIS = MRSF-TDDFT at 100% HF exchange, no XC (`functional=` empty).
Requires the `mrsfmntoia` fix (commit 82a7214a) — H2's ROHF triplet has no
doubly-occupied core, which crashed the unpatched MRSF Davidson.

Run:  openqp h2_mrsfcis_singlet.inp   (and ..._triplet.inp)

## H2 / aug-cc-pVDZ / R=1.4 bohr  (ROHF triplet = -0.7061829217 Ha)
  S0 = -1.1150791148 Ha   (<S^2>=0, singlet GS)
  S1 = -0.6941268390 Ha   ->  S1-S0 = 11.455 eV
  T1 = -0.7704108741 Ha   ->  T1-S0 =  9.379 eV  (<S^2>=2)

## H2 / cc-pVTZ / R=1.4 bohr  (ROHF triplet = -0.7740464716 Ha)
  S0 = -1.1347323726 Ha
  S1 = -0.6370624815 Ha   ->  S1-S0 = 13.542 eV

For comparison (independent, PySCF): in-basis FCI S0(aug-cc-pVDZ) = -1.164608,
T1-S0 = 10.509 eV, S1-S0 = 12.654 eV. Bare MRSF-CIS (CIS-level, no correlation)
lies ~49 mEh above FCI, as expected.

tau=0 gate: running either input with env OQP_PTC_MRSF=1 routes the reduced solve
through the non-Hermitian solver (tc_nonsym_tda_eig) and reproduces S0 bit-for-bit
(-1.1150791148), confirming the non-Hermitian path. The genuine TC effective
integrals (tau != 0) are NOT yet injected.
