# QMRSF benchmark results

## Table 1. Ground-state total energy (Hartree)

| system/basis | ref ROHF | CAS=DK | icPT2-EN | icPT2-Dyall | FCI | dyn.corr (Dyall) |
|---|---|---|---|---|---|---|
| H4/sto-3g | -1.519813 | -2.102608 | -2.102608 | -2.102608 | -2.102608 | 0.00000 |
| H4/6-31g | -1.706013 | -2.151743 | -2.187565 | -2.200845 | -2.181861 | -0.04910 |
| CBD/sto-3g | -151.492067 | -151.842968 | -151.844221 | -151.844322 | -- | -0.00135 |
| CBD/6-31g | -153.323540 | -153.632441 | -153.652072 | -153.657059 | -- | -0.02462 |

## Table 2. Lowest vertical SINGLET excitation energies S0->Sn (eV)

_Spin-matched: each S_n is the n-th <S^2>-labelled singlet. DK==CAS on HF integrals, including the spin labels (GATE 1)._

_DK-DFT(grid) = the genuine grid-derived collinear kernel (adiabatic f_xc) on a BHHLYP/ROKS reference, CSF spin-adapted (multiplicity-block projection -> spin-PURE states, <S^2> exact); bare DK==CAS on HF integrals._

| system/basis | state | CAS=DK | DK-DFT(grid) | icPT2-EN | icPT2-Dyall |
|---|---|---|---|---|---|
| H4/sto-3g | S0->S1 | 9.456 | -- | 9.456 | 9.456 |
| H4/sto-3g | S0->S2 | 13.219 | -- | 13.219 | 13.219 |
| H4/sto-3g | S0->S3 | 18.607 | -- | 18.607 | 18.607 |
| H4/6-31g | S0->S1 | 7.176 | -- | 7.545 | 7.729 |
| H4/6-31g | S0->S2 | 9.820 | -- | 8.750 | 8.955 |
| H4/6-31g | S0->S3 | 14.613 | -- | 13.236 | 13.731 |
| CBD/sto-3g | S0->S1 | 5.716 | 6.094 | 5.712 | 5.719 |
| CBD/sto-3g | S0->S2 | 7.706 | 7.639 | 7.645 | 7.660 |
| CBD/sto-3g | S0->S3 | 12.752 | 13.176 | 12.643 | 12.682 |
| CBD/6-31g | S0->S1 | 4.441 | 4.616 | 4.432 | 4.508 |
| CBD/6-31g | S0->S2 | 5.709 | 5.675 | 4.813 | 4.884 |
| CBD/6-31g | S0->S3 | 10.584 | 10.617 | 9.615 | 9.825 |

## Table 3. Spin-resolved: lowest triplet T1 vs lowest singlet S1 (eV)

_The naive 'state 1' (lowest root above ground) is the TRIPLET; comparing it to a singlet-only column overstates the DK--icPT2 gap. T1 < S1 in every row._

| system/basis | method | ground 2S+1 | S0->T1 (triplet) | S0->S1 (singlet) |
|---|---|---|---|---|
| H4/sto-3g | CAS=DK | 1 | 4.057 | 9.456 |
| H4/sto-3g | icPT2-EN | 1 | 4.057 | 9.456 |
| H4/sto-3g | icPT2-Dyall | 1 | 4.057 | 9.456 |
| H4/6-31g | CAS=DK | 1 | 3.087 | 7.176 |
| H4/6-31g | icPT2-EN | 1 | 3.345 | 7.545 |
| H4/6-31g | icPT2-Dyall | 1 | 3.407 | 7.729 |
| CBD/sto-3g | CAS=DK | 1 | 2.158 | 5.716 |
| CBD/sto-3g | icPT2-EN | 1 | 2.164 | 5.712 |
| CBD/sto-3g | icPT2-Dyall | 1 | 2.163 | 5.719 |
| CBD/6-31g | CAS=DK | 1 | 1.560 | 4.441 |
| CBD/6-31g | icPT2-EN | 1 | 1.610 | 4.432 |
| CBD/6-31g | icPT2-Dyall | 1 | 1.594 | 4.508 |
