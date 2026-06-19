# oqp.quantum — Quantum-computing bridge

Export an OpenQP mean-field calculation as a **second-quantized molecular
Hamiltonian / FCIDUMP**, the universal entry point for quantum-computing
electronic-structure workflows (Qiskit Nature, OpenFermion, Block2/DMRG, …).

```python
from oqp.quantum import from_openqp

# `mol` is an OpenQP Molecule after a converged HF/DFT single point.
ham = from_openqp(mol, eri_ao=eri)        # eri = AO two-electron integrals
ham.to_fcidump("molecule.FCIDUMP")
```

## What it provides

| Piece | Status |
|-------|--------|
| One-electron MO integrals `h_pq` (from `OQP::Hcore` + MOs) | ✅ works today |
| Core / nuclear-repulsion energy, `nelec`, `MS2` metadata    | ✅ works today |
| AO→MO 1- and 2-index/4-index transforms                     | ✅ pure NumPy, tested |
| FCIDUMP write + read (8-fold symmetry, chemist notation)    | ✅ pure NumPy, tested |
| Two-electron MO integrals `(pq\|rs)`                         | ⏳ needs AO ERIs |

## The one remaining hook

OpenQP builds the two-electron repulsion integrals (ERIs) on the fly inside
the Fortran SCF Fock construction; they are **not yet exposed to Python**.
`from_openqp` therefore takes the AO ERIs explicitly via `eri_ao=` or a
`eri_provider=` callable.

The natural completion is a small CFFI getter — `oqp.int2e(mol)` populating an
`OQP::ERI_AO` tag, mirroring the existing `oqp.int1e` — after which
`from_openqp(mol)` yields a full FCIDUMP with no external integral source.

## Conventions

* Two-electron integrals use **chemist notation** `(pq|rs)` (FCIDUMP / PySCF
  convention), tensor index order `[p, q, r, s]`.
* Symmetric one-electron OpenQP matrices are stored packed-triangular and are
  unpacked with `unpack_triangular`.

## Modules

* `integrals.py` — `unpack_triangular`, `ao_to_mo_1body`, `ao_to_mo_2body`
  (no dependency on the compiled `oqp` extension).
* `fcidump.py` — `write_fcidump`, `read_fcidump`.
* `hamiltonian.py` — `MolecularHamiltonian`, `from_openqp`.

See `examples/QUANTUM/export_fcidump.py` and `tests/test_quantum_fcidump.py`.
