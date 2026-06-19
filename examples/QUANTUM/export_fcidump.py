#!/usr/bin/env python3
"""Export an OpenQP mean-field calculation as a FCIDUMP for quantum computing.

FCIDUMP is the standard hand-off point from a classical quantum-chemistry
mean field to a quantum-computing electronic-structure workflow. The file this
script writes can be loaded by, e.g.:

* Qiskit Nature  -- ``qiskit_nature.second_q.formats.fcidump.FCIDump.from_file``
* OpenFermion    -- via PySCF's ``MolecularData`` / FCIDUMP loaders
* Block2 / DMRG, Dice/SHCI, and most active-space front ends

Usage
-----
    openqp h2.inp            # run an HF/DFT single point first
    python export_fcidump.py h2.inp h2.FCIDUMP

Two-electron integrals
----------------------
OpenQP builds the two-electron repulsion integrals (ERIs) inside the Fortran
SCF and does not yet expose them to Python. Until the ``oqp.int2e`` getter is
added, supply the AO ERIs through ``eri_provider``. This example shows the
wiring with a placeholder provider; replace it with your integral source (or
OpenQP's native ERIs once exposed).
"""

import sys

from oqp.pyoqp import Runner
from oqp.quantum import from_openqp


def eri_provider(mol):
    """Return AO two-electron integrals (nao,nao,nao,nao) in chemist notation.

    Placeholder. Replace with OpenQP's native ERIs once ``oqp.int2e`` lands,
    or with an external engine that uses the SAME basis and geometry.
    """
    raise NotImplementedError(
        "Provide AO two-electron integrals here. The one-electron Hamiltonian "
        "and all metadata are already available from OpenQP; only the ERIs "
        "are pending the oqp.int2e hook (see oqp/quantum/hamiltonian.py).")


def main():
    if len(sys.argv) < 3:
        sys.exit(f"usage: {sys.argv[0]} <input.inp> <out.FCIDUMP>")
    input_file, out = sys.argv[1], sys.argv[2]

    runner = Runner(project=input_file.rsplit(".", 1)[0], input_file=input_file)
    runner.run()
    mol = runner.mol

    # One-electron-only Hamiltonian works today (no ERIs needed):
    ham = from_openqp(mol)
    print(f"norb = {ham.n_orbitals}, nelec = {ham.n_electrons}, "
          f"ms2 = {ham.ms2}, E_core = {ham.core_energy:.10f}")
    print("one-body h_pq (MO basis):")
    print(ham.one_body)

    # Full FCIDUMP requires two-electron integrals:
    try:
        ham = from_openqp(mol, eri_provider=eri_provider)
        ham.to_fcidump(out)
        print(f"wrote {out}")
    except NotImplementedError as exc:
        print(f"\n[two-electron integrals not supplied] {exc}")


if __name__ == "__main__":
    main()
