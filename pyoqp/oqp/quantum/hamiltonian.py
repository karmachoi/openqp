"""Build the second-quantized molecular Hamiltonian from an OpenQP molecule.

This is the bridge that turns a converged OpenQP mean-field calculation into
the input expected by quantum-computing electronic-structure workflows
(Qiskit Nature, OpenFermion, PennyLane-via-OpenFermion, Block2, ...). The
output is either a :class:`MolecularHamiltonian` of NumPy tensors in the MO
basis, or a FCIDUMP file.

What OpenQP exposes today
-------------------------
The Python data container (``mol.data[...]``) provides everything needed for
the *one-electron* part of the Hamiltonian and all metadata:

* ``OQP::Hcore`` -- core Hamiltonian in the AO basis (packed triangular)
* ``OQP::SM``    -- overlap (packed triangular)
* ``OQP::VEC_MO_A`` / ``OQP::VEC_MO_B`` -- MO coefficients
* ``enuc`` -- nuclear repulsion energy
* ``nelec_A`` / ``nelec_B`` -- electron counts

What is still missing (the one follow-up hook)
-----------------------------------------------
The two-electron repulsion integrals (ERIs) are built on the fly inside the
Fortran SCF Fock construction and are **not currently exposed to Python**.
:func:`from_openqp` therefore accepts the AO ERIs explicitly (``eri_ao=``) or
via a callable ``eri_provider``. The natural completion of this feature is a
small CFFI getter -- ``oqp.int2e(mol)`` populating an ``OQP::ERI_AO`` tag --
mirroring the existing ``oqp.int1e``; once that lands, :func:`from_openqp`
will produce a full FCIDUMP with no external integral source.
"""

from dataclasses import dataclass, field

import numpy as np

from oqp.quantum.integrals import (
    unpack_triangular,
    ao_to_mo_1body,
    ao_to_mo_2body,
)
from oqp.quantum.fcidump import write_fcidump


@dataclass
class MolecularHamiltonian:
    """Second-quantized electronic Hamiltonian in the MO basis.

    Attributes
    ----------
    one_body : numpy.ndarray, shape (norb, norb)
        ``h_pq`` one-electron integrals.
    two_body : numpy.ndarray or None, shape (norb, norb, norb, norb)
        ``(pq|rs)`` two-electron integrals (chemist notation). ``None`` when
        ERIs were unavailable (one-electron-only export).
    core_energy : float
        Scalar energy (nuclear repulsion + any frozen-core contribution).
    n_electrons : int
        Number of correlated electrons.
    ms2 : int
        ``2 * S_z`` = (n_alpha - n_beta).
    orbsym : list of int
        Per-orbital symmetry labels (defaults to all 1).
    """

    one_body: np.ndarray
    core_energy: float
    n_electrons: int
    two_body: np.ndarray = None
    ms2: int = 0
    orbsym: list = field(default_factory=list)

    @property
    def n_orbitals(self):
        return self.one_body.shape[0]

    def to_fcidump(self, filename, tol=1.0e-12):
        """Write this Hamiltonian to a FCIDUMP file.

        Requires two-electron integrals to be present.
        """
        if self.two_body is None:
            raise ValueError(
                "Cannot write a FCIDUMP without two-electron integrals. "
                "Pass eri_ao=/eri_provider= to from_openqp (see module docs "
                "for the planned oqp.int2e hook).")
        orbsym = self.orbsym or [1] * self.n_orbitals
        write_fcidump(
            filename, self.one_body, self.two_body, self.core_energy,
            self.n_electrons, ms2=self.ms2, orbsym=orbsym, tol=tol)
        return filename


def _get_nbf(mol):
    basis = mol.data.get_basis()
    if not basis:
        raise RuntimeError("No basis available; run/apply a basis first.")
    return int(basis["nbf"])


def from_openqp(mol, eri_ao=None, eri_provider=None, mo_coeff=None,
                spin="alpha"):
    """Construct a :class:`MolecularHamiltonian` from an OpenQP ``Molecule``.

    Parameters
    ----------
    mol : oqp.molecule.molecule.Molecule
        A molecule whose SCF has completed (MO coefficients populated).
    eri_ao : array_like, optional
        Two-electron AO integrals ``(mu nu|la si)`` in chemist notation,
        shape ``(nao, nao, nao, nao)``. If given, the two-body MO tensor is
        built and a full FCIDUMP can be written.
    eri_provider : callable, optional
        ``eri_provider(mol) -> eri_ao``; used when ``eri_ao`` is not supplied.
        Lets callers plug in OpenQP's native ERIs (once exposed) or an
        external engine without changing this code.
    mo_coeff : array_like, optional
        Override the MO coefficient matrix (AO rows, MO columns). Defaults to
        the converged restricted/alpha MOs (``OQP::VEC_MO_A``).
    spin : {"alpha", "beta"}
        Which set of converged MOs to use when ``mo_coeff`` is not given.

    Returns
    -------
    MolecularHamiltonian
    """
    nbf = _get_nbf(mol)

    # --- one-electron (always available) ---------------------------------
    hcore_ao = unpack_triangular(mol.data["OQP::Hcore"], nbf)
    if mo_coeff is None:
        tag = "OQP::VEC_MO_A" if spin == "alpha" else "OQP::VEC_MO_B"
        mo_coeff = mol._mo_coefficients(tag, nbf)
    mo_coeff = np.asarray(mo_coeff, dtype=float)

    h1_mo = ao_to_mo_1body(hcore_ao, mo_coeff)

    # --- metadata ---------------------------------------------------------
    na = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
    nb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])
    n_electrons = na + nb
    ms2 = na - nb

    try:
        ecore = float(mol.data["enuc"])
    except Exception:
        ecore = float(getattr(mol.data._data.mol_energy, "enuc", 0.0))

    # --- two-electron (optional until oqp.int2e exists) -------------------
    h2_mo = None
    if eri_ao is None and eri_provider is not None:
        eri_ao = eri_provider(mol)
    if eri_ao is not None:
        h2_mo = ao_to_mo_2body(eri_ao, mo_coeff)

    return MolecularHamiltonian(
        one_body=h1_mo,
        two_body=h2_mo,
        core_energy=ecore,
        n_electrons=n_electrons,
        ms2=ms2,
        orbsym=[1] * h1_mo.shape[0],
    )
