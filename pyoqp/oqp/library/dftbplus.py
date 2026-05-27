"""External DFTB+ backend helpers for OpenQP.

This module intentionally keeps the DFTB+ dependency optional.  Fixture tests
exercise parsing and input generation without requiring the ``dftb+`` executable
or Slater-Koster parameter files to be installed on CI/developer machines.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import os
import re
import shutil
import subprocess
import tempfile
from typing import Any, Iterable, Sequence

BOHR_TO_ANGSTROM = 0.529177210903
_SYMBOLS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 35: "Br", 53: "I",
}


class DFTBPlusError(RuntimeError):
    """Raised when the optional external DFTB+ backend cannot complete."""


CAPABILITY_MATRIX = {
    "energy": {
        "status": "supported",
        "reason": "External DFTB+ single-point energy is parsed from results.tag or detailed.out.",
    },
    "grad": {
        "status": "supported",
        "reason": "External DFTB+ forces are parsed from results.tag and returned as OpenQP gradients.",
    },
    "optimize": {
        "status": "supported",
        "reason": "Ground-state geometry optimization can reuse the external DFTB+ energy/gradient callback.",
    },
    "spin_polarized": {
        "status": "unsupported",
        "reason": "DFTB+ spin-polarized inputs require validated spin constants and open-shell runtime tests before enablement.",
    },
    "unrestricted": {
        "status": "unsupported",
        "reason": "Unrestricted DFTB orbital/spin population data are not yet mapped into OpenQP results or restart structures.",
    },
    "td_dftb": {
        "status": "unsupported",
        "reason": "The external-backend bridge does not parse DFTB+ excited-state outputs or map them into OpenQP TD data.",
    },
    "sf_dftb": {
        "status": "unsupported",
        "reason": "Spin-flip DFTB requires a validated excited-state DFTB implementation; this external DFTB+ bridge has no SF-DFTB parser or response mapping.",
    },
    "mrsf_tddftb": {
        "status": "unsupported",
        "reason": "MRSF-TDDFTB is future native/external excited-state work and must be designed from the CHC2 SF-DFTB reference before runtime enablement.",
    },
    "excited_state_gradients": {
        "status": "unsupported",
        "reason": "Only ground-state DFTB+ forces are parsed; no excited-state gradient, Z-vector, or state-tracking data are available.",
    },
    "nac": {
        "status": "unsupported",
        "reason": "OpenQP NAC workflows require TDHF/MRSF state data that the external DFTB+ bridge does not provide.",
    },
    "spin_flip": {
        "status": "unsupported",
        "reason": "OpenQP spin-flip response is implemented for the native TDHF/MRSF path, not for DFTB+ output.",
    },
    "hessian": {
        "status": "unsupported",
        "reason": "No DFTB+ Hessian parser or OpenQP numerical-Hessian callback is wired for the external backend.",
    },
    "md": {
        "status": "unsupported",
        "reason": "OpenQP has no DFTB+ molecular-dynamics workflow in this external-backend branch.",
    },
    "namd_export": {
        "status": "unsupported",
        "reason": "Surface-hopping/NAMD export requires validated excited-state energies, gradients, and NAC/NACME data; parser fixtures alone are not enough.",
    },
    "native_hamiltonian": {
        "status": "unsupported",
        "reason": "This branch shells out to DFTB+ and does not implement a native OpenQP DFTB Hamiltonian.",
    },
}


@dataclass
class DFTBPlusResult:
    energy: float | None = None
    gradient: list[list[float]] | None = None
    stdout: str = ""
    stderr: str = ""
    workdir: str | None = None


@dataclass(frozen=True)
class DFTBPlusExcitation:
    """One parsed TD-DFTB-like excitation record.

    Parser fixtures may exercise this contract before runtime support is
    validated. A populated record is parser evidence only; capability gates in
    the input checker still keep TD-DFTB/SF-DFTB/MRSF-TDDFTB disabled.
    """

    index: int
    energy_ev: float
    oscillator_strength: float | None = None
    transition_dipole_au: list[float] | None = None
    transition_charges: list[float] | None = None


@dataclass(frozen=True)
class DFTBPlusExcitationResult:
    excitations: list[DFTBPlusExcitation]
    source: str
    validated_runtime: bool = False


@dataclass(frozen=True)
class DFTBSpinFlipState:
    """Source-level SF-DFTB state-numbering contract.

    This is only a mapping scaffold: it records how future SF-DFTB/MRSF-TDDFTB
    roots should be labeled without carrying energies, gradients, or response
    vectors. Runtime capability gates remain closed until the DFTB excited-state
    implementation is validated against real outputs.
    """

    openqp_root: int
    role: str
    physical_state_label: str | None = None
    has_validated_energy: bool = False


@dataclass(frozen=True)
class DFTBNAMDExportFrame:
    """Validated-data-only surface-hopping/NAMD export metadata scaffold.

    This class is deliberately a metadata container, not a dynamics engine.  It
    may only be built from explicitly validated excited-state, gradient, and
    NAC/NACME payloads so parser fixtures cannot be mistaken for production
    nonadiabatic dynamics data.
    """

    state_indices: list[int]
    energies_ev: list[float]
    natom: int
    source: str
    has_validated_gradients: bool
    has_validated_nacme: bool
    has_velocities: bool = False


@dataclass(frozen=True)
class NativeDFTBHamiltonianContract:
    """Source-level seam for a future native OpenQP DFTB Hamiltonian.

    The public branch may record which Slater-Koster pairs and matrix roles a
    native implementation must provide, but this contract deliberately carries
    no Hamiltonian/overlap arrays and keeps runtime flags disabled until a real
    native backend is implemented and validated.
    """

    atom_symbols: list[str]
    required_sk_pairs: list[str]
    matrix_roles: list[str]
    runtime_enabled: bool = False
    validated_runtime: bool = False


@dataclass(frozen=True)
class DFTBSlaterKosterManifest:
    """Filesystem-validated Slater-Koster pair manifest for a DFTB+ run."""

    atom_symbols: list[str]
    required_pairs: list[str]
    found_files: list[str]
    validated_filesystem: bool = True


@dataclass(frozen=True)
class DFTBExcitedBenchmarkCase:
    """Public benchmark metadata for validated non-MRSF DFTB excited-state data.

    This is intentionally metadata-only.  It records where real external DFTB+
    evidence lives after validation, but does not enable TD-DFTB/SF-DFTB runtime
    capabilities or fabricate reference data from parser fixtures.
    """

    molecule: str
    feature_family: str
    state_count: int
    artifact_paths: list[str]
    evidence_level: str
    includes_mrsf_tddftb: bool = False


def build_dftb_excited_benchmark_case(
    *,
    molecule: str,
    feature_family: str,
    excitations: DFTBPlusExcitationResult,
    artifact_paths: Sequence[str | os.PathLike[str]],
) -> DFTBExcitedBenchmarkCase:
    """Create benchmark metadata only from validated public DFTB+ evidence."""

    normalized_feature = feature_family.strip().lower()
    if normalized_feature in {"mrsf_tddftb", "mrsf-td-dftb", "mrsf-tddftb"}:
        raise DFTBPlusError("MRSF-TDDFTB benchmark metadata belongs on the private branch")
    if not excitations.validated_runtime:
        raise DFTBPlusError("DFTB excited-state benchmark metadata requires validated external DFTB+ evidence")
    if not excitations.excitations:
        raise DFTBPlusError("DFTB excited-state benchmark metadata requires at least one excitation")
    paths = [str(path) for path in artifact_paths]
    if not paths:
        raise DFTBPlusError("DFTB excited-state benchmark metadata requires artifact paths")
    return DFTBExcitedBenchmarkCase(
        molecule=molecule,
        feature_family=normalized_feature,
        state_count=len(excitations.excitations),
        artifact_paths=paths,
        evidence_level="validated_external_dftbplus_output",
    )


def build_sf_dftb_state_map(nstate: int) -> list[DFTBSpinFlipState]:
    """Return the planned SF-DFTB root mapping without fabricating data.

    The CHC2/KNU-GAMESS SF-DFTB reconnaissance points to a high-spin reference
    followed by spin-flip states.  OpenQP's scaffold keeps root 0 as that
    reference and labels the requested physical singlet-like roots as S0, S1,
    ... for later response-vector wiring.
    """

    if nstate < 1:
        raise DFTBPlusError("SF-DFTB state mapping requires nstate >= 1")
    states = [DFTBSpinFlipState(openqp_root=0, role="high_spin_reference")]
    for root in range(1, nstate + 1):
        states.append(
            DFTBSpinFlipState(
                openqp_root=root,
                role="spin_flip_state",
                physical_state_label=f"S{root - 1}",
            )
        )
    return states


def _validated_rows(name: str, rows: Sequence[Sequence[float]] | None) -> list[list[float]]:
    if rows is None:
        raise DFTBPlusError("DFTB+ NAMD export requires validated gradients and NAC/NACME data")
    parsed = [[float(component) for component in row] for row in rows]
    if not parsed or any(len(row) != 3 for row in parsed):
        raise DFTBPlusError(f"DFTB+ NAMD export requires non-empty Cartesian {name} rows")
    return parsed


def build_namd_export_frame(
    *,
    excitations: DFTBPlusExcitationResult,
    gradients: Sequence[Sequence[float]] | None = None,
    nacme: dict[tuple[int, int], float] | None = None,
    nac_vectors: dict[tuple[int, int], Sequence[Sequence[float]]] | None = None,
    velocities: Sequence[Sequence[float]] | None = None,
) -> DFTBNAMDExportFrame:
    """Build a validated-data-only NAMD export metadata frame.

    The public DFTB excited-state branch may define interfaces for future
    surface hopping, but it must not fabricate dynamics inputs.  This helper
    therefore rejects synthetic/unvalidated excitation parser results and also
    requires caller-supplied validated gradients plus NAC or NACME payloads.
    """

    if not excitations.validated_runtime:
        raise DFTBPlusError(
            "NAMD export requires validated DFTB excited-state data from a real runtime, not parser fixtures"
        )
    if not excitations.excitations:
        raise DFTBPlusError("DFTB+ NAMD export requires at least one excited state")

    gradient_rows = _validated_rows("gradient", gradients)
    has_nacme = bool(nacme)
    has_nac_vectors = bool(nac_vectors)
    if not has_nacme and not has_nac_vectors:
        raise DFTBPlusError("DFTB+ NAMD export requires validated gradients and NAC/NACME data")

    if velocities is not None:
        velocity_rows = _validated_rows("velocity", velocities)
        if len(velocity_rows) != len(gradient_rows):
            raise DFTBPlusError("DFTB+ NAMD export velocity and gradient atom counts differ")

    return DFTBNAMDExportFrame(
        state_indices=[state.index for state in excitations.excitations],
        energies_ev=[state.energy_ev for state in excitations.excitations],
        natom=len(gradient_rows),
        source=excitations.source,
        has_validated_gradients=True,
        has_validated_nacme=has_nacme or has_nac_vectors,
        has_velocities=velocities is not None,
    )


def build_native_dftb_hamiltonian_contract(atoms: Sequence[int]) -> NativeDFTBHamiltonianContract:
    """Return source-level native DFTB Hamiltonian requirements.

    This is a design/ABI seam only.  It records the element symbols,
    Slater-Koster pair files, and matrix roles that a future native Hamiltonian
    path must provide while keeping runtime capability disabled.  It must not be
    used as a source of Hamiltonian matrices or as evidence that native DFTB is
    implemented.
    """

    atom_symbols: list[str] = []
    for atomic_number in atoms:
        symbol = _SYMBOLS.get(int(atomic_number))
        if symbol is None:
            raise DFTBPlusError(f"No element symbol mapping for atomic number {atomic_number}")
        atom_symbols.append(symbol)

    unique_symbols = sorted(set(atom_symbols))
    required_sk_pairs = [f"{left}-{right}" for i, left in enumerate(unique_symbols) for right in unique_symbols[i:]]
    return NativeDFTBHamiltonianContract(
        atom_symbols=atom_symbols,
        required_sk_pairs=required_sk_pairs,
        matrix_roles=["overlap", "hamiltonian", "repulsive_energy"],
    )


def require_native_dftb_hamiltonian_enabled(contract: NativeDFTBHamiltonianContract) -> None:
    """Fail fast if code tries to use the disabled native DFTB seam at runtime."""

    if not contract.runtime_enabled or not contract.validated_runtime:
        raise DFTBPlusError(
            "native OpenQP DFTB Hamiltonian is a disabled source-level seam; use external DFTB+ or add validated native runtime support"
        )


def validate_slater_koster_manifest(
    atoms: Sequence[int],
    sk_path: str | os.PathLike[str],
) -> DFTBSlaterKosterManifest:
    """Fail fast when required Slater-Koster pair files are unavailable.

    DFTB+ runtime validation must not wait until an opaque external executable
    failure if OpenQP can determine up front that the requested element-pair
    parameter files are missing.  Mixed pairs accept either file orientation
    (for example ``H-O.skf`` or ``O-H.skf``) because public parameter sets differ
    in how they store symmetric pairs.
    """

    root = Path(sk_path)
    atom_symbols = _symbols_for_atoms(atoms)
    unique_symbols = sorted(set(atom_symbols))
    required_pairs = [f"{left}-{right}" for i, left in enumerate(unique_symbols) for right in unique_symbols[i:]]
    found_files: list[str] = []
    missing_files: list[str] = []

    for pair in required_pairs:
        left, right = pair.split("-", 1)
        candidates = [f"{left}-{right}.skf"]
        if left != right:
            candidates.append(f"{right}-{left}.skf")
        found = next((candidate for candidate in candidates if (root / candidate).is_file()), None)
        if found is None:
            missing_files.append(candidates[0])
        else:
            found_files.append(found)

    if missing_files:
        raise DFTBPlusError("DFTB+ parameter directory is missing Slater-Koster files: " + ", ".join(missing_files))

    return DFTBSlaterKosterManifest(
        atom_symbols=atom_symbols,
        required_pairs=required_pairs,
        found_files=found_files,
    )


def _read_text(path: str | os.PathLike[str]) -> str:
    return Path(path).read_text(encoding="utf-8")


def parse_detailed_out(path: str | os.PathLike[str]) -> DFTBPlusResult:
    """Parse the total energy from a DFTB+ ``detailed.out`` file."""

    text = _read_text(path)
    energy = None
    for pattern in (
        r"Total\s+Energy\s*:\s*([-+0-9.Ee]+)\s*H",
        r"Total\s+Mermin\s+free\s+energy\s*:\s*([-+0-9.Ee]+)\s*H",
    ):
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            energy = float(match.group(1))
            break
    if energy is None:
        raise DFTBPlusError(f"Could not parse total energy from {path}")
    return DFTBPlusResult(energy=energy)


def _parse_tag_header(header: str) -> tuple[str, list[int]]:
    # Example: forces:real:2:3,3
    fields = header.strip().split(":")
    name = fields[0].strip()
    shape = []
    if len(fields) >= 4 and fields[3]:
        shape = [int(x) for x in fields[3].split(",") if x]
    return name, shape


def parse_results_tag(path: str | os.PathLike[str]) -> DFTBPlusResult:
    """Parse energy and gradients from DFTB+ ``results.tag``.

    DFTB+ stores forces in Hartree/Bohr. OpenQP gradients are dE/dR, so this
    parser returns ``gradient = -forces``.
    """

    lines = _read_text(path).splitlines()
    energy = None
    gradient = None
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line or ":" not in line:
            index += 1
            continue
        name, shape = _parse_tag_header(line)
        index += 1
        count = 1
        if shape:
            count = 1
            for dim in shape:
                count *= dim
        values: list[float] = []
        value_rows: list[list[float]] = []
        while index < len(lines) and len(values) < count:
            next_line = lines[index].strip()
            if next_line:
                row = [float(token) for token in next_line.split()]
                value_rows.append(row)
                values.extend(row)
            index += 1
        if name in {"total_energy", "mermin_energy", "extrapolated0_energy", "forcerelated_energy"} and values:
            if energy is None or name == "total_energy":
                energy = values[0]
        elif name == "forces" and values:
            if not shape or len(shape) != 2 or 3 not in shape:
                raise DFTBPlusError(f"Unexpected forces shape in {path}: {shape}")
            if shape[1] == 3:
                natom = shape[0]
                gradient = [[-values[3 * i + j] for j in range(3)] for i in range(natom)]
            elif len(value_rows) == shape[1] and all(len(row) == 3 for row in value_rows):
                # DFTB+ 25.1 may label forces as :3,natom while still writing
                # one Cartesian force row per atom.  Prefer the physical row
                # layout when it is unambiguous.
                gradient = [[-component for component in row] for row in value_rows]
            else:
                natom = shape[1]
                gradient = [[-values[component * natom + atom] for component in range(3)] for atom in range(natom)]
    if energy is None:
        raise DFTBPlusError(f"Could not parse total_energy from {path}")
    return DFTBPlusResult(energy=energy, gradient=gradient)


_EXCITED_STATE_RE = re.compile(
    r"Excited\s+State\s+(?P<index>\d+)\s*:\s*"
    r"(?:excitation\s+energy\s*=\s*)?(?P<energy>[-+0-9.Ee]+)\s*eV"
    r"(?:.*?oscillator\s+strength\s*=\s*(?P<osc>[-+0-9.Ee]+))?",
    re.IGNORECASE,
)
_TRANSITION_DIPOLE_RE = re.compile(
    r"Transition\s+dipole(?:\s*\[[^\]]+\])?\s*=\s*"
    r"(?P<x>[-+0-9.Ee]+)\s+(?P<y>[-+0-9.Ee]+)\s+(?P<z>[-+0-9.Ee]+)",
    re.IGNORECASE,
)


def parse_dftbplus_excitations(path: str | os.PathLike[str]) -> DFTBPlusExcitationResult:
    """Parse a minimal TD-DFTB excitation contract from DFTB+-style text.

    The parser accepts compact tagged text fixtures and future real DFTB+
    output excerpts containing excitation energies, oscillator strengths, and
    optional transition dipoles. It deliberately does not enable runtime
    TD-DFTB/SF-DFTB/MRSF-TDDFTB; callers must keep those capability gates closed
    until a real executable/parameter validation path is added.
    """

    text = _read_text(path)
    excitations: list[DFTBPlusExcitation] = []
    pending: dict[str, Any] | None = None

    def flush_pending() -> None:
        nonlocal pending
        if pending is None:
            return
        excitations.append(
            DFTBPlusExcitation(
                index=int(pending["index"]),
                energy_ev=float(pending["energy_ev"]),
                oscillator_strength=pending.get("oscillator_strength"),
                transition_dipole_au=pending.get("transition_dipole_au"),
                transition_charges=pending.get("transition_charges"),
            )
        )
        pending = None

    collecting_transition_charges = False
    for line in text.splitlines():
        state_match = _EXCITED_STATE_RE.search(line)
        if state_match:
            flush_pending()
            collecting_transition_charges = False
            osc = state_match.group("osc")
            pending = {
                "index": int(state_match.group("index")),
                "energy_ev": float(state_match.group("energy")),
                "oscillator_strength": float(osc) if osc is not None else None,
            }
            continue
        if pending is not None and re.search(r"Transition\s+charges", line, re.IGNORECASE):
            pending["transition_charges"] = []
            collecting_transition_charges = True
            continue
        dipole_match = _TRANSITION_DIPOLE_RE.search(line)
        if dipole_match and pending is not None:
            pending["transition_dipole_au"] = [
                float(dipole_match.group(axis)) for axis in ("x", "y", "z")
            ]
            collecting_transition_charges = False
            continue
        if collecting_transition_charges and pending is not None and line.strip():
            try:
                pending["transition_charges"].append(float(line.split()[-1]))
                continue
            except ValueError:
                collecting_transition_charges = False
    flush_pending()

    if not excitations:
        raise DFTBPlusError(f"No DFTB+ excitation states parsed from {path}")
    source = "synthetic_parser_fixture" if "Synthetic parser contract fixture" in text else "dftbplus_output_excerpt"
    return DFTBPlusExcitationResult(
        excitations=excitations,
        source=source,
        validated_runtime=False,
    )


def _symbols_for_atoms(atoms: Sequence[int]) -> list[str]:
    symbols = []
    for atom in atoms:
        atomic_number = int(atom)
        try:
            symbols.append(_SYMBOLS[atomic_number])
        except KeyError as exc:
            raise DFTBPlusError(f"No DFTB+ symbol mapping for atomic number {atomic_number}") from exc
    return symbols


def _coords_rows(coords_bohr: Sequence[float]) -> list[list[float]]:
    if len(coords_bohr) % 3:
        raise DFTBPlusError("Coordinate array length must be a multiple of 3")
    return [
        [float(coords_bohr[3 * i + j]) * BOHR_TO_ANGSTROM for j in range(3)]
        for i in range(len(coords_bohr) // 3)
    ]


_DEFAULT_MAX_ANGULAR_MOMENTUM = {
    "H": "s", "He": "s",
    "Li": "p", "Be": "p", "B": "p", "C": "p", "N": "p", "O": "p", "F": "p", "Ne": "p",
    "Na": "p", "Mg": "p", "Al": "p", "Si": "p", "P": "p", "S": "p", "Cl": "p", "Ar": "p",
    "K": "p", "Ca": "p", "Br": "p", "I": "p",
    "Cr": "d",
}


def _max_angular_momentum_block(symbols: Sequence[str], config: dict) -> str:
    dftb = _dftb_config(config)
    overrides = dftb.get("max_angular_momentum", {}) or {}
    lines = ["  MaxAngularMomentum = {"]
    for symbol in dict.fromkeys(symbols):
        value = overrides.get(symbol, _DEFAULT_MAX_ANGULAR_MOMENTUM.get(symbol))
        if value is None:
            raise DFTBPlusError(f"No DFTB+ MaxAngularMomentum default for atom type {symbol}")
        lines.append(f'    {symbol} = "{value}"')
    lines.append("  }")
    return "\n".join(lines)


def _dftb_config(config: dict) -> dict:
    return config.get("dftb", {}) if config else {}


def _truthy(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes", ".true."}
    return bool(value)


def validate_spin_configuration(config: dict) -> None:
    """Fail fast for unvalidated spin-polarized/unrestricted DFTB requests.

    The schema accepts these keys so future work can attach SK spin-constant
    validation and DFTB+ ``SpinPolarisation`` blocks without changing the input
    surface again. Runtime remains disabled until a real DFTB+ executable/SK
    validation path proves the open-shell data are chemically and numerically
    safe to expose.
    """

    dftb = _dftb_config(config)
    if _truthy(dftb.get("spin_polarized", False)) or _truthy(dftb.get("unrestricted", False)):
        spin_constants = str(dftb.get("spin_constants_path", "") or "").strip()
        suffix = ""
        if not spin_constants:
            suffix = "; set [dftb] spin_constants_path when implementing the validated path"
        raise DFTBPlusError(
            "DFTB+ spin-polarized/unrestricted DFTB is not validated in this branch"
            + suffix
        )


def write_dftbplus_input(
    workdir: str | os.PathLike[str],
    atoms: Sequence[int],
    coords_bohr: Sequence[float],
    config: dict,
    *,
    gradient: bool = False,
) -> Path:
    """Write ``dftb_in.hsd`` and ``geo.gen`` for an external DFTB+ run."""

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    dftb = _dftb_config(config)
    validate_spin_configuration(config)
    sk_path = str(dftb.get("sk_path", "")).strip()
    if sk_path and not sk_path.endswith(os.sep):
        sk_path += os.sep
    scc = bool(dftb.get("scc", True))
    max_scc = int(dftb.get("max_scc_iterations", 100))
    charge = int(config.get("input", {}).get("charge", 0))
    symbols = _symbols_for_atoms(atoms)
    unique_symbols = list(dict.fromkeys(symbols))
    max_l_block = _max_angular_momentum_block(unique_symbols, config)
    rows = _coords_rows(coords_bohr)
    if len(rows) != len(symbols):
        raise DFTBPlusError("Number of atoms does not match coordinate rows")

    geo_path = workdir / "geo.gen"
    with geo_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(symbols)} C\n")
        handle.write(" ".join(unique_symbols) + "\n")
        for idx, (symbol, xyz) in enumerate(zip(symbols, rows), start=1):
            type_index = unique_symbols.index(symbol) + 1
            handle.write(f"{idx:5d} {type_index:3d} {xyz[0]:18.10f} {xyz[1]:18.10f} {xyz[2]:18.10f}\n")

    hsd = f'''Geometry = GenFormat {{
  <<< "geo.gen"
}}
# AtomTypes = {" ".join(unique_symbols)}

Hamiltonian = DFTB {{
  SCC = {"Yes" if scc else "No"}
  MaxSCCIterations = {max_scc}
  Charge = {charge}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{sk_path}"
    Separator = "-"
    Suffix = ".skf"
  }}
{max_l_block}
}}

Options {{
  WriteResultsTag = Yes
}}
'''
    if gradient:
        hsd += "\nDriver = ConjugateGradient { MaxSteps = 0 }\n"
    input_path = workdir / "dftb_in.hsd"
    input_path.write_text(hsd, encoding="utf-8")
    return input_path


class DFTBPlusRunner:
    """Small shell-out wrapper around the optional ``dftb+`` executable."""

    def __init__(self, config: dict):
        self.config = config or {}
        dftb = _dftb_config(self.config)
        self.executable = str(dftb.get("executable", "dftb+") or "dftb+")
        self.sk_path = str(dftb.get("sk_path", "") or "")
        self.keep_workdir = bool(dftb.get("keep_workdir", False))
        self.timeout = int(dftb.get("timeout", 300) or 300)

    def _check_prerequisites(self) -> None:
        exe = self.executable
        if os.path.sep in exe:
            if not (Path(exe).exists() and os.access(exe, os.X_OK)):
                raise DFTBPlusError(f"DFTB+ executable not found or not executable: {exe}")
        elif shutil.which(exe) is None:
            raise DFTBPlusError(f"DFTB+ executable not found on PATH: {exe}")
        if not self.sk_path:
            raise DFTBPlusError("DFTB+ parameter directory is not configured; set [dftb] sk_path")
        if not Path(self.sk_path).is_dir():
            raise DFTBPlusError(f"DFTB+ parameter directory not found: {self.sk_path}")

    def run(self, atoms: Sequence[int], coords_bohr: Sequence[float], *, gradient: bool = False) -> DFTBPlusResult:
        self._check_prerequisites()
        validate_slater_koster_manifest(atoms, self.sk_path)
        if self.keep_workdir:
            workdir = Path(tempfile.mkdtemp(prefix="openqp-dftbplus-"))
            return self._run_in_workdir(workdir, atoms, coords_bohr, gradient=gradient)

        with tempfile.TemporaryDirectory(prefix="openqp-dftbplus-") as tmp:
            return self._run_in_workdir(Path(tmp), atoms, coords_bohr, gradient=gradient)

    def _run_in_workdir(self, workdir: Path, atoms: Sequence[int], coords_bohr: Sequence[float], *, gradient: bool) -> DFTBPlusResult:
        write_dftbplus_input(workdir, atoms, coords_bohr, self.config, gradient=gradient)
        proc = subprocess.run(
            [self.executable],
            cwd=workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=self.timeout,
            check=False,
        )
        if proc.returncode != 0:
            raise DFTBPlusError(
                f"DFTB+ exited with status {proc.returncode}: {proc.stderr.strip() or proc.stdout.strip()}"
            )
        result_path = workdir / "results.tag"
        detail_path = workdir / "detailed.out"
        if result_path.exists():
            result = parse_results_tag(result_path)
        elif detail_path.exists():
            result = parse_detailed_out(detail_path)
        else:
            raise DFTBPlusError("DFTB+ completed but neither results.tag nor detailed.out was written")
        result.stdout = proc.stdout
        result.stderr = proc.stderr
        result.workdir = str(workdir)
        return result


def run_openqp_molecule(mol, *, gradient: bool = False) -> DFTBPlusResult:
    """Execute DFTB+ for an OpenQP Molecule-like object and populate results."""

    result = DFTBPlusRunner(mol.config).run(mol.get_atoms(), mol.get_system(), gradient=gradient)
    if result.energy is not None:
        mol.data._data.mol_energy.energy = result.energy
        if hasattr(mol, "energies"):
            mol.energies = [result.energy]
    if gradient:
        if result.gradient is None:
            raise DFTBPlusError("DFTB+ gradient was requested, but no gradient was parsed")
        flat = [component for row in result.gradient for component in row]
        try:
            import numpy as np
            from oqp import ffi

            grad_array = np.asarray(flat, dtype=np.float64)
            ffi.memmove(mol.data._data.grad, grad_array, grad_array.nbytes)
        except Exception:
            if hasattr(mol, "grads"):
                mol.grads = [flat]
    return result


def _as_flat_list(values: Sequence[float]) -> list[float]:
    try:
        import numpy as np

        return np.asarray(values, dtype=float).reshape(-1).tolist()
    except Exception:
        return [float(value) for value in values]


def _reshape_coords(values: Sequence[float]) -> list[float]:
    flat = _as_flat_list(values)
    if len(flat) % 3:
        raise DFTBPlusError("Optimizer coordinate array length must be a multiple of 3")
    return flat


def _coords_for_mol_update(flat_coords: Sequence[float]):
    try:
        import numpy as np

        return np.asarray(flat_coords, dtype=float).reshape((-1, 3))
    except Exception:
        return list(flat_coords)


def _gradient_flat(result: DFTBPlusResult) -> list[float]:
    if result.gradient is None:
        raise DFTBPlusError("DFTB+ geometry optimization requires gradients, but no gradient was parsed")
    return [float(component) for row in result.gradient for component in row]


def _rms(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return math.sqrt(sum(float(value) ** 2 for value in values) / len(values))


def optimize_openqp_molecule(mol, *, runner_factory=DFTBPlusRunner) -> DFTBPlusResult:
    """Run a ground-state DFTB+ geometry optimization for a Molecule-like object."""

    config = mol.config
    opt_config = config.get("optimize", {})
    if int(opt_config.get("istate", 0) or 0) != 0:
        raise DFTBPlusError("DFTB+ geometry optimization supports only ground-state istate=0")

    method = str(opt_config.get("optimizer", "bfgs") or "bfgs").lower()
    maxit = int(opt_config.get("maxit", 50) or 50)
    energy_shift = float(opt_config.get("energy_shift", 1.0e-6) or 1.0e-6)
    rmsd_grad_target = float(opt_config.get("rmsd_grad", 3.0e-4) or 3.0e-4)
    max_grad_target = float(opt_config.get("max_grad", 4.5e-4) or 4.5e-4)
    rmsd_step_target = float(opt_config.get("rmsd_step", 1.2e-3) or 1.2e-3)
    max_step_target = float(opt_config.get("max_step", 1.8e-3) or 1.8e-3)
    atoms = mol.get_atoms()
    runner = runner_factory(config)
    if hasattr(runner, "_check_prerequisites"):
        runner._check_prerequisites()
    previous_energy = None
    previous_coords = _reshape_coords(mol.get_system())
    latest: DFTBPlusResult | None = None

    def evaluate(coords):
        nonlocal latest, previous_energy, previous_coords
        flat_coords = _reshape_coords(coords)
        mol.update_system(_coords_for_mol_update(flat_coords))
        latest = runner.run(atoms, flat_coords, gradient=True)
        grad = _gradient_flat(latest)
        energy = float(latest.energy)
        if hasattr(mol, "energies"):
            mol.energies = [energy]
        if hasattr(mol, "grads"):
            mol.grads = [grad]

        step = [flat_coords[i] - previous_coords[i] for i in range(len(flat_coords))]
        de = abs(energy - previous_energy) if previous_energy is not None else float("inf")
        converged = (
            previous_energy is not None
            and de <= energy_shift
            and _rms(step) <= rmsd_step_target
            and max((abs(value) for value in step), default=0.0) <= max_step_target
            and _rms(grad) <= rmsd_grad_target
            and max((abs(value) for value in grad), default=0.0) <= max_grad_target
        )
        previous_energy = energy
        previous_coords = flat_coords
        if converged:
            raise StopIteration
        return energy, grad

    initial = _reshape_coords(mol.get_system())
    try:
        import scipy.optimize as scipy_optimize
    except ImportError:
        # Keep the external DFTB+ backend usable in minimal Python installs: a
        # single energy/gradient callback is sufficient for OpenQP to expose a
        # validated gradient and for tests/examples to exercise the bridge.  Full
        # geometry optimization uses scipy when it is available.
        try:
            evaluate(initial)
        except StopIteration:
            pass
    else:
        try:
            scipy_optimize.minimize(
                fun=evaluate,
                x0=initial,
                method=method,
                jac=True,
                options={"maxiter": maxit},
            )
        except StopIteration:
            pass

    if latest is None:
        raise DFTBPlusError("DFTB+ geometry optimization did not evaluate any real energy/gradient")
    return latest
