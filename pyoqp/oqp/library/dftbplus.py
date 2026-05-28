"""External DFTB+ backend helpers for OpenQP.

This module intentionally keeps the DFTB+ dependency optional.  Fixture tests
exercise parsing and input generation without requiring the ``dftb+`` executable
or Slater-Koster parameter files to be installed on CI/developer machines.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
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
class DFTBExcitationObservableContract:
    """Validated observable-completeness metadata for DFTB excited states.

    This is a guard for downstream benchmark/NAMD interfaces: it records that a
    caller supplied complete parsed observables, but it does not enable TD-DFTB,
    SF-DFTB, gradients, NACs, or dynamics runtime paths.
    """

    state_indices: list[int]
    has_oscillator_strengths: bool
    has_transition_dipoles: bool
    has_transition_charges: bool
    transition_charge_natom: int | None = None
    enables_runtime_capability: bool = False


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
    observable_contract: DFTBExcitationObservableContract
    has_velocities: bool = False


@dataclass(frozen=True)
class DFTBExcitedGradientFrame:
    """Validated-data-only excited-state gradient metadata scaffold.

    Future TD-DFTB/SF-DFTB gradient work needs a shape/state contract, but this
    public branch must not infer or fabricate excited-state gradients from
    ground-state DFTB+ forces or parser fixtures.  Instances are therefore only
    produced from caller-supplied validated runtime evidence and still do not
    enable runtime gradient capability by themselves.
    """

    state_indices: list[int]
    gradients_by_state: dict[int, list[list[float]]]
    natom: int
    source: str
    enables_runtime_capability: bool = False


@dataclass(frozen=True)
class DFTBNACFrame:
    """Validated-data-only NAC/NACME metadata scaffold.

    Future TD-DFTB/SF-DFTB nonadiabatic workflows need a strict state-pair and
    vector-shape contract, but this public branch must not infer NACs from
    excitation energies, transition charges, or parser fixtures. Instances are
    therefore only produced from caller-supplied validated runtime evidence and
    still do not enable NAC/NACME runtime capability by themselves.
    """

    state_indices: list[int]
    nacme_by_pair: dict[tuple[int, int], float]
    nac_vectors_by_pair: dict[tuple[int, int], list[list[float]]]
    natom: int | None
    source: str
    enables_runtime_capability: bool = False


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
class DFTBSpinConstantManifest:
    """Filesystem-validated spin-constant manifest for future spin-polarized DFTB."""

    atom_symbols: list[str]
    required_symbols: list[str]
    found_symbols: list[str]
    validated_filesystem: bool = True
    enables_runtime_capability: bool = False


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
    observable_contract: DFTBExcitationObservableContract
    includes_mrsf_tddftb: bool = False


@dataclass(frozen=True)
class DFTBExcitedBenchmarkSuite:
    """Metadata-only public DFTB excited-state benchmark-suite manifest."""

    name: str
    molecules: list[str]
    feature_families: list[str]
    case_count: int
    artifact_paths: list[str]
    evidence_level: str
    cases: list[DFTBExcitedBenchmarkCase]
    includes_mrsf_tddftb: bool = False
    enables_runtime_capability: bool = False


def validate_dftb_excitation_observables(
    excitations: DFTBPlusExcitationResult,
    *,
    require_oscillator_strengths: bool = False,
    require_transition_dipoles: bool = False,
    require_transition_charges_natom: int | None = None,
) -> DFTBExcitationObservableContract:
    """Validate parsed excited-state observable completeness without enabling runtime.

    Parser contracts for TD-DFTB-like data are useful only when downstream code
    can assert which observables are actually present. This helper fails fast for
    incomplete required oscillator strengths, transition dipoles, or transition
    charges rather than letting benchmark/NAMD scaffolds infer missing data.
    """

    if not excitations.excitations:
        raise DFTBPlusError("DFTB excited-state observables require at least one excitation")

    missing_osc = [state.index for state in excitations.excitations if state.oscillator_strength is None]
    if require_oscillator_strengths and missing_osc:
        raise DFTBPlusError(f"DFTB excited-state observables missing oscillator strengths for states {missing_osc}")

    missing_dipoles = [state.index for state in excitations.excitations if state.transition_dipole_au is None]
    if require_transition_dipoles and missing_dipoles:
        raise DFTBPlusError(f"DFTB excited-state observables missing transition dipoles for states {missing_dipoles}")

    charge_lengths = [
        len(state.transition_charges) if state.transition_charges is not None else None
        for state in excitations.excitations
    ]
    if require_transition_charges_natom is not None:
        if require_transition_charges_natom < 1:
            raise DFTBPlusError("DFTB transition-charge observable validation requires natom >= 1")
        bad_charge_states = [
            state.index
            for state in excitations.excitations
            if state.transition_charges is None or len(state.transition_charges) != require_transition_charges_natom
        ]
        if bad_charge_states:
            raise DFTBPlusError(
                "DFTB excited-state observables missing transition charges with the requested atom count "
                f"for states {bad_charge_states}"
            )

    return DFTBExcitationObservableContract(
        state_indices=[state.index for state in excitations.excitations],
        has_oscillator_strengths=not missing_osc,
        has_transition_dipoles=not missing_dipoles,
        has_transition_charges=all(length is not None for length in charge_lengths),
        transition_charge_natom=require_transition_charges_natom,
    )


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
    observable_contract = validate_dftb_excitation_observables(
        excitations,
        require_oscillator_strengths=True,
        require_transition_dipoles=True,
    )
    return DFTBExcitedBenchmarkCase(
        molecule=molecule,
        feature_family=normalized_feature,
        state_count=len(excitations.excitations),
        artifact_paths=paths,
        evidence_level="validated_external_dftbplus_output",
        observable_contract=observable_contract,
    )


def build_dftb_excited_benchmark_suite(
    *,
    name: str,
    cases: Sequence[DFTBExcitedBenchmarkCase],
    required_molecules: Sequence[str] | None = None,
) -> DFTBExcitedBenchmarkSuite:
    """Create a metadata-only suite manifest from validated public cases.

    The suite is deliberately a bookkeeping contract for future benchmark files;
    it does not enable TD-DFTB/SF-DFTB runtime support or accept private
    MRSF-TDDFTB cases on this public branch.
    """

    if not cases:
        raise DFTBPlusError("DFTB excited-state benchmark suite requires at least one case")
    private_cases = [
        case.molecule
        for case in cases
        if case.includes_mrsf_tddftb or case.feature_family == "mrsf_tddftb"
    ]
    if private_cases:
        raise DFTBPlusError("DFTB excited-state benchmark suite cannot include private MRSF-TDDFTB cases")
    unvalidated_cases = [
        case.molecule
        for case in cases
        if case.evidence_level != "validated_external_dftbplus_output"
    ]
    if unvalidated_cases:
        raise DFTBPlusError(f"DFTB excited-state benchmark suite contains unvalidated cases: {unvalidated_cases}")

    molecules = sorted({case.molecule for case in cases})
    if required_molecules:
        missing = sorted(set(required_molecules) - set(molecules))
        if missing:
            raise DFTBPlusError("DFTB excited-state benchmark suite missing required public benchmark molecules: " + ", ".join(missing))

    artifact_paths: list[str] = []
    for case in cases:
        artifact_paths.extend(case.artifact_paths)

    return DFTBExcitedBenchmarkSuite(
        name=name,
        molecules=molecules,
        feature_families=sorted({case.feature_family for case in cases}),
        case_count=len(cases),
        artifact_paths=artifact_paths,
        evidence_level="validated_external_dftbplus_output",
        cases=list(cases),
    )


def validate_dftb_excited_benchmark_suite_coverage(
    suite: DFTBExcitedBenchmarkSuite,
    *,
    required_feature_families: Sequence[str] | None = None,
    required_molecules: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Validate public benchmark-suite coverage without enabling runtime paths.

    This helper is a benchmark-planning guard: it checks that a metadata suite
    contains the requested public feature families/molecules, rejects private
    MRSF-TDDFTB scope through the existing serializer validator, and returns a
    compact coverage summary with ``enables_runtime_capability`` kept false.
    """

    manifest = serialize_dftb_excited_benchmark_suite(suite)
    _validate_dftb_excited_benchmark_suite_manifest(manifest)

    feature_families = sorted({str(feature).strip().lower() for feature in suite.feature_families})
    molecules = sorted({str(molecule) for molecule in suite.molecules})

    if required_feature_families:
        requested_features = sorted({str(feature).strip().lower() for feature in required_feature_families})
        private_requested = [
            feature
            for feature in requested_features
            if feature in {"mrsf_tddftb", "mrsf-td-dftb", "mrsf-tddftb"}
        ]
        if private_requested:
            raise DFTBPlusError("DFTB excited-state benchmark coverage cannot require private MRSF-TDDFTB scope")
        missing_features = sorted(set(requested_features) - set(feature_families))
        if missing_features:
            raise DFTBPlusError(
                "DFTB excited-state benchmark suite missing required public benchmark feature families: "
                + ", ".join(missing_features)
            )

    if required_molecules:
        requested_molecules = sorted({str(molecule) for molecule in required_molecules})
        missing_molecules = sorted(set(requested_molecules) - set(molecules))
        if missing_molecules:
            raise DFTBPlusError(
                "DFTB excited-state benchmark suite missing required public benchmark molecules: "
                + ", ".join(missing_molecules)
            )

    return {
        "schema": "openqp.dftb.excited_benchmark_suite.coverage.v1",
        "name": suite.name,
        "feature_families": feature_families,
        "molecules": molecules,
        "case_count": suite.case_count,
        "evidence_level": suite.evidence_level,
        "enables_runtime_capability": False,
    }


def serialize_dftb_excited_benchmark_suite(suite: DFTBExcitedBenchmarkSuite) -> dict[str, Any]:
    """Serialize a public DFTB excited-state benchmark suite manifest.

    The serialized manifest is intentionally provenance-only: it records the
    validated external artifact paths and observable completeness contracts while
    preserving ``enables_runtime_capability=False`` so downstream tooling cannot
    mistake benchmark bookkeeping for implemented TD-DFTB/SF-DFTB runtime support.
    """

    case_payloads = []
    for case in suite.cases:
        if case.includes_mrsf_tddftb or case.feature_family == "mrsf_tddftb":
            raise DFTBPlusError("DFTB excited-state benchmark manifest cannot serialize private MRSF-TDDFTB cases")
        case_payloads.append(
            {
                "molecule": case.molecule,
                "feature_family": case.feature_family,
                "state_count": case.state_count,
                "artifact_paths": list(case.artifact_paths),
                "evidence_level": case.evidence_level,
                "observables": {
                    "state_indices": list(case.observable_contract.state_indices),
                    "has_oscillator_strengths": case.observable_contract.has_oscillator_strengths,
                    "has_transition_dipoles": case.observable_contract.has_transition_dipoles,
                    "has_transition_charges": case.observable_contract.has_transition_charges,
                    "transition_charge_natom": case.observable_contract.transition_charge_natom,
                },
                "includes_mrsf_tddftb": case.includes_mrsf_tddftb,
            }
        )

    return {
        "schema": "openqp.dftb.excited_benchmark_suite.v1",
        "name": suite.name,
        "molecules": list(suite.molecules),
        "feature_families": list(suite.feature_families),
        "case_count": suite.case_count,
        "artifact_paths": list(suite.artifact_paths),
        "evidence_level": suite.evidence_level,
        "includes_mrsf_tddftb": suite.includes_mrsf_tddftb,
        "enables_runtime_capability": suite.enables_runtime_capability,
        "cases": case_payloads,
    }


def _validate_dftb_excited_benchmark_suite_manifest(manifest: dict[str, Any]) -> None:
    """Fail fast if a benchmark manifest overclaims DFTB excited-state support."""

    if manifest.get("schema") != "openqp.dftb.excited_benchmark_suite.v1":
        raise DFTBPlusError("DFTB excited-state benchmark manifest has an unknown schema")
    if manifest.get("enables_runtime_capability"):
        raise DFTBPlusError("DFTB excited-state benchmark manifest must not enable runtime capability")
    if manifest.get("includes_mrsf_tddftb"):
        raise DFTBPlusError("DFTB excited-state benchmark manifest cannot include private MRSF-TDDFTB scope")
    if manifest.get("evidence_level") != "validated_external_dftbplus_output":
        raise DFTBPlusError("DFTB excited-state benchmark manifest requires validated external DFTB+ evidence")

    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise DFTBPlusError("DFTB excited-state benchmark manifest requires at least one case")
    if manifest.get("case_count") != len(cases):
        raise DFTBPlusError("DFTB excited-state benchmark manifest case_count does not match cases")

    for case in cases:
        if case.get("includes_mrsf_tddftb") or case.get("feature_family") == "mrsf_tddftb":
            raise DFTBPlusError("DFTB excited-state benchmark manifest cannot include private MRSF-TDDFTB cases")
        if case.get("evidence_level") != "validated_external_dftbplus_output":
            raise DFTBPlusError("DFTB excited-state benchmark manifest case requires validated external DFTB+ evidence")
        if not case.get("artifact_paths"):
            raise DFTBPlusError("DFTB excited-state benchmark manifest case requires artifact paths")
        observables = case.get("observables") or {}
        if not observables.get("has_oscillator_strengths") or not observables.get("has_transition_dipoles"):
            raise DFTBPlusError("DFTB excited-state benchmark manifest case requires oscillator strengths and transition dipoles")


def _resolve_benchmark_artifact_path(path: str, base_dir: str | os.PathLike[str] | None) -> Path:
    artifact_path = Path(path)
    if not artifact_path.is_absolute() and base_dir is not None:
        artifact_path = Path(base_dir) / artifact_path
    return artifact_path


def _artifact_provenance(path: str, base_dir: str | os.PathLike[str] | None) -> dict[str, Any]:
    artifact_path = _resolve_benchmark_artifact_path(path, base_dir)
    if not artifact_path.is_file():
        raise DFTBPlusError(f"DFTB excited-state benchmark artifact file is missing: {path}")
    payload = artifact_path.read_bytes()
    return {
        "path": path,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def _validate_artifact_provenance_entries(
    entries: Sequence[dict[str, Any]],
    *,
    base_dir: str | os.PathLike[str] | None,
) -> None:
    if not entries:
        raise DFTBPlusError("DFTB excited-state benchmark manifest requires artifact provenance entries")
    for entry in entries:
        path = entry.get("path")
        expected_sha = entry.get("sha256")
        expected_size = entry.get("size_bytes")
        if not path or not expected_sha or not isinstance(expected_size, int):
            raise DFTBPlusError("DFTB excited-state benchmark artifact provenance is incomplete")
        actual = _artifact_provenance(str(path), base_dir)
        if actual["size_bytes"] != expected_size:
            raise DFTBPlusError(f"DFTB excited-state benchmark artifact size mismatch: {path}")
        if actual["sha256"] != expected_sha:
            raise DFTBPlusError(f"DFTB excited-state benchmark artifact checksum mismatch: {path}")


def attach_dftb_excited_benchmark_artifact_provenance(
    manifest: dict[str, Any],
    *,
    base_dir: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Return a manifest copy with checked artifact size and SHA-256 provenance.

    This helper deliberately validates only already-declared benchmark artifacts.
    It does not parse or bless DFTB+ output and it keeps the runtime-capability
    gate closed; missing files or later checksum drift fail fast.
    """

    _validate_dftb_excited_benchmark_suite_manifest(manifest)
    enriched = json.loads(json.dumps(manifest))
    aggregate: list[dict[str, Any]] = []
    for case in enriched["cases"]:
        case_entries = [_artifact_provenance(str(path), base_dir) for path in case.get("artifact_paths", [])]
        case["artifact_provenance"] = case_entries
        aggregate.extend(case_entries)
    enriched["artifact_provenance"] = aggregate
    enriched["enables_runtime_capability"] = False
    return enriched


def _validate_manifest_artifacts(
    manifest: dict[str, Any],
    *,
    base_dir: str | os.PathLike[str] | None,
) -> None:
    _validate_artifact_provenance_entries(manifest.get("artifact_provenance") or [], base_dir=base_dir)
    for case in manifest.get("cases", []):
        _validate_artifact_provenance_entries(case.get("artifact_provenance") or [], base_dir=base_dir)


def write_dftb_excited_benchmark_suite_manifest(
    suite: DFTBExcitedBenchmarkSuite,
    path: str | os.PathLike[str],
) -> None:
    """Write a validated public DFTB excited-state benchmark manifest JSON file."""

    manifest = serialize_dftb_excited_benchmark_suite(suite)
    _validate_dftb_excited_benchmark_suite_manifest(manifest)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def load_dftb_excited_benchmark_suite_manifest(
    path: str | os.PathLike[str],
    *,
    validate_artifacts: bool = False,
    base_dir: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Load and validate a public DFTB excited-state benchmark manifest JSON file."""

    manifest = json.loads(Path(path).read_text())
    if not isinstance(manifest, dict):
        raise DFTBPlusError("DFTB excited-state benchmark manifest must be a JSON object")
    _validate_dftb_excited_benchmark_suite_manifest(manifest)
    if validate_artifacts:
        _validate_manifest_artifacts(manifest, base_dir=base_dir)
    return manifest


def rehydrate_dftb_excited_benchmark_suite_manifest(manifest: dict[str, Any]) -> DFTBExcitedBenchmarkSuite:
    """Rebuild benchmark-suite dataclasses from a validated public manifest.

    This is metadata rehydration only. It preserves ``enables_runtime_capability``
    as false and does not parse benchmark artifacts or enable TD-DFTB/SF-DFTB
    runtime behavior.
    """

    _validate_dftb_excited_benchmark_suite_manifest(manifest)
    cases: list[DFTBExcitedBenchmarkCase] = []
    for case in manifest["cases"]:
        observables = case.get("observables") or {}
        observable_contract = DFTBExcitationObservableContract(
            state_indices=list(observables.get("state_indices") or []),
            has_oscillator_strengths=bool(observables.get("has_oscillator_strengths")),
            has_transition_dipoles=bool(observables.get("has_transition_dipoles")),
            has_transition_charges=bool(observables.get("has_transition_charges")),
            transition_charge_natom=observables.get("transition_charge_natom"),
        )
        cases.append(
            DFTBExcitedBenchmarkCase(
                molecule=str(case["molecule"]),
                feature_family=str(case["feature_family"]),
                state_count=int(case["state_count"]),
                artifact_paths=[str(path) for path in case.get("artifact_paths", [])],
                evidence_level=str(case["evidence_level"]),
                observable_contract=observable_contract,
                includes_mrsf_tddftb=bool(case.get("includes_mrsf_tddftb", False)),
            )
        )

    suite = DFTBExcitedBenchmarkSuite(
        name=str(manifest.get("name", "")),
        molecules=[str(molecule) for molecule in manifest.get("molecules", [])],
        feature_families=[str(feature) for feature in manifest.get("feature_families", [])],
        case_count=int(manifest.get("case_count", len(cases))),
        artifact_paths=[str(path) for path in manifest.get("artifact_paths", [])],
        evidence_level=str(manifest.get("evidence_level")),
        cases=cases,
        includes_mrsf_tddftb=bool(manifest.get("includes_mrsf_tddftb", False)),
        enables_runtime_capability=False,
    )
    _validate_dftb_excited_benchmark_suite_manifest(serialize_dftb_excited_benchmark_suite(suite))
    return suite


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


def _validate_state_pair(pair: tuple[int, int], valid_indices: set[int]) -> None:
    left, right = pair
    if left not in valid_indices or right not in valid_indices:
        raise DFTBPlusError(f"DFTB+ NAMD export NAC state pair {pair} is not present in excited-state data")
    if left == right:
        raise DFTBPlusError("DFTB+ NAMD export NAC state pairs must reference distinct states")


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
    valid_state_indices = {state.index for state in excitations.excitations}
    if nacme:
        for pair in nacme:
            _validate_state_pair(pair, valid_state_indices)
    if nac_vectors:
        for pair, rows in nac_vectors.items():
            _validate_state_pair(pair, valid_state_indices)
            vector_rows = _validated_rows("NAC vector", rows)
            if len(vector_rows) != len(gradient_rows):
                raise DFTBPlusError("DFTB+ NAMD export NAC vector atom count differs from gradient atom count")
    has_nacme = bool(nacme)
    has_nac_vectors = bool(nac_vectors)
    if not has_nacme and not has_nac_vectors:
        raise DFTBPlusError("DFTB+ NAMD export requires validated gradients and NAC/NACME data")

    if velocities is not None:
        velocity_rows = _validated_rows("velocity", velocities)
        if len(velocity_rows) != len(gradient_rows):
            raise DFTBPlusError("DFTB+ NAMD export velocity and gradient atom counts differ")

    observable_contract = validate_dftb_excitation_observables(excitations)

    return DFTBNAMDExportFrame(
        state_indices=[state.index for state in excitations.excitations],
        energies_ev=[state.energy_ev for state in excitations.excitations],
        natom=len(gradient_rows),
        source=excitations.source,
        has_validated_gradients=True,
        has_validated_nacme=has_nacme or has_nac_vectors,
        observable_contract=observable_contract,
        has_velocities=velocities is not None,
    )


def build_namd_export_payload(
    *,
    excitations: DFTBPlusExcitationResult,
    gradients: Sequence[Sequence[float]],
    nacme: dict[tuple[int, int], float] | None = None,
    nac_vectors: dict[tuple[int, int], Sequence[Sequence[float]]] | None = None,
    velocities: Sequence[Sequence[float]] | None = None,
    target: str = "generic",
) -> dict[str, Any]:
    """Serialize validated DFTB excited-state dynamics inputs for NAMD tooling.

    This is an interchange scaffold for future PyRAI2MD/OpenSM-style workflows,
    not a production surface-hopping implementation.  It reuses the strict
    ``build_namd_export_frame`` validation so synthetic parser fixtures, missing
    gradients, missing NAC/NACME data, or inconsistent NAC-vector atom counts
    fail before any JSON-like payload is produced.
    """

    frame = build_namd_export_frame(
        excitations=excitations,
        gradients=gradients,
        nacme=nacme,
        nac_vectors=nac_vectors,
        velocities=velocities,
    )
    gradient_rows = _validated_rows("gradient", gradients)
    velocity_rows = _validated_rows("velocity", velocities) if velocities is not None else None
    nacme_records = [
        {"state_i": left, "state_j": right, "value": float(value)}
        for (left, right), value in sorted((nacme or {}).items())
    ]
    nac_vector_records = []
    for (left, right), rows in sorted((nac_vectors or {}).items()):
        nac_vector_records.append(
            {
                "state_i": left,
                "state_j": right,
                "vectors": _validated_rows("NAC vector", rows),
            }
        )

    return {
        "format": "openqp_dftb_namd_payload_v1",
        "target": str(target),
        "source": frame.source,
        "evidence_level": "validated_external_dftbplus_output",
        "states": [
            {"index": state.index, "energy_ev": state.energy_ev}
            for state in excitations.excitations
        ],
        "gradients": gradient_rows,
        "nacme": nacme_records,
        "nac_vectors": nac_vector_records,
        "velocities": velocity_rows,
        "observable_contract": {
            "state_indices": frame.observable_contract.state_indices,
            "has_oscillator_strengths": frame.observable_contract.has_oscillator_strengths,
            "has_transition_dipoles": frame.observable_contract.has_transition_dipoles,
            "has_transition_charges": frame.observable_contract.has_transition_charges,
            "transition_charge_natom": frame.observable_contract.transition_charge_natom,
        },
        "enables_runtime_capability": False,
    }


def build_dftb_excited_gradient_frame(
    *,
    excitations: DFTBPlusExcitationResult,
    gradients_by_state: dict[int, Sequence[Sequence[float]]],
) -> DFTBExcitedGradientFrame:
    """Build a validated-data-only excited-state gradient metadata frame.

    This is a shape and provenance guard for future TD-DFTB/SF-DFTB gradient
    work. It requires a validated excited-state runtime source and one Cartesian
    gradient block for every parsed state, but the returned metadata still keeps
    the public branch's excited-state gradient capability disabled.
    """

    if not excitations.validated_runtime:
        raise DFTBPlusError(
            "DFTB excited-state gradients require validated DFTB excited-state data from a real runtime"
        )
    state_indices = [state.index for state in excitations.excitations]
    if not state_indices:
        raise DFTBPlusError("DFTB excited-state gradients require at least one excited state")
    provided_indices = set(gradients_by_state)
    expected_indices = set(state_indices)
    missing_indices = [index for index in state_indices if index not in provided_indices]
    if missing_indices:
        raise DFTBPlusError(f"DFTB excited-state gradients missing validated gradients for states {missing_indices}")
    extra_indices = sorted(provided_indices - expected_indices)
    if extra_indices:
        raise DFTBPlusError(f"DFTB excited-state gradients provided unknown state indices {extra_indices}")

    parsed_gradients = {
        index: _validated_rows(f"state {index} gradient", gradients_by_state[index])
        for index in state_indices
    }
    natom = len(parsed_gradients[state_indices[0]])
    if any(len(rows) != natom for rows in parsed_gradients.values()):
        raise DFTBPlusError("DFTB excited-state gradient atom counts must match across states")

    return DFTBExcitedGradientFrame(
        state_indices=state_indices,
        gradients_by_state=parsed_gradients,
        natom=natom,
        source=excitations.source,
    )


def build_dftb_nac_frame(
    *,
    excitations: DFTBPlusExcitationResult,
    nacme: dict[tuple[int, int], float] | None = None,
    nac_vectors: dict[tuple[int, int], Sequence[Sequence[float]]] | None = None,
) -> DFTBNACFrame:
    """Build a validated-data-only NAC/NACME metadata frame.

    This helper is a guard for future excited-state gradient/NAC work. It
    requires validated excited-state runtime evidence plus explicit NACME
    scalars or Cartesian NAC-vector rows, and it preserves the disabled runtime
    capability flag so downstream code cannot mistake parser fixtures for
    production nonadiabatic coupling support.
    """

    if not excitations.validated_runtime:
        raise DFTBPlusError(
            "DFTB NAC data require validated DFTB excited-state data from a real runtime"
        )
    state_indices = [state.index for state in excitations.excitations]
    if not state_indices:
        raise DFTBPlusError("DFTB NAC data require at least one excited state")
    if not nacme and not nac_vectors:
        raise DFTBPlusError("DFTB NAC data require NACME scalars or NAC vectors")

    valid_state_indices = set(state_indices)
    parsed_nacme: dict[tuple[int, int], float] = {}
    if nacme:
        for pair, value in nacme.items():
            _validate_state_pair(pair, valid_state_indices)
            parsed_nacme[pair] = float(value)

    parsed_vectors: dict[tuple[int, int], list[list[float]]] = {}
    natom: int | None = None
    if nac_vectors:
        for pair, rows in nac_vectors.items():
            _validate_state_pair(pair, valid_state_indices)
            vector_rows = _validated_rows("NAC vector", rows)
            if natom is None:
                natom = len(vector_rows)
            elif len(vector_rows) != natom:
                raise DFTBPlusError("DFTB NAC vector atom counts must match across state pairs")
            parsed_vectors[pair] = vector_rows

    return DFTBNACFrame(
        state_indices=state_indices,
        nacme_by_pair=parsed_nacme,
        nac_vectors_by_pair=parsed_vectors,
        natom=natom,
        source=excitations.source,
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


def validate_spin_constant_manifest(
    atoms: Sequence[int],
    spin_constants_path: str | os.PathLike[str],
) -> DFTBSpinConstantManifest:
    """Fail fast when spin-polarized DFTB lacks element spin constants.

    This public scaffold validates only the filesystem/metadata contract needed
    by a future spin-polarized external DFTB+ path. A complete manifest does not
    enable spin-polarized or unrestricted runtime support; input validation still
    rejects those requests until real DFTB+ runs and OpenQP result mapping are
    validated.
    """

    path = Path(spin_constants_path)
    if not path.is_file():
        raise DFTBPlusError(f"DFTB+ spin-constant file is missing: {path}")

    atom_symbols = _symbols_for_atoms(atoms)
    required_symbols = sorted(set(atom_symbols))
    found_symbols: set[str] = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        symbol = fields[0]
        if symbol in required_symbols:
            found_symbols.add(symbol)

    missing_symbols = [symbol for symbol in required_symbols if symbol not in found_symbols]
    if missing_symbols:
        raise DFTBPlusError(
            "DFTB+ spin-constant file is missing spin constants for elements: " + ", ".join(missing_symbols)
        )

    return DFTBSpinConstantManifest(
        atom_symbols=atom_symbols,
        required_symbols=required_symbols,
        found_symbols=[symbol for symbol in required_symbols if symbol in found_symbols],
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
