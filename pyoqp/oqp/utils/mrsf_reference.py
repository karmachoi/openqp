"""Prototype helpers for ensemble-reference MRSF.

This module defines the shared input parsing, SCF metadata, and lightweight
block-response bookkeeping that make the intended mixed-reference ROHF ensemble
explicit before the native response kernel is generalized.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math
import numpy as np
import re
from typing import Any


EV_PER_HARTREE = 27.211386245988
MRSF_REFERENCE_CANONICAL_MODES = {"off", "diagnostic", "ensemble"}
MRSF_REFERENCE_ALIASES = {"state_average": "ensemble"}
MRSF_REFERENCE_MODES = MRSF_REFERENCE_CANONICAL_MODES | set(MRSF_REFERENCE_ALIASES)
MRSF_REFERENCE_TRIAL_VECTOR_MODES = {"adaptive", "native"}
DEFAULT_MAX_REFS = 6
DEFAULT_TRIAL_SHIFT = 1.0e6


class MrsfReferenceError(ValueError):
    """Invalid ensemble-reference MRSF configuration."""


@dataclass(frozen=True)
class ParsedMrsfReferenceConfig:
    mode: str
    open_pairs: list[tuple[int, int]]
    pair_mode: str
    weight_mode: str
    weights: list[float]
    weight_temperature: float
    max_refs: int
    gap_threshold: float
    overlap_threshold: float
    trial_vectors: str
    trial_shift: float
    strict: bool


def parse_reference_pairs(raw: Any) -> list[tuple[int, int]]:
    """Parse one-based candidate open-shell MO pairs.

    Accepted forms include ``"12:13; 11:14"`` and ``"12,13;11,14"``.
    ``auto`` or an empty value means the current ROHF open pair will be inferred
    after SCF when electron counts are available.
    """

    if raw is None:
        return []
    if isinstance(raw, (list, tuple)):
        pairs = []
        for item in raw:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                pairs.append(_validate_pair(int(item[0]), int(item[1])))
            else:
                pairs.extend(parse_reference_pairs(item))
        return pairs

    text = str(raw).strip()
    if not text or text.lower() == "auto":
        return []

    pairs = []
    for chunk in re.split(r"[;|]", text):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = [part for part in re.split(r"[:,\s]+", chunk) if part]
        if len(parts) != 2:
            raise MrsfReferenceError(
                "mrsf_ref.open_pairs entries must contain exactly two one-based MO indices"
            )
        try:
            left, right = int(parts[0]), int(parts[1])
        except ValueError as exc:
            raise MrsfReferenceError("mrsf_ref.open_pairs must contain integer MO indices") from exc
        pairs.append(_validate_pair(left, right))
    return pairs


def parse_weights(raw: Any, nrefs: int) -> list[float]:
    """Parse ensemble weights for ``nrefs`` references.

    ``equal``, ``auto``, and empty values produce equal weights.  Explicit
    weights must be non-negative and sum to one, avoiding silent changes to the
    intended ensemble.
    """

    if nrefs <= 0:
        return []
    if raw is None:
        return _equal_weights(nrefs)

    if isinstance(raw, (list, tuple)):
        values = [float(item) for item in raw]
    else:
        text = str(raw).strip()
        if not text or text.lower() in {"equal", "auto"}:
            return _equal_weights(nrefs)
        try:
            values = [float(item) for item in re.split(r"[,;\s]+", text) if item]
        except ValueError as exc:
            raise MrsfReferenceError("mrsf_ref.weights must be numeric or equal") from exc

    if len(values) != nrefs:
        raise MrsfReferenceError(
            f"mrsf_ref.weights has {len(values)} values but {nrefs} reference(s) are expected"
        )
    if any(value < 0.0 for value in values):
        raise MrsfReferenceError("mrsf_ref.weights must be non-negative")
    total = sum(values)
    if total <= 0.0:
        raise MrsfReferenceError("mrsf_ref.weights must have a positive sum")
    if not math.isclose(total, 1.0, rel_tol=1.0e-10, abs_tol=1.0e-10):
        raise MrsfReferenceError("mrsf_ref.weights must sum to 1.0")
    return values


def parse_weight_request(raw: Any, nrefs: int) -> tuple[str, list[float]]:
    """Parse a weight request into a mode and optional explicit weights."""

    if raw is None:
        return "equal", _equal_weights(nrefs)
    if isinstance(raw, (list, tuple)):
        return "manual", parse_weights(raw, nrefs)

    text = str(raw).strip().lower()
    if not text or text in {"equal", "auto"}:
        return "equal", _equal_weights(nrefs)
    if text in {"gap_softmax", "softmax", "gap"}:
        return "gap_softmax", _equal_weights(nrefs)
    return "manual", parse_weights(raw, nrefs)


def parse_mrsf_reference_config(config: dict[str, Any]) -> ParsedMrsfReferenceConfig:
    section = config.get("mrsf_ref", {}) if isinstance(config, dict) else {}
    if not isinstance(section, dict):
        raise MrsfReferenceError("[mrsf_ref] must be a mapping")

    raw_mode = str(section.get("mode", "off")).strip().lower()
    mode = MRSF_REFERENCE_ALIASES.get(raw_mode, raw_mode)
    if mode not in MRSF_REFERENCE_CANONICAL_MODES:
        raise MrsfReferenceError(
            f"mrsf_ref.mode must be one of {', '.join(sorted(MRSF_REFERENCE_CANONICAL_MODES))}"
        )

    try:
        max_refs = int(section.get("max_refs", DEFAULT_MAX_REFS))
    except (TypeError, ValueError) as exc:
        raise MrsfReferenceError("mrsf_ref.max_refs must be an integer") from exc
    if max_refs < 1:
        raise MrsfReferenceError("mrsf_ref.max_refs must be at least 1")

    try:
        gap_threshold = float(section.get("gap_threshold", 0.01))
    except (TypeError, ValueError) as exc:
        raise MrsfReferenceError("mrsf_ref.gap_threshold must be numeric") from exc
    if gap_threshold <= 0.0:
        raise MrsfReferenceError("mrsf_ref.gap_threshold must be positive")

    try:
        overlap_threshold = float(section.get("overlap_threshold", 0.85))
    except (TypeError, ValueError) as exc:
        raise MrsfReferenceError("mrsf_ref.overlap_threshold must be numeric") from exc
    if not 0.0 <= overlap_threshold <= 1.0:
        raise MrsfReferenceError("mrsf_ref.overlap_threshold must be between 0 and 1")

    try:
        weight_temperature = float(section.get("weight_temperature", 0.05))
    except (TypeError, ValueError) as exc:
        raise MrsfReferenceError("mrsf_ref.weight_temperature must be numeric") from exc
    if weight_temperature <= 0.0:
        raise MrsfReferenceError("mrsf_ref.weight_temperature must be positive")

    trial_vectors = str(section.get("trial_vectors", "adaptive")).strip().lower()
    if trial_vectors not in MRSF_REFERENCE_TRIAL_VECTOR_MODES:
        raise MrsfReferenceError(
            f"mrsf_ref.trial_vectors must be one of {', '.join(sorted(MRSF_REFERENCE_TRIAL_VECTOR_MODES))}"
        )

    try:
        trial_shift = float(section.get("trial_shift", DEFAULT_TRIAL_SHIFT))
    except (TypeError, ValueError) as exc:
        raise MrsfReferenceError("mrsf_ref.trial_shift must be numeric") from exc
    if trial_shift <= 0.0:
        raise MrsfReferenceError("mrsf_ref.trial_shift must be positive")

    open_pairs = parse_reference_pairs(section.get("open_pairs", "auto"))
    pair_mode = "manual" if open_pairs else "auto"
    nrefs_for_weights = len(open_pairs) if open_pairs else max_refs
    weight_mode, weights = parse_weight_request(section.get("weights", "equal"), nrefs_for_weights)

    return ParsedMrsfReferenceConfig(
        mode=mode,
        open_pairs=open_pairs,
        pair_mode=pair_mode,
        weight_mode=weight_mode,
        weights=weights,
        weight_temperature=weight_temperature,
        max_refs=max_refs,
        gap_threshold=gap_threshold,
        overlap_threshold=overlap_threshold,
        trial_vectors=trial_vectors,
        trial_shift=trial_shift,
        strict=_parse_bool(section.get("strict", False)),
    )


def build_mrsf_reference_metadata(config: dict[str, Any], data: Any = None) -> dict[str, Any]:
    """Build a JSON-safe MRSF ensemble-reference metadata block."""

    parsed = parse_mrsf_reference_config(config)
    frontier = _frontier_metadata(data, parsed.gap_threshold)
    refs, pair_selection = _select_reference_pairs(parsed, frontier, data)

    weights, weight_model = _resolve_reference_weights(parsed, refs, data)

    ensemble = build_reference_ensemble(refs, weights, data)

    warnings = []
    if parsed.mode == "ensemble":
        warnings.append(
            "ensemble mixed-reference MRSF uses an energy-only state-interaction prototype; the full coupled ensemble-response kernel is not implemented"
        )
    if parsed.mode != "off" and not refs:
        warnings.append(
            "no explicit open_pairs and the current ROHF open-shell pair is not available yet"
        )
    if parsed.mode == "ensemble" and parsed.pair_mode == "auto" and len(refs) < 2:
        warnings.append(
            "automatic open-pair selection found fewer than two candidate ROHF configurations"
        )
    if weight_model.get("fallback") == "equal":
        warnings.append(
            f"mrsf_ref.weights={parsed.weight_mode} fell back to equal weights: {weight_model.get('reason')}"
        )
    if frontier.get("ambiguous"):
        warnings.append(
            "frontier MO gap is below mrsf_ref.gap_threshold; single-reference MRSF may switch character"
        )
    if parsed.mode != "off" and ensemble.get("status", "").startswith("invalid"):
        warnings.append(str(ensemble.get("reason", "invalid MRSF reference ensemble")))

    return {
        "status": _status(parsed.mode),
        "mode": parsed.mode,
        "implemented": parsed.mode in {"off", "diagnostic", "ensemble"},
        "scf_implemented": parsed.mode in {"off", "diagnostic", "ensemble"},
        "response_implemented": parsed.mode in {"off", "diagnostic", "ensemble"},
        "theory": {
            "reference_model": "mixed_rohf_triplet_reference_ensemble",
            "mean_field_target": "fractional occupation ensemble over candidate ROHF configurations",
            "response_model": "state-interaction MRSF response over reference-specific spin-flip block states",
            "target_response_model": "block-coupled MRSF response over all reference-specific spin-flip spaces",
            "coupled_response_required": parsed.mode == "ensemble",
            "inter_reference_coupling": False,
            "energy_only": parsed.mode == "ensemble",
        },
        "pair_selection": pair_selection,
        "weight_model": weight_model,
        "open_pairs": [[int(i), int(j)] for i, j in refs],
        "weights": [float(w) for w in weights],
        "weight_temperature_hartree": parsed.weight_temperature,
        "weight_temperature_ev": parsed.weight_temperature * EV_PER_HARTREE,
        "max_refs": parsed.max_refs,
        "gap_threshold_hartree": parsed.gap_threshold,
        "gap_threshold_ev": parsed.gap_threshold * EV_PER_HARTREE,
        "overlap_threshold": parsed.overlap_threshold,
        "trial_vector_model": {
            "mode": parsed.trial_vectors,
            "active_virtual_shift_hartree": parsed.trial_shift,
            "purpose": "avoid seeding active-space reference-changing intruders in the uncoupled block prototype",
        },
        "strict": parsed.strict,
        "frontier": frontier,
        "ensemble": ensemble,
        "scf": {
            "ensemble_occupations_applied": False,
            "occupation_tags": [],
        },
        "warnings": warnings,
    }


def build_reference_ensemble(
    open_pairs: list[tuple[int, int]],
    weights: list[float],
    data: Any = None,
) -> dict[str, Any]:
    """Build the mixed-reference ROHF occupation model for candidate references.

    Each candidate is a triplet ROHF determinant with the requested open-shell
    pair occupied by alpha electrons and the lowest available remaining orbitals
    doubly occupied.  The weighted sum of those determinants is the mean-field
    target for a discontinuity-free ensemble reference.
    """

    if not open_pairs:
        return {
            "available": False,
            "status": "no_references",
            "reason": "no candidate open-shell pairs are available",
        }

    nelec_alpha = _safe_int(_data_get(data, "nelec_A"))
    nelec_beta = _safe_int(_data_get(data, "nelec_B"))
    nmo = _infer_nmo(data, open_pairs)
    if nelec_alpha is None or nelec_beta is None:
        return {
            "available": False,
            "status": "missing_electron_counts",
            "reason": "electron counts are not available",
        }
    if nelec_alpha - nelec_beta != 2:
        return {
            "available": False,
            "status": "invalid_electron_count",
            "reason": "ensemble-reference MRSF currently requires a triplet ROHF pair with nelec_A = nelec_B + 2",
            "nelec_alpha": nelec_alpha,
            "nelec_beta": nelec_beta,
        }

    normalized_weights = weights if len(weights) == len(open_pairs) else _equal_weights(len(open_pairs))
    references = []
    alpha_ensemble = [0.0 for _ in range(nmo)]
    beta_ensemble = [0.0 for _ in range(nmo)]
    response_dim_total = 0
    response_dim_triplet_total = 0

    for ref_id, (pair, weight) in enumerate(zip(open_pairs, normalized_weights), start=1):
        reference = _build_rohf_reference(ref_id, pair, weight, nelec_alpha, nelec_beta, nmo)
        if not reference.get("valid", False):
            return {
                "available": False,
                "status": "invalid_reference",
                "reason": reference.get("reason", "invalid reference"),
                "reference": reference,
            }

        references.append(reference)
        response_dim_total += int(reference["response_space"]["raw_dimension"])
        response_dim_triplet_total += int(reference["response_space"]["triplet_dimension"])
        for mo in reference["alpha_occupied"]:
            alpha_ensemble[mo - 1] += float(weight)
        for mo in reference["beta_occupied"]:
            beta_ensemble[mo - 1] += float(weight)

    active_open_orbitals = sorted({mo for pair in open_pairs for mo in pair})
    alpha_nonzero = _nonzero_occupations(alpha_ensemble)
    beta_nonzero = _nonzero_occupations(beta_ensemble)

    return {
        "available": True,
        "status": "ready_for_coupled_response",
        "n_references": len(references),
        "nmo": nmo,
        "nelec_alpha": nelec_alpha,
        "nelec_beta": nelec_beta,
        "closed_shell_rule": "lowest MO indices not present in the candidate open-shell pair",
        "active_open_orbitals": active_open_orbitals,
        "references": references,
        "ensemble_occupations": {
            "alpha": alpha_nonzero,
            "beta": beta_nonzero,
            "alpha_sum": float(sum(alpha_ensemble)),
            "beta_sum": float(sum(beta_ensemble)),
        },
        "response_space": {
            "block_count": len(references),
            "raw_dimension": response_dim_total,
            "triplet_dimension": response_dim_triplet_total,
            "coupling": "state_interaction_between_reference_blocks",
        },
    }


def requires_coupled_response(config: dict[str, Any]) -> bool:
    """Return true when the requested mode needs the unfinished solver."""

    return parse_mrsf_reference_config(config).mode == "ensemble"


def ensemble_occupation_vectors(metadata: dict[str, Any]) -> tuple[list[float], list[float]]:
    """Return dense alpha and beta occupation vectors from metadata.

    The vectors use spin occupations, not ROHF closed/open labels: a closed MO
    has alpha=1 and beta=1, while ensemble-active MOs can be fractional.
    """

    ensemble = metadata.get("ensemble", {}) if isinstance(metadata, dict) else {}
    if not ensemble.get("available", False):
        raise MrsfReferenceError(str(ensemble.get("reason", "MRSF reference ensemble is not available")))

    nmo = _safe_int(ensemble.get("nmo"))
    if nmo is None or nmo <= 0:
        raise MrsfReferenceError("MRSF reference ensemble does not define a positive MO count")

    occupations = ensemble.get("ensemble_occupations", {})
    alpha = _dense_occupation_vector(occupations.get("alpha", []), nmo, "alpha")
    beta = _dense_occupation_vector(occupations.get("beta", []), nmo, "beta")

    nelec_alpha = float(ensemble.get("nelec_alpha", 0.0))
    nelec_beta = float(ensemble.get("nelec_beta", 0.0))
    if not math.isclose(sum(alpha), nelec_alpha, rel_tol=1.0e-10, abs_tol=1.0e-10):
        raise MrsfReferenceError("MRSF reference alpha occupations do not sum to nelec_A")
    if not math.isclose(sum(beta), nelec_beta, rel_tol=1.0e-10, abs_tol=1.0e-10):
        raise MrsfReferenceError("MRSF reference beta occupations do not sum to nelec_B")

    return alpha, beta


def reference_mo_permutation(reference: dict[str, Any], nmo: int) -> list[int]:
    """Return the one-based MO order expected by the native single-reference kernel.

    Native MRSF assumes the ROHF open pair is the last two alpha-occupied
    orbitals: closed orbitals first, then the two open orbitals, then the
    remaining virtual orbitals.  The ensemble references keep their original MO
    labels in metadata, so each response block uses this permutation before
    entering the native solver.
    """

    nmo_int = _safe_int(nmo)
    if nmo_int is None or nmo_int <= 0:
        raise MrsfReferenceError("reference MO permutation requires a positive nmo")

    try:
        closed = [int(item) for item in reference.get("closed_orbitals", [])]
        open_pair = [int(item) for item in reference.get("open_pair", [])]
    except (TypeError, ValueError) as exc:
        raise MrsfReferenceError("reference contains non-integer orbital labels") from exc

    if len(open_pair) != 2:
        raise MrsfReferenceError("reference MO permutation requires exactly two open orbitals")

    ordered = closed + open_pair
    if any(item < 1 or item > nmo_int for item in ordered):
        raise MrsfReferenceError("reference MO permutation contains out-of-range orbital labels")
    if len(set(ordered)) != len(ordered):
        raise MrsfReferenceError("reference MO permutation contains duplicate occupied orbitals")

    occupied = set(ordered)
    permutation = ordered + [mo for mo in range(1, nmo_int + 1) if mo not in occupied]
    if len(permutation) != nmo_int or len(set(permutation)) != nmo_int:
        raise MrsfReferenceError("reference MO permutation is not a complete MO ordering")
    return permutation


def collect_block_diagonal_response(
    blocks: list[dict[str, Any]],
    nstate: int,
) -> dict[str, Any]:
    """Sort uncoupled per-reference MRSF energies into a combined state list."""

    requested = _safe_int(nstate)
    if requested is None or requested < 1:
        raise MrsfReferenceError("block-diagonal response collection requires nstate >= 1")

    candidates: list[dict[str, Any]] = []
    raw_candidate_count = 0
    skipped_nonconverged_blocks = []
    for block_index, block in enumerate(blocks, start=1):
        energies = block.get("energies", [])
        reference_id = int(block.get("reference_id", block_index))
        raw_candidate_count += len(energies)
        if not bool(block.get("converged", True)):
            skipped_nonconverged_blocks.append(reference_id)
            continue
        for state_index, energy in enumerate(energies, start=1):
            value = float(energy)
            if not math.isfinite(value):
                continue
            candidates.append(
                {
                    "energy": value,
                    "reference_id": reference_id,
                    "block_index": block_index,
                    "state_index": state_index,
                    "weight": float(block.get("weight", 0.0)),
                    "open_pair": [int(item) for item in block.get("open_pair", [])],
                }
            )

    candidates.sort(key=lambda item: (item["energy"], item["reference_id"], item["state_index"]))
    selected = candidates[:requested]
    for rank, item in enumerate(selected, start=1):
        item["rank"] = rank

    return {
        "status": "ready" if selected else "no_states",
        "model": "block_diagonal_uncoupled",
        "block_count": len(blocks),
        "candidate_count": len(candidates),
        "raw_candidate_count": raw_candidate_count,
        "skipped_nonconverged_blocks": skipped_nonconverged_blocks,
        "requested_states": requested,
        "energies": [float(item["energy"]) for item in selected],
        "selected_states": selected,
    }


def mrsf_response_labels(reference: dict[str, Any], nmo: int) -> list[tuple[int, int]]:
    """Return common-basis labels for a reference-local MRSF vector.

    Native MRSF stores vectors in the local permuted order with beta-virtual
    labels running slowest and alpha-occupied labels fastest.  This helper maps
    those positions back to original MO labels so vectors from different
    references can be compared in one response basis.
    """

    permutation = reference_mo_permutation(reference, nmo)
    nalpha = len(reference.get("closed_orbitals", [])) + 2
    nbeta = len(reference.get("closed_orbitals", []))
    if nalpha <= 0 or nbeta < 0 or nalpha > len(permutation):
        raise MrsfReferenceError("reference does not define a valid MRSF occupation partition")

    labels = []
    for particle_position in range(nbeta + 1, len(permutation) + 1):
        particle_label = int(permutation[particle_position - 1])
        for hole_position in range(1, nalpha + 1):
            hole_label = int(permutation[hole_position - 1])
            labels.append((hole_label, particle_label))
    return labels


def collect_state_interaction_response(
    blocks: list[dict[str, Any]],
    block_vectors: dict[int, Any],
    references: list[dict[str, Any]],
    nstate: int,
    nmo: int,
    metric_threshold: float = 1.0e-8,
) -> dict[str, Any]:
    """Build a conservative state-interaction response from block MRSF states.

    The diagonal energies come from converged native MRSF blocks.  Off-diagonal
    Hamiltonian elements use the common-basis vector overlap,
    ``H_ij = 0.5 * S_ij * (E_i + E_j)``.  This is an energy-only state
    interaction approximation, not the full off-diagonal response kernel.
    """

    requested = _safe_int(nstate)
    nmo_int = _safe_int(nmo)
    if requested is None or requested < 1:
        raise MrsfReferenceError("state-interaction response collection requires nstate >= 1")
    if nmo_int is None or nmo_int <= 0:
        raise MrsfReferenceError("state-interaction response collection requires a positive nmo")

    reference_by_id = {
        int(reference.get("id", index)): reference
        for index, reference in enumerate(references, start=1)
    }

    candidates = []
    skipped_nonconverged_blocks = []
    skipped_vector_blocks = []
    common_labels: set[tuple[int, int]] = set()

    for block_index, block in enumerate(blocks, start=1):
        reference_id = int(block.get("reference_id", block_index))
        if not bool(block.get("converged", True)):
            skipped_nonconverged_blocks.append(reference_id)
            continue

        reference = reference_by_id.get(reference_id)
        vectors = block_vectors.get(reference_id)
        if reference is None or vectors is None:
            skipped_vector_blocks.append(reference_id)
            continue

        vector_array = np.asarray(vectors, dtype=float)
        if vector_array.ndim == 1:
            vector_array = vector_array.reshape((-1, 1))

        try:
            labels = mrsf_response_labels(reference, nmo_int)
        except MrsfReferenceError:
            skipped_vector_blocks.append(reference_id)
            continue

        if vector_array.shape[0] != len(labels):
            skipped_vector_blocks.append(reference_id)
            continue

        for state_index, energy in enumerate(block.get("energies", []), start=1):
            if state_index > vector_array.shape[1]:
                continue
            value = float(energy)
            if not math.isfinite(value):
                continue

            coeffs: dict[tuple[int, int], float] = {}
            for label, amplitude in zip(labels, vector_array[:, state_index - 1]):
                amp = float(amplitude)
                if amp == 0.0:
                    continue
                coeffs[label] = coeffs.get(label, 0.0) + amp
            norm = math.sqrt(sum(amp * amp for amp in coeffs.values()))
            if norm <= 0.0 or not math.isfinite(norm):
                continue
            coeffs = {label: amp / norm for label, amp in coeffs.items()}
            common_labels.update(coeffs)
            candidates.append(
                {
                    "energy": value,
                    "reference_id": reference_id,
                    "block_index": block_index,
                    "state_index": state_index,
                    "weight": float(block.get("weight", 0.0)),
                    "open_pair": [int(item) for item in block.get("open_pair", [])],
                    "coefficients": coeffs,
                }
            )

    if not candidates:
        return {
            "status": "no_states",
            "model": "state_interaction_overlap",
            "candidate_count": 0,
            "requested_states": requested,
            "energies": [],
            "selected_states": [],
            "skipped_nonconverged_blocks": skipped_nonconverged_blocks,
            "skipped_vector_blocks": skipped_vector_blocks,
        }

    label_list = sorted(common_labels)
    label_index = {label: index for index, label in enumerate(label_list)}
    vectors = np.zeros((len(candidates), len(label_list)), dtype=float)
    energies = np.asarray([candidate["energy"] for candidate in candidates], dtype=float)
    for candidate_index, candidate in enumerate(candidates):
        for label, amplitude in candidate["coefficients"].items():
            vectors[candidate_index, label_index[label]] = amplitude

    overlap = vectors @ vectors.T
    hamiltonian = 0.5 * overlap * (energies[:, None] + energies[None, :])
    np.fill_diagonal(hamiltonian, energies)

    try:
        metric_values, metric_vectors = np.linalg.eigh(overlap)
        keep = metric_values > metric_threshold
        if not np.any(keep):
            raise np.linalg.LinAlgError("state-interaction metric has no positive subspace")
        transform = metric_vectors[:, keep] / np.sqrt(metric_values[keep])
        orthogonal_hamiltonian = transform.T @ hamiltonian @ transform
        roots, root_vectors = np.linalg.eigh(orthogonal_hamiltonian)
        coefficients = transform @ root_vectors
    except np.linalg.LinAlgError as exc:
        return {
            "status": "failed",
            "model": "state_interaction_overlap",
            "reason": str(exc),
            "candidate_count": len(candidates),
            "requested_states": requested,
            "energies": [],
            "selected_states": [],
            "overlap_matrix": overlap.tolist(),
            "hamiltonian_matrix": hamiltonian.tolist(),
            "skipped_nonconverged_blocks": skipped_nonconverged_blocks,
            "skipped_vector_blocks": skipped_vector_blocks,
        }

    order = np.argsort(roots)
    selected = []
    for rank, root_index in enumerate(order[:requested], start=1):
        coeff_column = coefficients[:, root_index]
        components = []
        for candidate, coefficient in zip(candidates, coeff_column):
            components.append(
                {
                    "reference_id": int(candidate["reference_id"]),
                    "block_index": int(candidate["block_index"]),
                    "state_index": int(candidate["state_index"]),
                    "open_pair": [int(item) for item in candidate["open_pair"]],
                    "coefficient": float(coefficient),
                    "abs_coefficient": float(abs(coefficient)),
                    "energy": float(candidate["energy"]),
                }
            )
        components.sort(key=lambda item: item["abs_coefficient"], reverse=True)
        dominant = components[0] if components else {}
        selected.append(
            {
                "energy": float(roots[root_index]),
                "rank": rank,
                "dominant_reference_id": dominant.get("reference_id"),
                "dominant_open_pair": dominant.get("open_pair", []),
                "components": components[: min(6, len(components))],
            }
        )

    return {
        "status": "ready",
        "model": "state_interaction_overlap",
        "coupling": "overlap_averaged_energy",
        "full_response_kernel": False,
        "metric_threshold": float(metric_threshold),
        "candidate_count": len(candidates),
        "common_dimension": len(label_list),
        "requested_states": requested,
        "energies": [float(item["energy"]) for item in selected],
        "selected_states": selected,
        "overlap_matrix": overlap.tolist(),
        "hamiltonian_matrix": hamiltonian.tolist(),
        "skipped_nonconverged_blocks": skipped_nonconverged_blocks,
        "skipped_vector_blocks": skipped_vector_blocks,
    }


def _validate_pair(left: int, right: int) -> tuple[int, int]:
    if left < 1 or right < 1:
        raise MrsfReferenceError("mrsf_ref.open_pairs uses one-based positive MO indices")
    if left == right:
        raise MrsfReferenceError("mrsf_ref.open_pairs cannot repeat the same MO index")
    return tuple(sorted((left, right)))


def _equal_weights(nrefs: int) -> list[float]:
    return [1.0 / nrefs for _ in range(nrefs)]


def _resolve_reference_weights(
    parsed: ParsedMrsfReferenceConfig,
    refs: list[tuple[int, int]],
    data: Any,
) -> tuple[list[float], dict[str, Any]]:
    if not refs:
        return parsed.weights, {
            "mode": parsed.weight_mode,
            "resolved": False,
            "reason": "no reference pairs are available",
            "temperature_hartree": parsed.weight_temperature,
            "temperature_ev": parsed.weight_temperature * EV_PER_HARTREE,
        }

    if parsed.weight_mode in {"equal", "manual"}:
        weights = parsed.weights if len(parsed.weights) == len(refs) else _equal_weights(len(refs))
        return weights, {
            "mode": parsed.weight_mode,
            "resolved": True,
            "source": "explicit" if parsed.weight_mode == "manual" else "equal",
            "temperature_hartree": parsed.weight_temperature,
            "temperature_ev": parsed.weight_temperature * EV_PER_HARTREE,
        }

    if parsed.weight_mode == "gap_softmax":
        weights, scores, fallback_reason = _gap_softmax_weights(refs, data, parsed.weight_temperature)
        model = {
            "mode": "gap_softmax",
            "resolved": fallback_reason is None,
            "source": "orbital_energy_proxy",
            "temperature_hartree": parsed.weight_temperature,
            "temperature_ev": parsed.weight_temperature * EV_PER_HARTREE,
            "scores": [
                {
                    "pair": [int(pair[0]), int(pair[1])],
                    "energy_proxy": None if score is None else float(score),
                }
                for pair, score in zip(refs, scores)
            ],
        }
        if fallback_reason is not None:
            model["fallback"] = "equal"
            model["reason"] = fallback_reason
        return weights, model

    raise MrsfReferenceError(f"unsupported mrsf_ref.weights mode: {parsed.weight_mode}")


def _gap_softmax_weights(
    refs: list[tuple[int, int]],
    data: Any,
    temperature: float,
) -> tuple[list[float], list[float | None], str | None]:
    scores = _reference_weight_scores(refs, data)
    if len(scores) != len(refs) or any(score is None for score in scores):
        return _equal_weights(len(refs)), scores, "reference energy proxies are unavailable"

    min_score = min(float(score) for score in scores if score is not None)
    factors = [math.exp(-(float(score) - min_score) / temperature) for score in scores]
    total = sum(factors)
    if total <= 0.0 or not math.isfinite(total):
        return _equal_weights(len(refs)), scores, "softmax normalization failed"
    return [float(item / total) for item in factors], scores, None


def _reference_weight_scores(refs: list[tuple[int, int]], data: Any) -> list[float | None]:
    nelec_alpha = _safe_int(_data_get(data, "nelec_A"))
    nelec_beta = _safe_int(_data_get(data, "nelec_B"))
    if nelec_alpha is None or nelec_beta is None:
        return [None for _ in refs]
    nmo = _infer_nmo(data, refs)
    energies = _reference_energy_list(data)
    return [
        _reference_energy_proxy(pair, energies, nelec_alpha, nelec_beta, nmo)
        for pair in refs
    ]


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"true", ".true.", "t", "1", "yes", "on"}
    return False


def _status(mode: str) -> str:
    if mode == "off":
        return "disabled"
    if mode == "diagnostic":
        return "diagnostic"
    return "ensemble_requested"


def _select_reference_pairs(
    parsed: ParsedMrsfReferenceConfig,
    frontier: dict[str, Any],
    data: Any,
) -> tuple[list[tuple[int, int]], dict[str, Any]]:
    if parsed.mode == "off":
        return [], {
            "mode": "off",
            "strategy": "disabled",
            "selected_pairs": [],
            "candidate_pairs": [],
            "truncated": False,
        }

    if parsed.open_pairs:
        refs = list(parsed.open_pairs)
        return refs, {
            "mode": "manual",
            "strategy": "explicit_open_pairs",
            "selected_pairs": _pair_records(refs),
            "candidate_pairs": _pair_records(refs),
            "truncated": False,
        }

    if parsed.mode == "ensemble":
        return _auto_reference_pairs(frontier, data, parsed.max_refs, parsed.gap_threshold)

    current_pair = _current_open_pair(frontier)
    refs = [current_pair] if current_pair else []
    return refs, {
        "mode": "auto",
        "strategy": "current_rohf_open_pair",
        "selected_pairs": _pair_records(refs),
        "candidate_pairs": _pair_records(refs),
        "truncated": False,
    }


def _auto_reference_pairs(
    frontier: dict[str, Any],
    data: Any,
    max_refs: int,
    gap_threshold: float,
) -> tuple[list[tuple[int, int]], dict[str, Any]]:
    current_pair = _current_open_pair(frontier)
    if current_pair is None:
        return [], {
            "mode": "auto",
            "strategy": "frontier_window",
            "reason": "current ROHF open pair is unavailable",
            "selected_pairs": [],
            "candidate_pairs": [],
            "truncated": False,
        }

    nelec_alpha = _safe_int(frontier.get("nelec_alpha"))
    nelec_beta = _safe_int(frontier.get("nelec_beta"))
    if nelec_alpha is None or nelec_beta is None or nelec_alpha - nelec_beta != 2:
        return [current_pair], {
            "mode": "auto",
            "strategy": "frontier_window",
            "reason": "automatic state-average selection currently requires a triplet ROHF pair",
            "selected_pairs": _pair_records([current_pair]),
            "candidate_pairs": _pair_records([current_pair]),
            "truncated": False,
        }

    nmo = _infer_nmo(data, [current_pair])
    if nmo <= 0:
        return [current_pair], {
            "mode": "auto",
            "strategy": "frontier_window",
            "reason": "MO space size is unavailable",
            "selected_pairs": _pair_records([current_pair]),
            "candidate_pairs": _pair_records([current_pair]),
            "truncated": False,
        }

    energies = _reference_energy_list(data)
    active_orbitals, active_window = _auto_active_orbitals(
        nelec_alpha,
        nelec_beta,
        nmo,
        max_refs,
        energies,
        gap_threshold,
    )
    current_score = _reference_energy_proxy(current_pair, energies, nelec_alpha, nelec_beta, nmo)
    balanced_pair = _balanced_frontier_pair(nelec_alpha, nelec_beta, nmo)

    candidates = []
    excluded_candidates = []
    for pair in combinations(active_orbitals, 2):
        pair = _validate_pair(pair[0], pair[1])
        score = _reference_energy_proxy(pair, energies, nelec_alpha, nelec_beta, nmo)
        item = {
            "pair": pair,
            "role": _auto_pair_role(pair, current_pair, balanced_pair),
            "energy_proxy": score,
            "delta_energy_proxy": None
            if score is None or current_score is None
            else float(score - current_score),
        }
        if not _auto_pair_is_admissible(pair, current_pair):
            item["excluded_reason"] = "promoted_high_high_pair"
            excluded_candidates.append(item)
            continue
        candidates.append(item)

    candidates.sort(
        key=lambda item: (
            _auto_pair_priority(item["pair"], current_pair, balanced_pair),
            _auto_pair_energy_delta(item["energy_proxy"], current_score),
            _pair_distance(item["pair"], current_pair),
            item["pair"],
        )
    )

    selected = candidates[:max_refs]
    refs = [item["pair"] for item in selected]
    return refs, {
        "mode": "auto",
        "strategy": "frontier_window",
        "active_orbitals": [int(item) for item in active_orbitals],
        "active_window": active_window,
        "selected_pairs": _candidate_records(selected),
        "candidate_pairs": _candidate_records(candidates),
        "excluded_pairs": _candidate_records(excluded_candidates),
        "truncated": len(candidates) > len(selected),
    }


def _current_open_pair(frontier: dict[str, Any]) -> tuple[int, int] | None:
    current = frontier.get("current_open_pair", []) if isinstance(frontier, dict) else []
    if len(current) != 2:
        return None
    try:
        return _validate_pair(int(current[0]), int(current[1]))
    except (TypeError, ValueError, MrsfReferenceError):
        return None


def _auto_active_orbitals(
    nelec_alpha: int,
    nelec_beta: int,
    nmo: int,
    max_refs: int,
    energies: list[float],
    gap_threshold: float,
) -> tuple[list[int], dict[str, Any]]:
    active = {nelec_beta + 1, nelec_beta + 2}
    if len(energies) >= nmo:
        lower_open = nelec_beta + 1
        upper_open = nelec_beta + 2
        lower_energy = energies[lower_open - 1]
        upper_energy = energies[upper_open - 1]
        candidates = []
        for mo in range(1, nmo + 1):
            if mo in active:
                continue
            if mo < lower_open:
                gap = abs(float(energies[mo - 1] - lower_energy))
            elif mo > upper_open:
                gap = abs(float(energies[mo - 1] - upper_energy))
            else:
                continue
            if gap <= gap_threshold:
                candidates.append((gap, abs(mo - lower_open) + abs(mo - upper_open), mo))

        for _gap, _distance, mo in sorted(candidates):
            active.add(mo)
            if _n_pairs(active) >= max_refs and len(active) >= 4:
                break
        return sorted(mo for mo in active if 1 <= mo <= nmo), {
            "mode": "energy_gap",
            "gap_threshold_hartree": float(gap_threshold),
            "candidate_orbitals": [
                {"orbital": int(mo), "gap_hartree": float(gap)}
                for gap, _distance, mo in sorted(candidates)
            ],
        }

    for offset in range(1, nmo + 1):
        left = nelec_beta - offset + 1
        right = nelec_alpha + offset
        if left >= 1:
            active.add(left)
        if right <= nmo:
            active.add(right)
        if _n_pairs(active) >= max_refs and len(active) >= 4:
            break
    return sorted(mo for mo in active if 1 <= mo <= nmo), {
        "mode": "index_fallback",
        "reason": "MO energies unavailable for gap-threshold active-window selection",
    }


def _n_pairs(items: set[int]) -> int:
    nitems = len(items)
    return nitems * (nitems - 1) // 2


def _balanced_frontier_pair(
    nelec_alpha: int,
    nelec_beta: int,
    nmo: int,
) -> tuple[int, int] | None:
    if nelec_beta < 1 or nelec_alpha + 1 > nmo:
        return None
    return _validate_pair(nelec_beta, nelec_alpha + 1)


def _auto_pair_role(
    pair: tuple[int, int],
    current_pair: tuple[int, int],
    balanced_pair: tuple[int, int] | None,
) -> str:
    if pair == current_pair:
        return "current_open_pair"
    if balanced_pair is not None and pair == balanced_pair:
        return "balanced_frontier_pair"
    return "frontier_window_pair"


def _auto_pair_priority(
    pair: tuple[int, int],
    current_pair: tuple[int, int],
    balanced_pair: tuple[int, int] | None,
) -> int:
    if pair == current_pair:
        return 0
    if balanced_pair is not None and pair == balanced_pair:
        return 1
    return 2


def _auto_pair_is_admissible(
    pair: tuple[int, int],
    current_pair: tuple[int, int],
) -> bool:
    """Keep auto references near the current ROHF open-pair manifold.

    The upper-upper pair in a four-orbital frontier window, such as [9, 10]
    for a current [8, 9] triplet, is a promoted non-Aufbau reference in the
    block prototype rather than a competing ROHF open-pair switch.
    """

    return pair[0] <= current_pair[0]


def _auto_pair_energy_delta(score: float | None, current_score: float | None) -> float:
    if score is None or current_score is None:
        return 0.0
    return abs(float(score - current_score))


def _pair_distance(pair: tuple[int, int], current_pair: tuple[int, int]) -> int:
    return abs(pair[0] - current_pair[0]) + abs(pair[1] - current_pair[1])


def _reference_energy_list(data: Any) -> list[float]:
    alpha = _as_float_list(_data_get(data, "OQP::E_MO_A"))
    if alpha:
        return alpha
    return _as_float_list(_data_get(data, "OQP::E_MO_B"))


def _reference_energy_proxy(
    pair: tuple[int, int],
    energies: list[float],
    nelec_alpha: int,
    nelec_beta: int,
    nmo: int,
) -> float | None:
    if len(energies) < nmo:
        return None
    reference = _build_rohf_reference(0, pair, 1.0, nelec_alpha, nelec_beta, nmo)
    if not reference.get("valid", False):
        return None
    alpha_energy = sum(energies[mo - 1] for mo in reference["alpha_occupied"])
    beta_energy = sum(energies[mo - 1] for mo in reference["beta_occupied"])
    return float(alpha_energy + beta_energy)


def _pair_records(pairs: list[tuple[int, int]]) -> list[dict[str, list[int]]]:
    return [{"pair": [int(pair[0]), int(pair[1])]} for pair in pairs]


def _candidate_records(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    records = []
    for item in items:
        pair = item["pair"]
        records.append(
            {
                "pair": [int(pair[0]), int(pair[1])],
                "role": item.get("role", "reference_pair"),
                "energy_proxy": item.get("energy_proxy"),
                "delta_energy_proxy": item.get("delta_energy_proxy"),
                "excluded_reason": item.get("excluded_reason"),
            }
        )
    return records


def _frontier_metadata(data: Any, threshold: float) -> dict[str, Any]:
    nelec_alpha = _safe_int(_data_get(data, "nelec_A"))
    nelec_beta = _safe_int(_data_get(data, "nelec_B"))
    if nelec_alpha is None or nelec_beta is None:
        return {
            "available": False,
            "reason": "electron counts are not available",
            "ambiguous": False,
        }

    current_open_pair = []
    if nelec_alpha > nelec_beta:
        current_open_pair = list(range(nelec_beta + 1, nelec_alpha + 1))

    gaps = {}
    min_abs_gap = None
    for spin, key in (("alpha", "OQP::E_MO_A"), ("beta", "OQP::E_MO_B")):
        energies = _as_float_list(_data_get(data, key))
        if not energies:
            continue
        spin_gaps = _spin_frontier_gaps(energies, nelec_alpha, nelec_beta)
        gaps[spin] = spin_gaps
        for value in spin_gaps.values():
            abs_value = abs(value)
            min_abs_gap = abs_value if min_abs_gap is None else min(min_abs_gap, abs_value)

    return {
        "available": True,
        "nelec_alpha": nelec_alpha,
        "nelec_beta": nelec_beta,
        "current_open_pair": current_open_pair,
        "gaps_hartree": gaps,
        "min_abs_gap_hartree": min_abs_gap,
        "min_abs_gap_ev": None if min_abs_gap is None else min_abs_gap * EV_PER_HARTREE,
        "ambiguous": bool(min_abs_gap is not None and min_abs_gap <= threshold),
    }


def _infer_nmo(data: Any, open_pairs: list[tuple[int, int]]) -> int:
    sizes = []
    for key in ("OQP::E_MO_A", "OQP::E_MO_B"):
        values = _as_float_list(_data_get(data, key))
        if values:
            sizes.append(len(values))
    nbf = _safe_int(_data_get(data, "nbf"))
    if nbf is not None:
        sizes.append(nbf)
    if sizes:
        return max(sizes)
    pair_indices = [idx for pair in open_pairs for idx in pair]
    return max(pair_indices) if pair_indices else 0


def _build_rohf_reference(
    ref_id: int,
    open_pair: tuple[int, int],
    weight: float,
    nelec_alpha: int,
    nelec_beta: int,
    nmo: int,
) -> dict[str, Any]:
    if max(open_pair) > nmo:
        return {
            "id": ref_id,
            "valid": False,
            "open_pair": [int(open_pair[0]), int(open_pair[1])],
            "reason": "open-shell pair exceeds available MO count",
            "nmo": nmo,
        }

    closed = []
    open_set = set(open_pair)
    for mo in range(1, nmo + 1):
        if mo in open_set:
            continue
        closed.append(mo)
        if len(closed) == nelec_beta:
            break

    if len(closed) != nelec_beta:
        return {
            "id": ref_id,
            "valid": False,
            "open_pair": [int(open_pair[0]), int(open_pair[1])],
            "reason": "not enough non-open orbitals to build the closed ROHF shell",
            "nmo": nmo,
            "nelec_beta": nelec_beta,
        }

    alpha_occupied = sorted(closed + list(open_pair))
    beta_occupied = list(closed)
    if len(alpha_occupied) != nelec_alpha:
        return {
            "id": ref_id,
            "valid": False,
            "open_pair": [int(open_pair[0]), int(open_pair[1])],
            "reason": "candidate pair does not produce the expected alpha-electron count",
            "alpha_occupied": alpha_occupied,
            "nelec_alpha": nelec_alpha,
        }

    return {
        "id": ref_id,
        "valid": True,
        "weight": float(weight),
        "open_pair": [int(open_pair[0]), int(open_pair[1])],
        "closed_orbitals": [int(item) for item in closed],
        "alpha_occupied": [int(item) for item in alpha_occupied],
        "beta_occupied": [int(item) for item in beta_occupied],
        "response_space": _mrsf_response_space(alpha_occupied, beta_occupied, open_pair, nmo),
    }


def _mrsf_response_space(
    alpha_occupied: list[int],
    beta_occupied: list[int],
    open_pair: tuple[int, int],
    nmo: int,
) -> dict[str, Any]:
    beta_occ = set(beta_occupied)
    open_set = set(open_pair)
    beta_virtual = [mo for mo in range(1, nmo + 1) if mo not in beta_occ]
    role_counts = {
        "closed_to_open": 0,
        "closed_to_virtual": 0,
        "open_to_open": 0,
        "open_to_virtual": 0,
    }
    for hole in alpha_occupied:
        hole_role = "open" if hole in open_set else "closed"
        for particle in beta_virtual:
            particle_role = "open" if particle in open_set else "virtual"
            role_counts[f"{hole_role}_to_{particle_role}"] += 1

    raw_dimension = len(alpha_occupied) * len(beta_virtual)
    return {
        "alpha_holes": [int(item) for item in alpha_occupied],
        "beta_particles": [int(item) for item in beta_virtual],
        "raw_dimension": raw_dimension,
        "singlet_dimension": max(raw_dimension - 1, 0),
        "triplet_dimension": max(raw_dimension - 3, 0),
        "role_counts": role_counts,
        "spin_adaptation": "one open-open pair coordinate retained; three redundant triplet open-open coordinates removed",
    }


def _nonzero_occupations(occupations: list[float]) -> list[dict[str, float | int]]:
    result = []
    for idx, occupation in enumerate(occupations, start=1):
        if abs(occupation) > 1.0e-14:
            result.append({"mo": int(idx), "occupation": float(occupation)})
    return result


def _dense_occupation_vector(items: Any, nmo: int, spin: str) -> list[float]:
    vector = [0.0 for _ in range(nmo)]
    for item in items:
        try:
            mo = int(item["mo"])
            occupation = float(item["occupation"])
        except (KeyError, TypeError, ValueError) as exc:
            raise MrsfReferenceError(f"invalid {spin} MRSF reference occupation entry") from exc
        if mo < 1 or mo > nmo:
            raise MrsfReferenceError(f"{spin} MRSF reference occupation MO index is out of range")
        if occupation < -1.0e-12 or occupation > 1.0 + 1.0e-12:
            raise MrsfReferenceError(f"{spin} MRSF reference occupation must be between 0 and 1")
        vector[mo - 1] = occupation
    return vector


def _spin_frontier_gaps(energies: list[float], nelec_alpha: int, nelec_beta: int) -> dict[str, float]:
    boundaries = {
        "closed_to_open": (nelec_beta, nelec_beta + 1),
        "open_to_virtual": (nelec_alpha, nelec_alpha + 1),
    }
    if nelec_alpha - nelec_beta >= 2:
        boundaries["open_pair"] = (nelec_beta + 1, nelec_beta + 2)

    result = {}
    nmo = len(energies)
    for name, (left, right) in boundaries.items():
        if 1 <= left <= nmo and 1 <= right <= nmo:
            result[name] = float(energies[right - 1] - energies[left - 1])
    return result


def _data_get(data: Any, key: str) -> Any:
    if data is None:
        return None
    try:
        return data[key]
    except Exception:
        return getattr(data, key, None)


def _safe_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _as_float_list(value: Any) -> list[float]:
    if value is None:
        return []
    try:
        return [float(item) for item in value]
    except TypeError:
        return []
