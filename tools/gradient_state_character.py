#!/usr/bin/env python3
"""Parse lightweight OpenQP TDDFT/MRSF state-character signatures.

This diagnostic helper is intentionally log-only: it helps decide whether a
finite-difference mismatch is likely root/state-character instability before
changing analytic-gradient assembly code.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any


_MRSF_STATE_RE = re.compile(
    r"State\s+#\s+(?P<root>\d+)\s+Energy\s+=\s+(?P<energy>[-+0-9.]+)\s+eV"
)
_S2_RE = re.compile(r"<S\^2>\s*=\s*(?P<s2>[-+0-9.]+)")
_CONFIG_RE = re.compile(
    r"^\s*(?P<drf>\d+)\s+(?P<coeff>[-+0-9.]+)\s+"
    r"(?P<occ>\d+)\s*->\s*(?P<vir>\d+)\s*$"
)
_TDDFT_ROW_RE = re.compile(
    r"^\s*0\s*->\s*(?P<root>\d+)\s+"
    r"(?P<energy>[-+0-9.]+)\s+"
    r"(?P<dx>[-+0-9.]+)\s+"
    r"(?P<dy>[-+0-9.]+)\s+"
    r"(?P<dz>[-+0-9.]+)\s+"
    r"(?P<osc>[-+0-9.]+)\s*$"
)


def parse_mrsf_states(text: str) -> dict[int, dict[str, Any]]:
    """Return MRSF states parsed from the spin-adapted spin-flip block."""
    states: dict[int, dict[str, Any]] = {}
    current: dict[str, Any] | None = None
    in_block = False

    for line in text.splitlines():
        if "Spin-adapted spin-flip excitations" in line:
            in_block = True
            continue
        if in_block and "Summary table" in line:
            break
        if not in_block:
            continue

        state_match = _MRSF_STATE_RE.search(line)
        if state_match:
            root = int(state_match.group("root"))
            current = {
                "root": root,
                "energy_ev": float(state_match.group("energy")),
                "s2": None,
                "configurations": [],
            }
            states[root] = current
            continue

        if current is None:
            continue

        s2_match = _S2_RE.search(line)
        if s2_match:
            current["s2"] = float(s2_match.group("s2"))
            continue

        config_match = _CONFIG_RE.match(line)
        if config_match:
            coeff = float(config_match.group("coeff"))
            current["configurations"].append(
                {
                    "drf": int(config_match.group("drf")),
                    "coeff": coeff,
                    "occ": int(config_match.group("occ")),
                    "vir": int(config_match.group("vir")),
                    "weight": coeff * coeff,
                }
            )

    for state in states.values():
        configs = sorted(
            state["configurations"], key=lambda item: abs(item["coeff"]), reverse=True
        )
        state["configurations"] = configs
        state["top_configuration"] = configs[0] if configs else None
    return states


def parse_tddft_transitions(text: str) -> dict[int, dict[str, Any]]:
    """Return conventional TDDFT transition signatures from the summary table."""
    transitions: dict[int, dict[str, Any]] = {}
    for line in text.splitlines():
        match = _TDDFT_ROW_RE.match(line)
        if not match:
            continue
        dipoles = {
            "x": float(match.group("dx")),
            "y": float(match.group("dy")),
            "z": float(match.group("dz")),
        }
        axis, value = max(dipoles.items(), key=lambda item: abs(item[1]))
        root = int(match.group("root"))
        transitions[root] = {
            "root": root,
            "energy_ev": float(match.group("energy")),
            "dipole": dipoles,
            "dominant_dipole_axis": axis,
            "dominant_dipole_abs": abs(value),
            "oscillator_strength": float(match.group("osc")),
        }
    return transitions


def _mrsf_signature(path: Path, root: int) -> dict[str, Any]:
    states = parse_mrsf_states(path.read_text(errors="replace"))
    if root not in states:
        raise SystemExit(f"{path}: MRSF root {root} was not found")
    state = states[root]
    top = state["top_configuration"]
    top_occ_vir = [top["occ"], top["vir"]] if top else None
    top_weight = top["weight"] if top else None
    return {
        "path": str(path),
        "root": root,
        "energy_ev": state["energy_ev"],
        "s2": state["s2"],
        "top_occ_vir": top_occ_vir,
        "top_drf": top["drf"] if top else None,
        "top_coeff": top["coeff"] if top else None,
        "top_weight": top_weight,
    }


def _tddft_signature(path: Path, root: int) -> dict[str, Any]:
    transitions = parse_tddft_transitions(path.read_text(errors="replace"))
    if root not in transitions:
        raise SystemExit(f"{path}: TDDFT root {root} was not found")
    transition = transitions[root]
    return {
        "path": str(path),
        "root": root,
        "energy_ev": transition["energy_ev"],
        "dominant_dipole_axis": transition["dominant_dipole_axis"],
        "dominant_dipole_abs": transition["dominant_dipole_abs"],
        "oscillator_strength": transition["oscillator_strength"],
        "dipole": transition["dipole"],
    }


def summarize_mrsf(paths: list[Path], target_root: int) -> dict[str, Any]:
    signatures = [_mrsf_signature(path, target_root) for path in paths]
    first = signatures[0]
    s2_values = [sig["s2"] for sig in signatures if sig["s2"] is not None]
    weight_values = [sig["top_weight"] for sig in signatures if sig["top_weight"] is not None]
    top_configs_match = all(sig["top_occ_vir"] == first["top_occ_vir"] for sig in signatures)
    max_s2_delta = max(s2_values) - min(s2_values) if s2_values else math.inf
    max_top_weight_delta = (
        max(weight_values) - min(weight_values) if weight_values else math.inf
    )
    return {
        "method": "mrsf",
        "target_root": target_root,
        "stable": bool(top_configs_match and max_s2_delta <= 0.05),
        "top_configs_match": top_configs_match,
        "max_s2_delta": max_s2_delta,
        "max_top_weight_delta": max_top_weight_delta,
        "signatures": signatures,
    }


def summarize_tddft(paths: list[Path], target_root: int) -> dict[str, Any]:
    signatures = [_tddft_signature(path, target_root) for path in paths]
    first = signatures[0]
    axes_match = all(
        sig["dominant_dipole_axis"] == first["dominant_dipole_axis"]
        for sig in signatures
    )
    energies = [sig["energy_ev"] for sig in signatures]
    osc = [sig["oscillator_strength"] for sig in signatures]
    max_energy_delta_ev = max(energies) - min(energies)
    max_oscillator_delta = max(osc) - min(osc)
    return {
        "method": "tddft",
        "target_root": target_root,
        "stable": bool(axes_match and max_energy_delta_ev <= 0.25),
        "dominant_axes_match": axes_match,
        "max_energy_delta_ev": max_energy_delta_ev,
        "max_oscillator_delta": max_oscillator_delta,
        "signatures": signatures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="method", required=True)
    for name in ("mrsf", "tddft"):
        sub = subparsers.add_parser(name)
        sub.add_argument("--target-root", type=int, required=True)
        sub.add_argument("logs", nargs="+", type=Path)
    args = parser.parse_args()

    if args.method == "mrsf":
        result = summarize_mrsf(args.logs, args.target_root)
    else:
        result = summarize_tddft(args.logs, args.target_root)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
