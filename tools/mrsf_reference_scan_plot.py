#!/usr/bin/env python3
"""Plot MRSF ensemble-reference scan summaries."""

from __future__ import annotations

import argparse
import ast
import csv
from pathlib import Path
from typing import Any


def parse_list(text: str) -> Any:
    if not text:
        return None
    try:
        return ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return text


def load_rows(summary: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with summary.open(newline="") as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, Any] = dict(row)
            parsed["coordinate"] = float(row.get("scale") or row.get("coordinate") or 0.0)
            parsed["scf_energy"] = _optional_float(row.get("scf_energy", ""))
            parsed["response_energy"] = _optional_float(row.get("response_energy", ""))
            parsed["state_energy"] = _optional_float(row.get("state_energy", ""))
            parsed["scf_iterations"] = _optional_int(row.get("scf_iterations", ""))
            parsed["mrsf_converged_blocks"] = _optional_int(row.get("mrsf_converged_blocks", ""))
            parsed["scf_escalated"] = str(row.get("scf_escalated", "")).lower() == "true"
            parsed["applied_weights"] = parse_list(row.get("applied_weights", ""))
            parsed["applied_pairs"] = parse_list(row.get("applied_pairs", ""))
            parsed["open_pairs"] = parse_list(row.get("open_pairs", ""))
            parsed["dominant_open_pair"] = parse_list(row.get("dominant_open_pair", ""))
            rows.append(parsed)
    return rows


def plot_summary(summary: Path, outdir: Path, prefix: str) -> list[Path]:
    rows = load_rows(summary)
    if not rows:
        raise ValueError(f"no rows found in {summary}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    coordinate_label = rows[0].get("coordinate_label") or "coordinate"

    energy_path = outdir / f"{prefix}_energy.pdf"
    energy_column = _preferred_energy_column(rows)
    _plot_relative_energy(rows, coordinate_label, energy_column, energy_path, plt)

    weights_path = outdir / f"{prefix}_weights_iterations.pdf"
    _plot_weights_and_iterations(rows, coordinate_label, weights_path, plt)

    return [energy_path, weights_path]


def _plot_relative_energy(
    rows: list[dict[str, Any]],
    coordinate_label: str,
    energy_column: str,
    path: Path,
    plt: Any,
) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    for variant in _ordered_variants(rows):
        series = _series(rows, variant)
        if not series:
            continue
        coordinates = [row["coordinate"] for row in series]
        energies = [_energy_for_row(row, energy_column) for row in series]
        if any(energy is None for energy in energies):
            continue
        center = min(range(len(coordinates)), key=lambda idx: abs(coordinates[idx] - _midpoint(coordinates)))
        reference = energies[center]
        relative_mhartree = [(energy - reference) * 1000.0 for energy in energies]
        ax.plot(coordinates, relative_mhartree, marker="o", linewidth=1.8, label=variant)

    ax.axhline(0.0, color="0.75", linewidth=0.8)
    ax.set_xlabel(_axis_label(coordinate_label))
    ax.set_ylabel(f"Relative {_energy_label(energy_column)} (mEh)")
    ax.legend(frameon=False)
    ax.grid(True, linewidth=0.4, alpha=0.35)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_weights_and_iterations(rows: list[dict[str, Any]], coordinate_label: str, path: Path, plt: Any) -> None:
    fig, (weight_ax, iter_ax) = plt.subplots(2, 1, figsize=(6.5, 5.2), sharex=True)

    for variant in _ordered_variants(rows):
        series = _series(rows, variant)
        coordinates = [row["coordinate"] for row in series]
        first_weights = []
        for row in series:
            weights = row.get("applied_weights")
            if isinstance(weights, list) and weights:
                first_weights.append(float(weights[0]))
            else:
                first_weights.append(None)
        if any(weight is not None for weight in first_weights):
            weight_ax.plot(
                coordinates,
                first_weights,
                marker="o",
                linewidth=1.8,
                label=f"{variant} ref 1",
            )

        iterations = [row.get("scf_iterations") for row in series]
        if any(iteration is not None for iteration in iterations):
            iter_ax.plot(coordinates, iterations, marker="s", linewidth=1.5, label=variant)
            for coordinate, iteration, row in zip(coordinates, iterations, series):
                if iteration is not None and row.get("scf_escalated"):
                    iter_ax.annotate(
                        "SOSCF",
                        (coordinate, iteration),
                        textcoords="offset points",
                        xytext=(0, 6),
                        ha="center",
                        fontsize=8,
                    )

    weight_ax.set_ylabel("Applied weight")
    weight_ax.set_ylim(-0.05, 1.05)
    weight_ax.legend(frameon=False, ncol=2)
    weight_ax.grid(True, linewidth=0.4, alpha=0.35)

    iter_ax.set_xlabel(_axis_label(coordinate_label))
    iter_ax.set_ylabel("SCF iterations")
    iter_ax.legend(frameon=False, ncol=3)
    iter_ax.grid(True, linewidth=0.4, alpha=0.35)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _ordered_variants(rows: list[dict[str, Any]]) -> list[str]:
    preferred = ["rohf", "mrsf", "equal", "gap_softmax"]
    variants = {str(row.get("variant", "")) for row in rows}
    ordered = [variant for variant in preferred if variant in variants]
    ordered.extend(sorted(variants.difference(ordered)))
    return ordered


def _series(rows: list[dict[str, Any]], variant: str) -> list[dict[str, Any]]:
    return sorted(
        [row for row in rows if row.get("variant") == variant],
        key=lambda row: row["coordinate"],
    )


def _midpoint(values: list[float]) -> float:
    return 0.5 * (min(values) + max(values))


def _preferred_energy_column(rows: list[dict[str, Any]]) -> str:
    if any(row.get("state_energy") is not None for row in rows):
        return "state_energy"
    if any(row.get("response_energy") is not None for row in rows):
        return "response_energy"
    return "scf_energy"


def _energy_for_row(row: dict[str, Any], energy_column: str) -> float | None:
    value = row.get(energy_column)
    if value is not None:
        return value
    if energy_column == "state_energy":
        return row.get("scf_energy")
    return None


def _energy_label(energy_column: str) -> str:
    if energy_column == "state_energy":
        return "state energy"
    if energy_column == "response_energy":
        return "response energy"
    return "SCF energy"


def _axis_label(label: str) -> str:
    if label == "torsion_degrees":
        return "Torsion angle (degrees)"
    if label == "oh_scale":
        return "O-H scale factor"
    return label.replace("_", " ")


def _optional_float(text: str | None) -> float | None:
    if text in {None, "", "None"}:
        return None
    return float(text)


def _optional_int(text: str | None) -> int | None:
    if text in {None, "", "None"}:
        return None
    return int(text)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot MRSF reference scan summary CSVs.")
    parser.add_argument("summary", type=Path, help="summary.csv from tools/mrsf_reference_scan.py")
    parser.add_argument("--outdir", type=Path, default=None, help="output directory; default is summary directory")
    parser.add_argument("--prefix", default=None, help="output filename prefix; default is scan name")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = load_rows(args.summary)
    if not rows:
        raise SystemExit(f"no rows found in {args.summary}")
    outdir = args.outdir or args.summary.parent
    prefix = args.prefix or str(rows[0].get("scan") or args.summary.parent.name)
    for path in plot_summary(args.summary, outdir, prefix):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
