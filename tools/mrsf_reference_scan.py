#!/usr/bin/env python3
"""Run small MRSF ensemble-reference SCF continuity scans.

The first target is intentionally modest: a triplet H2O OH-stretch sanity scan
that compares ordinary ROHF against ensemble-reference SCF with equal and
gap-softmax weights.  The state-averaged MRSF response is still intentionally
guarded, so a NotImplementedError after SCF is counted as the expected endpoint
for ensemble variants.
"""

from __future__ import annotations

import argparse
import ast
import csv
from dataclasses import asdict, dataclass
from datetime import datetime
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = ROOT / "tools" / "_mrsf_reference_scan_scratch"

H2O_TRIPLET_BASE = (
    (8, 0.000000000, 0.000000000, -0.041061554),
    (1, -0.533194329, 0.533194329, -0.614469223),
    (1, 0.533194329, -0.533194329, -0.614469223),
)


@dataclass(frozen=True)
class Variant:
    key: str
    method: str
    mrsf_ref_mode: str
    weights: str
    weight_temperature: float | None = None
    expect_response_guard: bool = False


@dataclass
class ScanResult:
    point: int
    scale: float
    variant: str
    status: str
    returncode: int
    input: str
    log: str
    scf_energy: float | None = None
    scf_iterations: int | None = None
    scf_converged: bool = False
    scf_escalated: bool = False
    pair_selection: str | None = None
    open_pairs: Any = None
    reference_weights: Any = None
    applied_pairs: Any = None
    applied_weights: Any = None
    weight_model: str | None = None
    weight_temperature_hartree: float | None = None
    min_frontier_gap_hartree: float | None = None


VARIANTS = {
    "rohf": Variant("rohf", method="hf", mrsf_ref_mode="off", weights="equal"),
    "equal": Variant(
        "equal",
        method="tdhf",
        mrsf_ref_mode="state_average",
        weights="equal",
        expect_response_guard=True,
    ),
    "gap_softmax": Variant(
        "gap_softmax",
        method="tdhf",
        mrsf_ref_mode="state_average",
        weights="gap_softmax",
        weight_temperature=0.05,
        expect_response_guard=True,
    ),
}


def h2o_triplet_geometry(scale: float) -> list[tuple[int, float, float, float]]:
    """Scale both OH vectors around the oxygen atom."""

    oxygen = H2O_TRIPLET_BASE[0]
    geometry = [oxygen]
    ox, oy, oz = oxygen[1:]
    for atomic_number, x, y, z in H2O_TRIPLET_BASE[1:]:
        geometry.append(
            (
                atomic_number,
                ox + scale * (x - ox),
                oy + scale * (y - oy),
                oz + scale * (z - oz),
            )
        )
    return geometry


def format_geometry(geometry: list[tuple[int, float, float, float]]) -> str:
    return "\n".join(
        f"{atomic_number:2d} {x:16.9f} {y:16.9f} {z:16.9f}"
        for atomic_number, x, y, z in geometry
    )


def render_input(geometry: list[tuple[int, float, float, float]], variant: Variant) -> str:
    sections = [
        "[input]",
        "system=",
        format_geometry(geometry),
        "charge=0",
        f"method={variant.method}",
        "basis=6-31g",
        "runtype=energy",
        "functional=",
        "d4=False",
        "",
        "[guess]",
        "type=huckel",
        "save_mol=False",
        "",
        "[scf]",
        "type=rohf",
        "maxit=50",
        "forced_attempt=1",
        "maxdiis=5",
        "multiplicity=3",
        "conv=1.0e-8",
        "save_molden=False",
        "",
    ]

    if variant.method == "tdhf":
        sections.extend(
            [
                "[dftgrid]",
                "rad_npts=96",
                "ang_npts=302",
                "pruned=",
                "",
                "[tdhf]",
                "type=mrsf",
                "maxit=2",
                "multiplicity=1",
                "conv=1.0e-8",
                "nstate=1",
                "zvconv=1.0e-8",
                "",
            ]
        )

    if variant.mrsf_ref_mode != "off":
        sections.extend(
            [
                "[mrsf_ref]",
                f"mode={variant.mrsf_ref_mode}",
                "open_pairs=auto",
                f"weights={variant.weights}",
                f"weight_temperature={variant.weight_temperature or 0.05}",
                "max_refs=2",
                "gap_threshold=0.01",
                "overlap_threshold=0.85",
                "strict=False",
                "",
            ]
        )

    return "\n".join(sections).rstrip() + "\n"


def parse_scan_values(raw: str) -> list[float]:
    values = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if not values:
        raise ValueError("at least one scan point is required")
    return values


def parse_list_value(text: str) -> Any:
    try:
        return ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return text


def parse_optional_float(text: str) -> float | None:
    text = text.strip()
    if not text or text == "not available":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_log(log_path: Path) -> dict[str, Any]:
    parsed: dict[str, Any] = {}
    if not log_path.exists():
        return parsed

    final_energy_matches: list[tuple[float, int]] = []
    for line in log_path.read_text(errors="replace").splitlines():
        if "SCF convergence achieved" in line or "PyOQP: SCF converged" in line:
            parsed["scf_converged"] = True
        if "SCF not converged; escalating" in line:
            parsed["scf_escalated"] = True

        match = re.search(
            r"Final ROHF energy is\s+([-+0-9.Ee]+)\s+after\s+(\d+)\s+iterations",
            line,
        )
        if match:
            final_energy_matches.append((float(match.group(1)), int(match.group(2))))

        if "PyOQP MRSF pair selection:" in line:
            parsed["pair_selection"] = line.split(":", 1)[1].strip()
        elif "PyOQP MRSF reference open pairs:" in line:
            parsed["open_pairs"] = parse_list_value(line.split(":", 1)[1].strip())
        elif "PyOQP MRSF reference weights:" in line:
            parsed["reference_weights"] = parse_list_value(line.split(":", 1)[1].strip())
        elif "PyOQP MRSF SCF applied pairs:" in line:
            parsed["applied_pairs"] = parse_list_value(line.split(":", 1)[1].strip())
        elif "PyOQP MRSF SCF applied weights:" in line:
            parsed["applied_weights"] = parse_list_value(line.split(":", 1)[1].strip())
        elif "PyOQP MRSF weight model:" in line:
            parsed["weight_model"] = line.split(":", 1)[1].strip()
        elif "PyOQP MRSF weight temperature (Eh):" in line:
            parsed["weight_temperature_hartree"] = parse_optional_float(
                line.split(":", 1)[1]
            )
        elif "PyOQP MRSF min frontier gap (Eh):" in line:
            parsed["min_frontier_gap_hartree"] = parse_optional_float(
                line.split(":", 1)[1]
            )

    if final_energy_matches:
        parsed["scf_energy"], parsed["scf_iterations"] = final_energy_matches[-1]
    return parsed


def classify_status(proc: subprocess.CompletedProcess[str], variant: Variant) -> str:
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if proc.returncode == 0:
        return "ok"
    if variant.expect_response_guard and "mrsf_ref.mode=state_average" in combined:
        return "expected_response_guard"
    return "failed"


def run_case(
    openqp_cmd: list[str],
    point: int,
    scale: float,
    variant: Variant,
    case_dir: Path,
    dry_run: bool,
) -> ScanResult:
    case_dir.mkdir(parents=True, exist_ok=True)
    input_path = case_dir / f"h2o_triplet_s{point:03d}_{variant.key}.inp"
    log_path = input_path.with_suffix(".log")
    input_path.write_text(render_input(h2o_triplet_geometry(scale), variant))

    if dry_run:
        return ScanResult(
            point=point,
            scale=scale,
            variant=variant.key,
            status="dry_run",
            returncode=0,
            input=str(input_path),
            log=str(log_path),
        )

    proc = subprocess.run(
        [*openqp_cmd, input_path.name],
        cwd=case_dir,
        text=True,
        capture_output=True,
        check=False,
    )
    (case_dir / f"{input_path.stem}.stdout").write_text(proc.stdout or "")
    (case_dir / f"{input_path.stem}.stderr").write_text(proc.stderr or "")

    fields = parse_log(log_path)
    return ScanResult(
        point=point,
        scale=scale,
        variant=variant.key,
        status=classify_status(proc, variant),
        returncode=proc.returncode,
        input=str(input_path),
        log=str(log_path),
        **fields,
    )


def write_summary(results: list[ScanResult], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = [asdict(result) for result in results]
    (outdir / "summary.json").write_text(json.dumps(rows, indent=2, sort_keys=True))

    fields = list(rows[0].keys()) if rows else [field.name for field in ScanResult.__dataclass_fields__.values()]
    with (outdir / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_table(results: list[ScanResult]) -> None:
    print("point,scale,variant,status,energy,iter,converged,escalated,pairs,applied_weights")
    for item in results:
        print(
            f"{item.point},{item.scale:.6f},{item.variant},{item.status},"
            f"{item.scf_energy},{item.scf_iterations},{item.scf_converged},"
            f"{item.scf_escalated},{item.open_pairs},{item.applied_weights}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a small MRSF ensemble-reference SCF scan."
    )
    parser.add_argument(
        "--points",
        default="0.98,1.00,1.02",
        help="comma-separated OH stretch scale factors for the H2O triplet sanity scan",
    )
    parser.add_argument(
        "--variants",
        default="rohf,equal,gap_softmax",
        help=f"comma-separated variants from: {', '.join(sorted(VARIANTS))}",
    )
    parser.add_argument(
        "--outdir",
        default="",
        help="output directory; default is a timestamped scratch directory",
    )
    parser.add_argument(
        "--openqp",
        default="openqp --nompi --omp 1",
        help="OpenQP command used for each run",
    )
    parser.add_argument("--dry-run", action="store_true", help="write inputs without running OpenQP")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    points = parse_scan_values(args.points)
    variant_keys = [item.strip() for item in args.variants.split(",") if item.strip()]
    unknown = [key for key in variant_keys if key not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir) if args.outdir else DEFAULT_SCRATCH / f"h2o_triplet_{timestamp}"
    outdir.mkdir(parents=True, exist_ok=True)

    openqp_cmd = shlex.split(args.openqp)
    results: list[ScanResult] = []
    for point, scale in enumerate(points):
        for variant_key in variant_keys:
            result = run_case(
                openqp_cmd=openqp_cmd,
                point=point,
                scale=scale,
                variant=VARIANTS[variant_key],
                case_dir=outdir / f"point_{point:03d}",
                dry_run=args.dry_run,
            )
            results.append(result)
            print(
                f"[{result.status}] point={point} scale={scale:.6f} "
                f"variant={variant_key} energy={result.scf_energy}"
            )
            if result.status == "failed":
                print(f"  input: {result.input}", file=sys.stderr)
                print(f"  log:   {result.log}", file=sys.stderr)

    write_summary(results, outdir)
    print_table(results)
    print(f"summary: {outdir / 'summary.csv'}")
    return 1 if any(result.status == "failed" for result in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
