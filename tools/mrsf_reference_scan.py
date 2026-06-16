#!/usr/bin/env python3
"""Run small MRSF ensemble-reference continuity scans.

The first target was intentionally modest: a triplet H2O OH-stretch sanity scan
that compares ordinary ROHF against ensemble-reference SCF with equal and
gap-softmax weights.  The ethylene torsion target is the next, more relevant
near-degeneracy probe.  Ensemble variants now use the energy-only
state-interaction MRSF response prototype.
"""

from __future__ import annotations

import argparse
import ast
import csv
from dataclasses import asdict, dataclass
from datetime import datetime
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = ROOT / "tools" / "_mrsf_reference_scan_scratch"

H2O_TRIPLET_BASE = (
    (8, 0.000000000, 0.000000000, -0.041061554),
    (1, -0.533194329, 0.533194329, -0.614469223),
    (1, 0.533194329, -0.533194329, -0.614469223),
)

ETHYLENE_BASE = (
    (6, 0.000000000, 0.000000000, 0.665000000),
    (6, 0.000000000, 0.000000000, -0.665000000),
    (1, 0.000000000, 0.920000000, 1.230000000),
    (1, 0.000000000, -0.920000000, 1.230000000),
    (1, 0.000000000, 0.920000000, -1.230000000),
    (1, 0.000000000, -0.920000000, -1.230000000),
)


@dataclass(frozen=True)
class Variant:
    key: str
    method: str
    mrsf_ref_mode: str
    weights: str
    weight_temperature: float | None = None
    coupling: str = "overlap_offdiagonal"
    expect_response_guard: bool = False


@dataclass(frozen=True)
class ScanTarget:
    key: str
    label: str
    coordinate_name: str
    default_points: str
    filename_stem: str
    geometry: Callable[[float], list[tuple[int, float, float, float]]]


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
    response_energy: float | None = None
    state_energy: float | None = None
    scf_iterations: int | None = None
    scf_converged: bool = False
    scf_escalated: bool = False
    mrsf_converged_blocks: int = 0
    mrsf_block_iterations: Any = None
    pair_selection: str | None = None
    open_pairs: Any = None
    reference_weights: Any = None
    applied_pairs: Any = None
    applied_weights: Any = None
    weight_model: str | None = None
    weight_temperature_hartree: float | None = None
    min_frontier_gap_hartree: float | None = None
    response_status: str | None = None
    response_model: str | None = None
    response_coupling: str | None = None
    response_coupled: bool | None = None
    full_response_kernel: bool | None = None
    response_energy_only: bool | None = None
    response_offdiagonal_count: int | None = None
    response_max_abs_offdiagonal_hamiltonian: float | None = None
    response_candidate_count: int | None = None
    response_raw_candidate_count: int | None = None
    response_skipped_blocks: Any = None
    response_si_common_dimension: int | None = None
    response_si_kept_states: int | None = None
    response_si_dropped_redundant: int | None = None
    response_si_dropped_floor: int | None = None
    response_si_metric_min_eig: float | None = None
    response_selected_states: Any = None
    dominant_open_pair: Any = None
    scan: str | None = None
    coordinate_label: str | None = None


VARIANTS = {
    "rohf": Variant("rohf", method="hf", mrsf_ref_mode="off", weights="equal"),
    "mrsf": Variant("mrsf", method="tdhf", mrsf_ref_mode="off", weights="equal"),
    "equal": Variant(
        "equal",
        method="tdhf",
        mrsf_ref_mode="ensemble",
        weights="equal",
    ),
    "equal_block": Variant(
        "equal_block",
        method="tdhf",
        mrsf_ref_mode="ensemble",
        weights="equal",
        coupling="block_diagonal",
    ),
    "gap_softmax": Variant(
        "gap_softmax",
        method="tdhf",
        mrsf_ref_mode="ensemble",
        weights="gap_softmax",
        weight_temperature=0.05,
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


def _rotate_xy(x: float, y: float, angle_degrees: float) -> tuple[float, float]:
    radians = math.radians(angle_degrees)
    cosine = math.cos(radians)
    sine = math.sin(radians)
    return x * cosine - y * sine, x * sine + y * cosine


def ethylene_torsion_geometry(angle_degrees: float) -> list[tuple[int, float, float, float]]:
    """Twist the two CH2 planes symmetrically around the C=C axis.

    The carbon atoms sit on the z axis.  The hydrogens on the +z carbon rotate
    by +angle/2 and the hydrogens on the -z carbon rotate by -angle/2, so the
    relative H-C-C-H torsion is the requested angle.
    """

    geometry: list[tuple[int, float, float, float]] = []
    half_angle = 0.5 * angle_degrees
    for atom_index, (atomic_number, x, y, z) in enumerate(ETHYLENE_BASE):
        if atom_index in {2, 3}:
            x, y = _rotate_xy(x, y, half_angle)
        elif atom_index in {4, 5}:
            x, y = _rotate_xy(x, y, -half_angle)
        geometry.append((atomic_number, x, y, z))
    return geometry


def o2_dissociation_geometry(distance_angstrom: float) -> list[tuple[int, float, float, float]]:
    """Place triplet O2 along the z axis at the requested bond distance."""

    half_distance = 0.5 * distance_angstrom
    return [
        (8, 0.000000000, 0.000000000, -half_distance),
        (8, 0.000000000, 0.000000000, half_distance),
    ]


def allene_torsion_geometry(dihedral_degrees: float) -> list[tuple[int, float, float, float]]:
    """Allene (H2C=C=CH2) with the two terminal CH2 groups at a chosen dihedral.

    Central C at the origin, terminal carbons on the +/-z axis.  The +z CH2 is
    fixed in the xz plane; the -z CH2 is rotated by ``dihedral_degrees`` about the
    C=C=C axis.  90 deg is the D2d ground state (perpendicular CH2 groups); 0 deg
    is the planar D2h rotation barrier -- a four-orbital pi diradical and the
    CAS(4,4) multireference point.
    """

    cc = 1.308  # C=C bond length (A)
    ch = 1.087  # C-H bond length (A)
    out = math.radians(180.0 - 121.0)  # C-H angle from the outward C=C axis
    cz = ch * math.cos(out)  # C-H component along the outward z axis
    cp = ch * math.sin(out)  # C-H component in the CH2 plane
    phi = math.radians(dihedral_degrees)
    cphi, sphi = math.cos(phi), math.sin(phi)
    return [
        (6, 0.0, 0.0, 0.0),            # central C
        (6, 0.0, 0.0, cc),            # +z terminal C
        (6, 0.0, 0.0, -cc),           # -z terminal C
        (1, cp, 0.0, cc + cz),        # +z CH2 fixed in the xz plane
        (1, -cp, 0.0, cc + cz),
        (1, cp * cphi, cp * sphi, -cc - cz),   # -z CH2 rotated by the dihedral
        (1, -cp * cphi, -cp * sphi, -cc - cz),
    ]


SCAN_TARGETS = {
    "h2o_triplet": ScanTarget(
        key="h2o_triplet",
        label="triplet H2O OH-stretch sanity scan",
        coordinate_name="oh_scale",
        default_points="0.98,1.00,1.02",
        filename_stem="h2o_triplet",
        geometry=h2o_triplet_geometry,
    ),
    "ethylene_torsion": ScanTarget(
        key="ethylene_torsion",
        label="triplet ethylene torsion near-degeneracy scan",
        coordinate_name="torsion_degrees",
        default_points="85,90,95",
        filename_stem="ethylene_torsion",
        geometry=ethylene_torsion_geometry,
    ),
    "o2_dissociation": ScanTarget(
        key="o2_dissociation",
        label="triplet O2 dissociation scan",
        coordinate_name="bond_distance_angstrom",
        default_points="1.10,1.21,1.40,1.70,2.10,2.60,3.20",
        filename_stem="o2_dissociation",
        geometry=o2_dissociation_geometry,
    ),
    "allene_torsion": ScanTarget(
        key="allene_torsion",
        label="triplet allene torsion CAS(4,4) scan",
        coordinate_name="dihedral_degrees",
        default_points="0,30,45,60,90",
        filename_stem="allene_torsion",
        geometry=allene_torsion_geometry,
    ),
}


def format_geometry(geometry: list[tuple[int, float, float, float]]) -> str:
    return "\n".join(
        f"{atomic_number:2d} {x:16.9f} {y:16.9f} {z:16.9f}"
        for atomic_number, x, y, z in geometry
    )


def render_input(
    geometry: list[tuple[int, float, float, float]],
    variant: Variant,
    open_pairs: str = "auto",
    max_refs: int = 6,
) -> str:
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
                "maxit=60",
                "multiplicity=1",
                "conv=1.0e-6",
                "nstate=1",
                "zvconv=1.0e-6",
                "",
            ]
        )

    if variant.mrsf_ref_mode != "off":
        sections.extend(
            [
                "[mrsf_ref]",
                f"mode={variant.mrsf_ref_mode}",
                f"open_pairs={open_pairs}",
                f"weights={variant.weights}",
                f"weight_temperature={variant.weight_temperature or 0.05}",
                f"max_refs={max_refs}",
                "gap_threshold=0.01",
                "overlap_threshold=0.85",
                "trial_vectors=adaptive",
                "trial_shift=1.0e6",
                f"coupling={variant.coupling}",
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


def parse_optional_int(text: str) -> int | None:
    text = text.strip()
    if not text or text == "not available":
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_yes_no(text: str) -> bool | None:
    value = text.strip().lower()
    if value in {"yes", "true", "1"}:
        return True
    if value in {"no", "false", "0"}:
        return False
    return None


def first_selected_response_energy(selected_states: Any) -> float | None:
    if not isinstance(selected_states, list) or not selected_states:
        return None
    first = selected_states[0]
    if not isinstance(first, dict):
        return None
    try:
        return float(first["energy"])
    except (KeyError, TypeError, ValueError):
        return None


def first_dominant_open_pair(selected_states: Any) -> Any:
    if not isinstance(selected_states, list) or not selected_states:
        return None
    first = selected_states[0]
    if isinstance(first, dict):
        return first.get("dominant_open_pair")
    return None


def parse_log(log_path: Path) -> dict[str, Any]:
    parsed: dict[str, Any] = {}
    if not log_path.exists():
        return parsed

    final_energy_matches: list[tuple[float, int]] = []
    state_energies: dict[int, float] = {}
    mrsf_iterations: list[int] = []
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

        match = re.search(r"MRSF-TD-DFT energies converged in\s+(\d+)\s+iterations", line)
        if match:
            mrsf_iterations.append(int(match.group(1)))

        match = re.search(r"PyOQP state\s+(\d+)\s+([-+0-9.Ee]+)", line)
        if match:
            state_energies[int(match.group(1))] = float(match.group(2))

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
        elif "PyOQP MRSF response status:" in line:
            parsed["response_status"] = line.split(":", 1)[1].strip()
        elif "PyOQP MRSF response model:" in line:
            parsed["response_model"] = line.split(":", 1)[1].strip()
        elif "PyOQP MRSF response coupling:" in line:
            parsed["response_coupling"] = line.split(":", 1)[1].strip()
        elif "PyOQP MRSF response coupled:" in line:
            parsed["response_coupled"] = parse_yes_no(line.split(":", 1)[1])
        elif "PyOQP MRSF full response kernel:" in line:
            parsed["full_response_kernel"] = parse_yes_no(line.split(":", 1)[1])
        elif "PyOQP MRSF response energy only:" in line:
            parsed["response_energy_only"] = parse_yes_no(line.split(":", 1)[1])
        elif "PyOQP MRSF offdiag couplings:" in line:
            parsed["response_offdiagonal_count"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF max abs offdiag H (Eh):" in line:
            parsed["response_max_abs_offdiagonal_hamiltonian"] = parse_optional_float(
                line.split(":", 1)[1]
            )
        elif "PyOQP MRSF selected states:" in line:
            selected_states = parse_list_value(line.split(":", 1)[1].strip())
            parsed["response_selected_states"] = selected_states
            parsed["response_energy"] = first_selected_response_energy(selected_states)
            parsed["dominant_open_pair"] = first_dominant_open_pair(selected_states)
        elif "PyOQP MRSF candidate states:" in line:
            parsed["response_candidate_count"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF raw candidate states:" in line:
            parsed["response_raw_candidate_count"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF skipped blocks:" in line:
            parsed["response_skipped_blocks"] = parse_list_value(line.split(":", 1)[1].strip())
        elif "PyOQP MRSF SI common dimension:" in line:
            parsed["response_si_common_dimension"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF SI kept states:" in line:
            parsed["response_si_kept_states"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF SI dropped redundant:" in line:
            parsed["response_si_dropped_redundant"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF SI dropped below floor:" in line:
            parsed["response_si_dropped_floor"] = parse_optional_int(line.split(":", 1)[1])
        elif "PyOQP MRSF SI metric min eig:" in line:
            parsed["response_si_metric_min_eig"] = parse_optional_float(line.split(":", 1)[1])

    if final_energy_matches:
        parsed["scf_energy"], parsed["scf_iterations"] = final_energy_matches[-1]
    if mrsf_iterations:
        parsed["mrsf_block_iterations"] = mrsf_iterations
        parsed["mrsf_converged_blocks"] = len(mrsf_iterations)
    if state_energies:
        parsed["state_energy"] = state_energies.get(1, state_energies.get(max(state_energies)))
        if "response_energy" not in parsed and 0 in state_energies and parsed["state_energy"] is not None:
            parsed["response_energy"] = parsed["state_energy"] - state_energies[0]
    return parsed


def classify_status(proc: subprocess.CompletedProcess[str], variant: Variant) -> str:
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if proc.returncode == 0:
        return "ok"
    if variant.expect_response_guard and (
        "mrsf_ref.mode=ensemble" in combined or "mrsf_ref.mode=state_average" in combined
    ):
        return "expected_response_guard"
    return "failed"


def run_case(
    openqp_cmd: list[str],
    target: ScanTarget,
    point: int,
    coordinate: float,
    variant: Variant,
    case_dir: Path,
    dry_run: bool,
    open_pairs: str,
    max_refs: int,
) -> ScanResult:
    case_dir.mkdir(parents=True, exist_ok=True)
    input_path = case_dir / f"{target.filename_stem}_s{point:03d}_{variant.key}.inp"
    log_path = input_path.with_suffix(".log")
    input_path.write_text(
        render_input(
            target.geometry(coordinate),
            variant,
            open_pairs=open_pairs,
            max_refs=max_refs,
        )
    )

    if dry_run:
        return ScanResult(
            point=point,
            scale=coordinate,
            variant=variant.key,
            status="dry_run",
            returncode=0,
            input=str(input_path),
            log=str(log_path),
            scan=target.key,
            coordinate_label=target.coordinate_name,
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
        scale=coordinate,
        variant=variant.key,
        status=classify_status(proc, variant),
        returncode=proc.returncode,
        input=str(input_path),
        log=str(log_path),
        scan=target.key,
        coordinate_label=target.coordinate_name,
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
    coordinate_label = results[0].coordinate_label if results else "coordinate"
    print(
        f"point,{coordinate_label},variant,status,scf_energy,response_energy,"
        "state_energy,iter,converged,escalated,pairs,applied_weights,response_model"
    )
    for item in results:
        print(
            f"{item.point},{item.scale:.6f},{item.variant},{item.status},"
            f"{item.scf_energy},{item.response_energy},{item.state_energy},"
            f"{item.scf_iterations},{item.scf_converged},{item.scf_escalated},"
            f"{item.open_pairs},{item.applied_weights},{item.response_model}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a small MRSF ensemble-reference scan."
    )
    parser.add_argument(
        "--scan",
        default="h2o_triplet",
        choices=sorted(SCAN_TARGETS),
        help="scan target to run",
    )
    parser.add_argument(
        "--points",
        default="",
        help="comma-separated scan coordinates; defaults depend on --scan",
    )
    parser.add_argument(
        "--variants",
        default="rohf,mrsf,equal,gap_softmax",
        help=f"comma-separated variants from: {', '.join(sorted(VARIANTS))}",
    )
    parser.add_argument(
        "--open-pairs",
        default="auto",
        help="open-shell MO pairs for ensemble variants, e.g. '8:9;7:10'; default auto",
    )
    parser.add_argument(
        "--max-refs",
        type=int,
        default=6,
        help="maximum automatic references for ensemble variants",
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
    target = SCAN_TARGETS[args.scan]
    points = parse_scan_values(args.points or target.default_points)
    variant_keys = [item.strip() for item in args.variants.split(",") if item.strip()]
    unknown = [key for key in variant_keys if key not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variant(s): {', '.join(unknown)}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir) if args.outdir else DEFAULT_SCRATCH / f"{target.key}_{timestamp}"
    outdir.mkdir(parents=True, exist_ok=True)

    openqp_cmd = shlex.split(args.openqp)
    results: list[ScanResult] = []
    for point, coordinate in enumerate(points):
        for variant_key in variant_keys:
            result = run_case(
                openqp_cmd=openqp_cmd,
                target=target,
                point=point,
                coordinate=coordinate,
                variant=VARIANTS[variant_key],
                case_dir=outdir / f"point_{point:03d}",
                dry_run=args.dry_run,
                open_pairs=args.open_pairs,
                max_refs=args.max_refs,
            )
            results.append(result)
            print(
                f"[{result.status}] {target.coordinate_name}={coordinate:.6f} "
                f"variant={variant_key} state={result.state_energy} scf={result.scf_energy}"
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
