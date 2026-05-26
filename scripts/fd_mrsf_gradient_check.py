#!/usr/bin/env python3
"""Finite-difference diagnostic for MRSF-TDDFT analytic gradients.

The script is intentionally standalone so it can be copied into benchmark
artifact directories.  It writes OpenQP inputs, runs central finite-difference
energy jobs, runs analytic gradient jobs, and summarizes component-wise errors.

State numbering follows OpenQP/MRSF response-root numbering:
  root 2 = physical S1, root 3 = physical S2.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable

BOHR_PER_ANG = 1.8897261254576558
VARIANTS = ("default", "spc_ovov_0", "flip_ovov_sign")

BUILTIN_GEOMETRIES = {
    "h2o": [
        ("O", 0.000000000, 0.000000000, -0.041061554),
        ("H", -0.533194329, 0.533194329, -0.614469223),
        ("H", 0.533194329, -0.533194329, -0.614469223),
    ],
    "nh3": [
        ("N", 0.000000000, 0.000000000, 0.116489000),
        ("H", 0.000000000, 0.939731000, -0.271808000),
        ("H", 0.813831000, -0.469866000, -0.271808000),
        ("H", -0.813831000, -0.469866000, -0.271808000),
    ],
}


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"empty XYZ: {path}")
    natom = int(lines[0].strip())
    atoms = []
    for line in lines[2 : 2 + natom]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    if len(atoms) != natom:
        raise ValueError(f"expected {natom} atoms in {path}, found {len(atoms)}")
    return atoms


def load_geometry(name_or_path: str) -> list[tuple[str, float, float, float]]:
    if name_or_path in BUILTIN_GEOMETRIES:
        return list(BUILTIN_GEOMETRIES[name_or_path])
    return read_xyz(Path(name_or_path))


def physical_state(root: int) -> str:
    return f"S{root - 1}"


def spc_lines_for_variant(variant: str) -> str:
    if variant == "default":
        return ""
    if variant == "spc_ovov_0":
        return "spc_ovov=0\n"
    if variant == "flip_ovov_sign":
        return ""
    raise ValueError(f"unknown variant: {variant}")


def render_input(
    *,
    geom: list[tuple[str, float, float, float]],
    runtype: str,
    target_root: int,
    variant: str,
    basis: str,
    functional: str,
    nstate: int,
    z_solver: int,
    huckel: bool = True,
) -> str:
    geom_text = "".join(f"   {sym:2s} {x: .9f} {y: .9f} {z: .9f}\n" for sym, x, y, z in geom)
    guess = "[guess]\ntype=huckel\n\n" if huckel else ""
    props = f"\n[properties]\ngrad={target_root}\n" if runtype == "grad" else ""
    return f"""[input]
system=
{geom_text}charge=0
runtype={runtype}
basis={basis}
functional={functional}
method=tdhf

{guess}[scf]
multiplicity=3
type=rohf
maxit=50
converger_type=diis
alternative_scf=trah
forced_attempt=2

[tdhf]
type=mrsf
nstate={nstate}
zvconv=1e-7
z_solver={z_solver}
{spc_lines_for_variant(variant)}{props}"""


def central_fd_ha_per_bohr(e_plus_ha: float, e_minus_ha: float, h_ang: float) -> float:
    return (e_plus_ha - e_minus_ha) / (2.0 * h_ang) / BOHR_PER_ANG


def displaced(
    geom: list[tuple[str, float, float, float]], component: int, delta_ang: float
) -> list[tuple[str, float, float, float]]:
    atom = component // 3
    axis = component % 3
    rows: list[list[str | float]] = [list(row) for row in geom]
    rows[atom][1 + axis] = float(rows[atom][1 + axis]) + delta_ang
    return [(str(row[0]), float(row[1]), float(row[2]), float(row[3])) for row in rows]


def openqp_command(command_template: str | None, input_name: str) -> list[str]:
    if command_template:
        return command_template.format(input=input_name).split()
    return ["openqp", "--nompi", input_name]


def run_openqp_job(
    *,
    job_dir: Path,
    name: str,
    geom: list[tuple[str, float, float, float]],
    runtype: str,
    target_root: int,
    variant: str,
    basis: str,
    functional: str,
    nstate: int,
    z_solver: int,
    threads: int,
    command_template: str | None,
    force: bool,
    timeout_s: int,
) -> dict:
    job_dir.mkdir(parents=True, exist_ok=True)
    parsed_json = job_dir / "parsed.json"
    if parsed_json.exists() and not force:
        return json.loads(parsed_json.read_text())

    inp = job_dir / f"{name}.inp"
    log = job_dir / f"{name}.log"
    inp.write_text(
        render_input(
            geom=geom,
            runtype=runtype,
            target_root=target_root,
            variant=variant,
            basis=basis,
            functional=functional,
            nstate=nstate,
            z_solver=z_solver,
        )
    )
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": str(threads),
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "OMP_PROC_BIND": "false",
        }
    )
    t0 = time.time()
    with log.open("w") as handle:
        proc = subprocess.run(
            openqp_command(command_template, inp.name),
            cwd=job_dir,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            timeout=timeout_s,
        )
    elapsed = time.time() - t0
    text = log.read_text(errors="replace")
    parsed = {
        "name": name,
        "variant": variant,
        "runtype": runtype,
        "root": target_root,
        "returncode": proc.returncode,
        "elapsed_s": elapsed,
        "log": str(log),
        "trah": bool(re.search(r"SCF did not converge|TRAH / Trust-Region", text)),
        "crash": bool(re.search(r"Segmentation fault|SIGSEGV|Bus error|SIGBUS|Fatal Python error", text, re.I)),
    }
    energies = {int(m.group(1)): float(m.group(2)) for m in re.finditer(r"PyOQP state\s+(\d+)\s+([-+0-9.]+)", text)}
    parsed["energies_ha"] = energies
    parsed["target_energy_ha"] = energies.get(target_root)

    if runtype == "grad":
        parsed["grad_ha_per_bohr"] = parse_last_gradient(text, target_root, len(geom))

    parsed["ok"] = parsed["returncode"] == 0 and parsed.get("target_energy_ha") is not None and (
        runtype != "grad" or parsed.get("grad_ha_per_bohr") is not None
    )
    if not parsed["ok"]:
        parsed["tail"] = "\n".join(text.splitlines()[-80:])
    parsed_json.write_text(json.dumps(parsed, indent=2))
    return parsed


def parse_last_gradient(text: str, target_root: int, natom: int) -> list[float] | None:
    marker = f"PyOQP state {target_root}"
    lines = text.splitlines()
    grad = None
    atom_line = re.compile(r"\s*[A-Z][a-z]?\s+[-+0-9.]+\s+[-+0-9.]+\s+[-+0-9.]+")
    for idx, line in enumerate(lines):
        if marker not in line or idx + natom >= len(lines):
            continue
        block = lines[idx + 1 : idx + 1 + natom]
        if not all(atom_line.match(row) for row in block):
            continue
        vals = []
        for row in block:
            parts = row.split()
            vals.extend(float(x) for x in parts[1:4])
        grad = vals
    return grad


def summarize_variant(
    *,
    workdir: Path,
    molecule: str,
    geom: list[tuple[str, float, float, float]],
    variant: str,
    roots: Iterable[int],
    h_ang: float,
    basis: str,
    functional: str,
    nstate: int,
    z_solver: int,
    threads: int,
    command_template: str | None,
    force: bool,
    timeout_s: int,
) -> dict:
    variant_dir = workdir / molecule / variant / f"h_{h_ang:g}"
    natom = len(geom)
    ncomp = natom * 3
    energy_jobs = {}
    for comp in range(ncomp):
        for sign, label in ((1, "plus"), (-1, "minus")):
            name = f"{variant}_e_c{comp:02d}_{label}"
            energy_jobs[(comp, sign)] = run_openqp_job(
                job_dir=variant_dir / name,
                name=name,
                geom=displaced(geom, comp, sign * h_ang),
                runtype="energy",
                target_root=max(roots),
                variant=variant,
                basis=basis,
                functional=functional,
                nstate=nstate,
                z_solver=z_solver,
                threads=threads,
                command_template=command_template,
                force=force,
                timeout_s=timeout_s,
            )

    roots_summary = {}
    component_rows = []
    for root in roots:
        grad_name = f"{variant}_grad_r{root}"
        grad_job = run_openqp_job(
            job_dir=variant_dir / grad_name,
            name=grad_name,
            geom=geom,
            runtype="grad",
            target_root=root,
            variant=variant,
            basis=basis,
            functional=functional,
            nstate=nstate,
            z_solver=z_solver,
            threads=threads,
            command_template=command_template,
            force=force,
            timeout_s=timeout_s,
        )
        root_rows = []
        for comp in range(ncomp):
            eplus = energy_jobs[(comp, 1)].get("energies_ha", {}).get(root)
            eminus = energy_jobs[(comp, -1)].get("energies_ha", {}).get(root)
            fd = central_fd_ha_per_bohr(eplus, eminus, h_ang) if eplus is not None and eminus is not None else None
            analytic = grad_job.get("grad_ha_per_bohr", [None] * ncomp)[comp] if grad_job.get("grad_ha_per_bohr") else None
            row = {
                "molecule": molecule,
                "variant": variant,
                "h_ang": h_ang,
                "root": root,
                "physical_state": physical_state(root),
                "component": comp,
                "atom": comp // 3,
                "axis": "xyz"[comp % 3],
                "analytic_ha_per_bohr": analytic,
                "fd_ha_per_bohr": fd,
                "diff": None if analytic is None or fd is None else analytic - fd,
                "abs_diff": None if analytic is None or fd is None else abs(analytic - fd),
            }
            root_rows.append(row)
            component_rows.append(row)
        valid = [row for row in root_rows if row["abs_diff"] is not None]
        roots_summary[str(root)] = {
            "physical_state": physical_state(root),
            "max_abs_diff": max((row["abs_diff"] for row in valid), default=None),
            "rms_diff": math.sqrt(sum(row["diff"] ** 2 for row in valid) / len(valid)) if valid else None,
            "worst_component": max(valid, key=lambda row: row["abs_diff"]) if valid else None,
            "trah_count": sum(1 for job in list(energy_jobs.values()) + [grad_job] if job.get("trah")),
            "failed_jobs": [job["name"] for job in list(energy_jobs.values()) + [grad_job] if not job.get("ok")],
        }

    summary = {"molecule": molecule, "variant": variant, "h_ang": h_ang, "roots": roots_summary}
    variant_dir.mkdir(parents=True, exist_ok=True)
    (variant_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_component_csv(variant_dir / "component_matrix.csv", component_rows)
    return summary


def write_component_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", default="h2o", help="built-in name or XYZ path")
    parser.add_argument("--workdir", default="mrsf_fd_diagnostic_runs")
    parser.add_argument("--variant", choices=VARIANTS + ("all",), default="default")
    parser.add_argument("--roots", nargs="+", type=int, default=[2, 3])
    parser.add_argument("--steps", nargs="+", type=float, default=[0.002])
    parser.add_argument("--basis", default="3-21g")
    parser.add_argument("--functional", default="bhhlyp")
    parser.add_argument("--nstate", type=int, default=4)
    parser.add_argument("--z-solver", type=int, choices=[0, 1], default=0, help="0=CG, 1=GMRES")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--command", help="OpenQP command template, e.g. 'python -m oqp.pyoqp --nompi {input}'")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args(argv)

    geom = load_geometry(args.molecule)
    molecule = Path(args.molecule).stem if Path(args.molecule).exists() else args.molecule
    variants = VARIANTS if args.variant == "all" else (args.variant,)
    workdir = Path(args.workdir).resolve()
    summaries = []
    for variant in variants:
        for h_ang in args.steps:
            summaries.append(
                summarize_variant(
                    workdir=workdir,
                    molecule=molecule,
                    geom=geom,
                    variant=variant,
                    roots=args.roots,
                    h_ang=h_ang,
                    basis=args.basis,
                    functional=args.functional,
                    nstate=args.nstate,
                    z_solver=args.z_solver,
                    threads=args.threads,
                    command_template=args.command,
                    force=args.force,
                    timeout_s=args.timeout,
                )
            )
    (workdir / "summary.json").write_text(json.dumps(summaries, indent=2))
    print(json.dumps(summaries, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
