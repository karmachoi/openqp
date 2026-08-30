#!/usr/bin/env python3
"""Generate matched OpenQP uracil NAMD inputs from the archived NX conditions."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


BOHR_TO_ANGSTROM = 0.529177210903
SOURCE_COMMIT = "deb378d6e062495fddda0b4391f4d7821b7f6c1b"
METHODS = {
    "baseline_npi_iso": {
        "tdc": "npi",
        "rescale": "isotropic",
        "environment": {},
    },
    "full_analytic": {
        "tdc": "analytic",
        "rescale": "analytic_nac",
        "environment": {},
    },
    "rlzt10": {
        "tdc": "analytic",
        "rescale": "analytic_nac",
        "environment": {
            "OQP_MRSF_NAC_ZV_PREDICTOR": "linear_approx",
            "OQP_MRSF_NAC_ZV_ETA": "1.0",
            "OQP_MRSF_NAC_ZV_MAX_DISP": "0.25",
            "OQP_MRSF_NAC_ZV_WARMUP_EXACT": "2",
            "OQP_MRSF_NAC_ZV_EXACT_EVERY": "10",
        },
    },
    "ht_nac": {
        "tdc": "npi",
        "rescale": "hop_analytic_nac",
        "environment": {},
    },
}
PHASE_STEPS = {"smoke": 1, "pilot": 100, "production": 1000}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_nx(path: Path, wanted: set[int]) -> dict[int, dict[str, object]]:
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    records: dict[int, dict[str, object]] = {}
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        if not stripped.startswith("Initial condition ="):
            index += 1
            continue
        condition = int(stripped.split("=")[-1])
        index += 1
        if condition not in wanted:
            continue
        while index < len(lines) and "Geometry in COLUMBUS" not in lines[index]:
            index += 1
        index += 1
        geometry = []
        while index < len(lines) and len(geometry) < 12:
            fields = lines[index].split()
            index += 1
            if len(fields) >= 6 and fields[0][0].isalpha():
                geometry.append((fields[0], *(float(value) for value in fields[2:5])))
        while index < len(lines) and "Velocity in NX input format" not in lines[index]:
            index += 1
        index += 1
        velocity = []
        while index < len(lines) and len(velocity) < 12:
            fields = lines[index].split()
            index += 1
            if len(fields) == 3:
                velocity.append(tuple(float(value) for value in fields))
        if len(geometry) != 12 or len(velocity) != 12:
            raise ValueError(f"initial condition {condition} is incomplete")
        records[condition] = {"geometry_bohr": geometry, "velocity_au": velocity}
    missing = sorted(wanted - records.keys())
    if missing:
        raise ValueError(f"missing initial conditions: {missing}")
    return records


def input_text(record: dict[str, object], active: int, method: dict[str, object],
               nstep: int, stream: int) -> str:
    geometry = "\n".join(
        f"   {symbol:<2s} {x * BOHR_TO_ANGSTROM: .12f} "
        f"{y * BOHR_TO_ANGSTROM: .12f} {z * BOHR_TO_ANGSTROM: .12f}"
        for symbol, x, y, z in record["geometry_bohr"]
    )
    return f"""[input]
system=
{geometry}
charge=0
runtype=namd
basis=6-31g*
functional=bhhlyp
method=tdhf
ispher=auto
perf=1
omp_threads=1

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf
maxit=200
conv=1e-10

[tdhf]
type=mrsf
nstate=4
multiplicity=1
conv=1e-10

[properties]
grad={active}

[md]
nstep={nstep}
dt=0.5
active={active}
substep=50000
init_temp=300
velocity=velocity.au
seed=20260830
rng_stream={stream}
first_hop_step=1
tdc={method['tdc']}
rescale={method['rescale']}
nacme_check=off
thrshe=0.36749322175655
decoherence=edc
edc_c=0.1
trivial=False
ensemble=nve
thermostat=off
nve_gate=warn
trajectory_interval=1
restart_interval=20
trajectory_file=trajectory.namd.trj
restart_file=restart.npz
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--phase", choices=PHASE_STEPS, required=True)
    parser.add_argument("--method", choices=METHODS, action="append")
    args = parser.parse_args()

    selection = json.loads(args.selection.read_text(encoding="utf-8"))
    s1 = set(selection["initial_s1"])
    s2 = set(selection["initial_s2"])
    if args.phase == "pilot":
        chosen = set(selection["pilot"])
    elif args.phase == "smoke":
        chosen = {selection["pilot"][0]}
    else:
        chosen = s1 | s2
    methods = args.method or list(METHODS)
    records = parse_nx(args.source, chosen)
    source_hash = sha256(args.source)
    selection_hash = sha256(args.selection)
    manifest = {
        "schema": "openqp-uracil-namd-input-v1",
        "source_commit": SOURCE_COMMIT,
        "phase": args.phase,
        "nstep": PHASE_STEPS[args.phase],
        "dt_fs": 0.5,
        "electronic_substeps": 50000,
        "source_initial_conditions": str(args.source.resolve()),
        "source_sha256": source_hash,
        "selection_sha256": selection_hash,
        "methods": METHODS,
        "calculations": [],
    }
    job_rows = []
    for method_name in methods:
        method = METHODS[method_name]
        for condition in sorted(chosen):
            active = 2 if condition in s1 else 3
            calc_root = args.output / args.phase / method_name / f"ic-{condition:04d}"
            calc_root.mkdir(parents=True, exist_ok=False)
            velocity_path = calc_root / "velocity.au"
            velocity_path.write_text("".join(
                f"{x: .12e} {y: .12e} {z: .12e}\n"
                for x, y, z in records[condition]["velocity_au"]
            ), encoding="ascii")
            input_path = calc_root / "uracil.inp"
            input_path.write_text(input_text(
                records[condition], active, method,
                PHASE_STEPS[args.phase], condition), encoding="ascii")
            raw_path = calc_root / "initial-condition.json"
            raw_path.write_text(json.dumps({
                "condition": condition,
                "initial_state": "S1" if condition in s1 else "S2",
                "active_root": active,
                **records[condition],
            }, indent=2) + "\n", encoding="ascii")
            environment_path = calc_root / "run.env"
            environment_path.write_text("".join(
                f"{key}={value}\n"
                for key, value in sorted(method["environment"].items())
            ), encoding="ascii")
            relative_root = calc_root.relative_to(args.output)
            array_index = len(job_rows)
            job_rows.append(
                f"{array_index}\t{args.phase}:{method_name}:ic-{condition:04d}"
                f"\t{method_name}\t{relative_root}\n"
            )
            manifest["calculations"].append({
                "identity": f"{args.phase}:{method_name}:ic-{condition:04d}",
                "method": method_name,
                "condition": condition,
                "initial_state": "S1" if condition in s1 else "S2",
                "active_root": active,
                "environment": method["environment"],
                "input": str(input_path.relative_to(args.output)),
                "input_sha256": sha256(input_path),
                "velocity_sha256": sha256(velocity_path),
                "initial_condition_sha256": sha256(raw_path),
                "environment_sha256": sha256(environment_path),
            })
    jobs_path = args.output / f"{args.phase}-jobs.tsv"
    jobs_path.write_text("".join(job_rows), encoding="ascii")
    manifest["jobs_tsv"] = str(jobs_path.relative_to(args.output))
    manifest["jobs_tsv_sha256"] = sha256(jobs_path)
    manifest_path = args.output / f"{args.phase}-input-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="ascii")
    print(manifest_path)
    print(sha256(manifest_path))


if __name__ == "__main__":
    main()
