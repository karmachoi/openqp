#!/usr/bin/env python3
"""Polyradical / quintet-relevant test systems for QMRSF-icPT2.

TMM (trimethylenemethane, C4H6, D3h): 4 pi orbitals / 4 electrons -> fits CAS(4,4) and the
    quintet (S=2) reference natively. The canonical non-Kekule triplet-ground hydrocarbon.
TME (tetramethyleneethane, C6H8): 6 pi electrons -> genuinely wants CAS(6,6). Run here only
    as an explicitly OUT-OF-(4,4)-SCOPE demonstration (asterisked in the paper): the quintet
    CAS(4,4) treats the 4 frontier electrons, deliberately truncating the 6-pi manifold.

Idealized planar geometries built programmatically (CH2 hydrogens placed automatically).
Quintet ROHF reference (multiplicity=5) -> 4 singly-occupied frontier orbitals.

Run from stageB:
  OPENQP_ROOT=/tmp/qmrsf_root PYTHONPATH=<wt>/pyoqp python3 run_polyrad.py
"""
import os, json
import numpy as np
import run_benchmarks as RB

EV = RB.EV
HERE = RB.HERE
OUT = "/Volumes/External_Storage/claude/sessions/20260620_223000_icpt2_bench_figs"
RCH = 1.08      # C-H
HALF = np.radians(59.0)


def ch2(Cpos, inward_unit):
    """Two H on a planar CH2: bisector points OUTWARD (away from heavy neighbour)."""
    out = -np.asarray(inward_unit, float); out /= np.hypot(*out[:2]) or 1.0
    n = np.array([-out[1], out[0], 0.0])
    H = []
    for s in (+1, -1):
        dirn = np.cos(HALF) * out + s * np.sin(HALF) * n
        H.append(tuple(np.asarray(Cpos) + RCH * dirn))
    return H


def fmt(atoms):
    return "\n".join(" %s  %.6f  %.6f  %.6f" % (s, x, y, z) for s, (x, y, z) in atoms)


def tmm_geom(rcc=1.40):
    C0 = np.array([0.0, 0.0, 0.0])
    angs = [90.0, 210.0, 330.0]
    atoms = [("6", tuple(C0))]
    Hs = []
    for a in angs:
        u = np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])
        Ct = C0 + rcc * u
        atoms.append(("6", tuple(Ct)))
        Hs += ch2(Ct, inward_unit=-u)        # inward = toward centre = -u
    atoms += [("1", h) for h in Hs]
    return fmt(atoms)


def tme_geom(rcc_core=1.47, rcc=1.40):
    C2 = np.array([-rcc_core / 2, 0.0, 0.0]); C3 = np.array([rcc_core / 2, 0.0, 0.0])
    atoms = [("6", tuple(C2)), ("6", tuple(C3))]
    Hs = []
    # C2 substituents at 150 and 210 deg; C3 at 30 and -30 deg
    for C, angs, sgn in ((C2, (150.0, 210.0), -1), (C3, (30.0, -30.0), +1)):
        for a in angs:
            u = np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])
            Ct = C + rcc * u
            atoms.append(("6", tuple(Ct)))
            Hs += ch2(Ct, inward_unit=-u)
    atoms += [("1", h) for h in Hs]
    return fmt(atoms)


def run_sys(name, geom, basis, mult=5):
    tag = "%s_%s" % (name, basis.replace("-", ""))
    # custom template: polyradical multiplicity
    inp = os.path.join(HERE, "poly_%s.inp" % tag)
    open(inp, "w").write(RB.TEMPLATE.format(geom=geom, basis=basis, td="qmrsf_icpt2",
                                            extra="qmrsf_icpt2_shift=0.1").replace(
        "multiplicity=5", "multiplicity=%d" % mult))
    j = RB.run(inp)
    if not j:
        print("  %s: NO JSON (run failed)" % tag); return None
    st = j["states"]
    def sg(ek, mk):
        s = [x[ek] for x in st if x.get(mk) == 1]
        return min(s) if s else min(x[ek] for x in st)
    rec = dict(system=name, basis=basis, ref=j.get("reference_energy"),
               cas=sg("E_CAS", "mult_cas"), en=sg("E_icPT2_EN", "mult_en"),
               dy=sg("E_icPT2_Dyall", "mult_dy"), nstates=len(st))
    print("  %s: ref=%.6f CAS=%.6f EN=%.6f Dy=%.6f (%d states)" %
          (tag, rec["ref"], rec["cas"], rec["en"], rec["dy"], rec["nstates"]))
    return rec


def main():
    res = {}
    print("== TMM (fits CAS(4,4)) ==")
    for basis in ("sto-3g", "6-31g"):
        r = run_sys("TMM", tmm_geom(), basis)
        if r:
            res["TMM/%s" % basis] = r
    print("== TME (asterisked, out-of-(4,4)-scope) ==")
    for basis in ("sto-3g", "6-31g"):
        r = run_sys("TME", tme_geom(), basis)
        if r:
            res["TME/%s" % basis] = r
    json.dump(res, open(os.path.join(OUT, "polyrad_results.json"), "w"), indent=2)
    print("wrote polyrad_results.json")


if __name__ == "__main__":
    main()
