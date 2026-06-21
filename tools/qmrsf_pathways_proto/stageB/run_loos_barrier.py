#!/usr/bin/env python3
"""CBD automerization barrier at the REFERENCE (Loos/Monino 2022) optimized geometries,
cc-pVDZ and cc-pVTZ, for QMRSF-icPT2 (CAS/EN/Dyall). Geometries reconstructed from the
CASPT2(12,12)/aug-cc-pVTZ internal coordinates (Table 1 of Monino et al. 2022):
  D2h (1^1Ag): C-C 1.354 / 1.566, C-H 1.076 A, ring 90deg -> H-C-C 135deg
  D4h (1^1B1g): C-C 1.449, C-H 1.076 A
This makes our barrier directly comparable to the literature (same geometry, basis-converged).
"""
import os, json, math
import run_benchmarks as RB

OUT = "/Volumes/External_Storage/claude/sessions/20260620_223000_icpt2_bench_figs"
CH = 1.076
B = CH * math.cos(math.radians(45.0))   # H offset along each axis (external bisector, 135deg)


def cbd_geom(cc_short, cc_long):
    xc, yc = cc_short / 2.0, cc_long / 2.0
    C = [(xc, yc), (-xc, yc), (-xc, -yc), (xc, -yc)]
    lines = []
    for (x, y) in C:
        lines.append(" 6  %.6f  %.6f  0.0" % (x, y))
    for (x, y) in C:
        hx = x + math.copysign(B, x); hy = y + math.copysign(B, y)
        lines.append(" 1  %.6f  %.6f  0.0" % (hx, hy))
    return "\n".join(lines)


D2H = cbd_geom(1.354, 1.566)   # rectangle minimum
D4H = cbd_geom(1.449, 1.449)   # square TS
KCAL = 627.509474


def run(name, geom, basis):
    tag = "loos_%s_%s" % (name, basis.replace("-", ""))
    j = RB.run(RB.write_inp(tag, geom, basis, "qmrsf_icpt2", 0.1))
    if not j:
        print("  %s FAILED" % tag); return None
    st = j["states"]
    def sg(ek, mk):
        s = [x[ek] for x in st if x.get(mk) == 1]
        return min(s) if s else min(x[ek] for x in st)
    rec = dict(cas=sg("E_CAS", "mult_cas"), en=sg("E_icPT2_EN", "mult_en"),
               dy=sg("E_icPT2_Dyall", "mult_dy"), ref=j.get("reference_energy"))
    print("  %s: CAS=%.6f EN=%.6f Dyall=%.6f" % (tag, rec["cas"], rec["en"], rec["dy"]), flush=True)
    return rec


def main():
    out = {}
    for basis in ("cc-pvdz", "cc-pvtz"):
        d2 = run("D2h", D2H, basis)
        d4 = run("D4h", D4H, basis)
        if d2 and d4:
            bar = {k: (d4[k] - d2[k]) * KCAL for k in ("cas", "en", "dy")}
            out[basis] = dict(D2h=d2, D4h=d4, barrier_kcal=bar)
            print("  >> %s barrier (kcal/mol): CAS=%+.2f EN=%+.2f Dyall=%+.2f"
                  % (basis, bar["cas"], bar["en"], bar["dy"]), flush=True)
    json.dump(out, open(os.path.join(OUT, "loos_barrier.json"), "w"), indent=2)
    print("wrote loos_barrier.json")


if __name__ == "__main__":
    main()
