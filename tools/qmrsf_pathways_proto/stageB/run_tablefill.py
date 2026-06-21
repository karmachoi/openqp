#!/usr/bin/env python3
"""Fill-in QMRSF-icPT2 runs for the revised manuscript tables:
 tab:ground -> H4/cc-pVDZ, CBD rectangular + CBD square at cc-pVDZ
 tab:poly   -> TMM, TME at cc-pVDZ
Collects ROHF reference, ground-singlet CAS/EN/Dyall totals, and the full <S^2>-labelled
state ladder (for the spin table). Saves table_fill.json.
"""
import os, json
import run_benchmarks as RB
import run_polyrad as RP

EV = RB.EV
OUT = "/Volumes/External_Storage/claude/sessions/20260620_223000_icpt2_bench_figs"

H4 = RB.H4
CBD_RECT = RB.CBD
CBD_SQUARE = """ 6   0.728   0.728   0.0
 6  -0.728   0.728   0.0
 6  -0.728  -0.728   0.0
 6   0.728  -0.728   0.0
 1   1.498   1.498   0.0
 1  -1.498   1.498   0.0
 1  -1.498  -1.498   0.0
 1   1.498  -1.498   0.0"""


def run_cfg(name, geom, basis):
    tag = "tf_%s_%s" % (name, basis.replace("-", ""))
    inp = RB.write_inp(tag, geom, basis, "qmrsf_icpt2", 0.1)
    j = RB.run(inp)
    if not j:
        print("  %s: FAILED" % tag); return None
    st = j["states"]
    def sg(ek, mk):
        s = [x[ek] for x in st if x.get(mk) == 1]
        return min(s) if s else min(x[ek] for x in st)
    rec = dict(system=name, basis=basis, ref=j.get("reference_energy"),
               cas=sg("E_CAS", "mult_cas"), en=sg("E_icPT2_EN", "mult_en"),
               dy=sg("E_icPT2_Dyall", "mult_dy"), states=st)
    print("  %-12s ref=%.5f CAS=%.5f EN=%.5f Dy=%.5f (%d states)" %
          (tag, rec["ref"], rec["cas"], rec["en"], rec["dy"], len(st)))
    return rec


def main():
    res = {}
    cfgs = [
        ("H4", H4, "cc-pvdz"),
        ("CBDrect", CBD_RECT, "cc-pvdz"),
        ("CBDsquare", CBD_SQUARE, "cc-pvdz"),
    ]
    for name, geom, basis in cfgs:
        r = run_cfg(name, geom, basis)
        if r:
            res["%s/%s" % (name, basis)] = r
    # polyradicals at cc-pVDZ
    for name, geom in (("TMM", RP.tmm_geom()), ("TME", RP.tme_geom())):
        r = run_cfg(name, geom, "cc-pvdz")
        if r:
            res["%s/cc-pvdz" % name] = r
    json.dump(res, open(os.path.join(OUT, "table_fill.json"), "w"), indent=2)
    print("wrote table_fill.json")


if __name__ == "__main__":
    main()
