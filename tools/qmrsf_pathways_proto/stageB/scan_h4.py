#!/usr/bin/env python3
"""Geometry SCANS for the QMRSF-icPT2 physics tier (answers the 'single-point toy' criticism).

(A) Square -> rectangle H4 distortion: 4 H on a rectangle of sides a(1+d) x a(1-d).
    d=0 is the D4h square (strongly correlated, near-degenerate) -- the same square->rectangle
    topology as the cyclobutadiene automerization, but exactly FCI-solvable. Shows the icPT2
    S0 surface is SMOOTH through the degeneracy and tracks in-window FCI.
(B) Linear H4 symmetric stretch: z = [0,R,2R,3R], a simple dissociation sanity curve.

For each geometry we run the live OpenQP QMRSF-icPT2 (quintet ROHF reference) and read the
singlet-ground CAS/EN/Dyall totals from the JSON, plus the exact in-window FCI from the live
dump. Serial (each run overwrites the shared live dump, read immediately).

Run from stageB:
  OPENQP_ROOT=/tmp/qmrsf_root PYTHONPATH=<wt>/pyoqp python3 scan_h4.py
"""
import os, json
import numpy as np
import run_benchmarks as RB     # reuse write_inp/run/fci_ground/EV/HERE

EV = RB.EV
HERE = RB.HERE
OUT = "/Volumes/External_Storage/claude/sessions/20260620_223000_icpt2_bench_figs"


def h4_rect(a, d):
    x = 0.5 * a * (1.0 + d); y = 0.5 * a * (1.0 - d)
    return "\n".join(" 1  %.6f  %.6f  0.0" % p for p in
                     [(x, y), (-x, y), (-x, -y), (x, -y)])


def h4_linear(R):
    return "\n".join(" 1  0.0  0.0  %.6f" % (i * R) for i in range(4))


def singlet_ground(states, ek, mk):
    s = [x[ek] for x in states if x.get(mk) == 1]
    return min(s) if s else min(x[ek] for x in states)


def first_singlet_exc(states, ek, mk):
    s = sorted(x[ek] for x in states if x.get(mk) == 1)
    return (s[1] - s[0]) * EV if len(s) > 1 else None


def run_point(tag, geom, basis):
    inp = RB.write_inp(tag, geom, basis, "qmrsf_icpt2", 0.1)
    j = RB.run(inp)
    if not j:
        return None
    st = j["states"]
    rec = dict(
        cas=singlet_ground(st, "E_CAS", "mult_cas"),
        en=singlet_ground(st, "E_icPT2_EN", "mult_en"),
        dy=singlet_ground(st, "E_icPT2_Dyall", "mult_dy"),
        s1_dy=first_singlet_exc(st, "E_icPT2_Dyall", "mult_dy"),
    )
    try:
        rec["fci"] = RB.fci_ground(os.path.join(HERE, "qmrsf_icpt2_full_live.dat"))
    except Exception as e:
        rec["fci"] = None
        print("  fci fail:", e)
    return rec


def scan_square_rect(basis="6-31g", a=1.40, npts=16):
    ds = np.linspace(0.0, 0.45, npts)
    rows = []
    for i, d in enumerate(ds):
        r = run_point("scanSR_%02d" % i, h4_rect(a, d), basis)
        if r:
            r["delta"] = float(d)
            rows.append(r)
            print("  SR d=%.3f  CAS=%.6f EN=%.6f Dy=%.6f FCI=%s" %
                  (d, r["cas"], r["en"], r["dy"], ("%.6f" % r["fci"]) if r["fci"] else "--"))
    return dict(kind="square_rect", basis=basis, a=a, rows=rows)


def scan_linear(basis="6-31g", npts=14):
    Rs = np.linspace(0.90, 2.20, npts)
    rows = []
    for i, R in enumerate(Rs):
        r = run_point("scanLIN_%02d" % i, h4_linear(R), basis)
        if r:
            r["R"] = float(R)
            rows.append(r)
            print("  LIN R=%.3f  CAS=%.6f EN=%.6f Dy=%.6f FCI=%s" %
                  (R, r["cas"], r["en"], r["dy"], ("%.6f" % r["fci"]) if r["fci"] else "--"))
    return dict(kind="linear", basis=basis, rows=rows)


def main():
    print("== square->rectangle H4 scan (6-31G) ==")
    sr = scan_square_rect("6-31g")
    print("== linear H4 stretch scan (6-31G) ==")
    lin = scan_linear("6-31g")
    out = {"square_rect": sr, "linear": lin}
    p = os.path.join(OUT, "scan_results.json")
    json.dump(out, open(p, "w"), indent=2)
    print("wrote", p)


if __name__ == "__main__":
    main()
