#!/usr/bin/env python3
"""Scalability tier for QMRSF-icPT2.

For each system/basis: the number of STREAMED external contracted perturbers N_pert
(validated closed form count_model, = manuscript's 4,064,220 for CBD/cc-pVDZ), the size of
the in-window full-CI determinant list that is NEVER materialized, and the measured live
wall-time + peak resident memory. The selling point: peak RAM tracks the orbital window
(nvirt), NOT the perturber/determinant count.

Run from stageB:
  OPENQP_ROOT=/tmp/qmrsf_root PYTHONPATH=<wt>/pyoqp python3 run_scalability.py
"""
import os, sys, json, subprocess, re, time
from math import comb
import run_benchmarks as RB

HERE = RB.HERE
PYOQP = RB.PYOQP
ENV = RB.ENV
OUT = "/Volumes/External_Storage/claude/sessions/20260620_223000_icpt2_bench_figs"
NACT = 4
TIME = "/usr/bin/time"

# nbf per (system, basis): per-element contracted-function counts
NBF_EL = {
    "sto-3g": {"H": 1, "C": 5},
    "6-31g":  {"H": 2, "C": 9},
    "cc-pvdz": {"H": 5, "C": 14},
}
SYS = {  # (n_H, n_C, n_electrons)
    "H4":  (4, 0, 4),
    "CBD": (4, 4, 28),
}


def count_model(nvirt, na=2, nb=2):
    def sc(nsig):
        return sum(comb(NACT, nsig - nv) * comb(nvirt, nv)
                   for nv in range(0, nsig + 1) if nsig - nv <= NACT and nv <= nvirt)
    return sc(na) * sc(nb) - comb(NACT, na) * comb(NACT, nb)


def nbf_of(system, basis):
    nH, nC, _ = SYS[system]
    t = NBF_EL[basis]
    return nH * t["H"] + nC * t["C"]


def timed_run(tag, geom, basis):
    inp = RB.write_inp(tag, geom, basis, "qmrsf_icpt2", 0.1)
    # clear stale json so we can confirm success
    for suf in (".qmrsf.json", ".qmrsf_dk.json"):
        j = inp[:-4] + suf
        if os.path.exists(j):
            os.remove(j)
    cmd = [TIME, "-l", sys.executable, PYOQP, inp]
    t0 = time.perf_counter()
    p = subprocess.run(cmd, env=ENV, cwd=HERE, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE, timeout=3600)
    wall = time.perf_counter() - t0
    err = p.stderr.decode()
    m = re.search(r"(\d+)\s+maximum resident set size", err)
    rss = int(m.group(1)) if m else None     # bytes on macOS
    ok = os.path.exists(inp[:-4] + ".qmrsf.json")
    # read live dump for the in-window orbital count (rigorous, tied to this run)
    norb = None
    try:
        with open(os.path.join(HERE, "qmrsf_icpt2_full_live.dat")) as f:
            norb = int(f.readline().split()[0])
    except Exception:
        pass
    return dict(ok=ok, wall=wall, rss=rss, norb=norb)


GEOM = {"H4": RB.H4, "CBD": RB.CBD}


def main():
    rows = []
    for system in ("H4", "CBD"):
        for basis in ("sto-3g", "6-31g", "cc-pvdz"):
            nbf = nbf_of(system, basis)
            _, _, nel = SYS[system]
            ncore = (nel - 4) // 2
            nvirt_an = nbf - ncore - NACT
            npert = count_model(nvirt_an)
            r = timed_run("scal_%s_%s" % (system, basis.replace("-", "")), GEOM[system], basis)
            nvirt_live = (r["norb"] - NACT) if r["norb"] else None
            fci_list = comb(nbf - ncore, 2) ** 2     # in-window FCI det list (na=nb=2)
            rec = dict(system=system, basis=basis, nbf=nbf, ncore=ncore,
                       nvirt=nvirt_an, nvirt_live=nvirt_live, n_pert=npert,
                       fci_window_dets=fci_list, wall_s=r["wall"],
                       peak_mb=(r["rss"] / 1e6 if r["rss"] else None), ok=r["ok"])
            rows.append(rec)
            print("  %-4s/%-8s nbf=%-3d nvirt=%-3d N_pert=%-10d  wall=%6.2fs  "
                  "peakMB=%-7s ok=%s  (nvirt_live=%s)" %
                  (system, basis, nbf, nvirt_an, npert, r["wall"],
                   ("%.0f" % (r["rss"]/1e6)) if r["rss"] else "NA", r["ok"], nvirt_live))
    json.dump(rows, open(os.path.join(OUT, "scalability_results.json"), "w"), indent=2)
    print("wrote scalability_results.json")


if __name__ == "__main__":
    main()
