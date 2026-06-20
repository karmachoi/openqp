#!/usr/bin/env python3
"""QMRSF benchmark driver: H4 and CBD, STO-3G and 6-31G, for QMRSF-icPT2 (CAS / EN / Dyall)
and QMRSF-DK. Runs each live, collects the per-state totals from the JSON output, computes an
exact FCI reference for H4 (full det space over the active+virtual window from the live dump),
and emits markdown + LaTeX benchmark tables.

Run from tools/qmrsf_pathways_proto/stageB after building the merged liboqp:
  OPENQP_ROOT=/tmp/qmrsf_root PYTHONPATH=<worktree>/pyoqp python3 run_benchmarks.py
"""
import os, sys, json, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import qmrsf_icpt2_ppp_proto as P
WORKTREE = "/Users/cheolhochoi/Documents/openqp-private-qmrsf-pathways"
PYOQP = os.path.join(WORKTREE, "pyoqp", "oqp", "pyoqp.py")
EV = 27.211386245988

ENV = dict(os.environ)
ENV.setdefault("OPENQP_ROOT", "/tmp/qmrsf_root")
ENV["PYTHONPATH"] = os.path.join(WORKTREE, "pyoqp")

H4 = """ 1   0.0  0.0  0.0
 1   0.0  0.0  1.2
 1   0.0  0.0  2.4
 1   0.0  0.0  3.6"""
CBD = """ 6   0.6735   0.7825   0.0
 6  -0.6735   0.7825   0.0
 6  -0.6735  -0.7825   0.0
 6   0.6735  -0.7825   0.0
 1   1.5235   1.4625   0.0
 1  -1.5235   1.4625   0.0
 1  -1.5235  -1.4625   0.0
 1   1.5235  -1.4625   0.0"""

TEMPLATE = """[input]
system=
{geom}
charge=0
runtype=energy
basis={basis}
method=tdhf

[guess]
type=huckel
save_mol=False

[scf]
save_molden=False
multiplicity=5
type=rohf
converger_type=diis
maxit=200
conv=1.0e-8

[tdhf]
type={td}
{extra}
"""


def write_inp(tag, geom, basis, td, shift):
    extra = "qmrsf_icpt2_shift=%g" % shift if td == "qmrsf_icpt2" else ""
    path = os.path.join(HERE, "bench_%s.inp" % tag)
    open(path, "w").write(TEMPLATE.format(geom=geom, basis=basis, td=td, extra=extra))
    return path


def run(inp):
    log = inp[:-4] + ".log"
    for j in (log[:-4] + ".qmrsf.json", log[:-4] + ".qmrsf_dk.json"):
        if os.path.exists(j):
            os.remove(j)
    subprocess.run([sys.executable, PYOQP, inp], env=ENV, cwd=HERE,
                   stdout=open("/tmp/bench_%s.out" % os.path.basename(inp), "w"),
                   stderr=subprocess.STDOUT, timeout=900)
    for suf in (".qmrsf.json", ".qmrsf_dk.json"):
        j = log[:-4] + suf
        if os.path.exists(j):
            return json.load(open(j))
    return None


def spectrum(states, key):
    tot = sorted(s[key] for s in states)
    g = tot[0]
    return g, [(t - g) * EV for t in tot]


def fci_ground(dump_path):
    """Exact FCI ground (electronic+ecore) over the live window (na=nb=2)."""
    f = open(dump_path); norb, nPd = map(int, f.readline().split())
    h = np.array([[float(x) for x in f.readline().split()] for _ in range(norb)])
    eri = np.zeros((norb,) * 4)
    for a in range(norb):
        for b in range(norb):
            for c in range(norb):
                eri[a, b, c, :] = [float(x) for x in f.readline().split()]
    ecore = float(f.readline())
    H1, g, _ = P.spinorb(h, eri)
    dets = P.gen_dets(norb, 2, 2)
    Hf = P.build_H(dets, H1, g)
    return float(np.linalg.eigvalsh(Hf)[0]) + ecore


def main():
    bench = {}
    for sysname, geom in (("H4", H4), ("CBD", CBD)):
        for basis in ("sto-3g", "6-31g"):
            bkey = basis.replace("-", "")
            jic = run(write_inp("%s_%s_icpt2" % (sysname, bkey), geom, basis, "qmrsf_icpt2", 0.1))
            fci = None
            if sysname == "H4":
                dpath = os.path.join(HERE, "qmrsf_icpt2_full_live.dat")  # written by the icPT2 run
                try:
                    fci = fci_ground(dpath)
                except Exception as e:
                    print("FCI ref failed:", e)
            jdk = run(write_inp("%s_%s_dk" % (sysname, bkey), geom, basis, "qmrsf_dk", 0.0))
            rec = {"system": sysname, "basis": basis, "fci": fci}
            if jic:
                rec["ref"] = jic["reference_energy"]
                rec["cas_g"], rec["cas_exc"] = spectrum(jic["states"], "E_CAS")
                rec["en_g"], rec["en_exc"] = spectrum(jic["states"], "E_icPT2_EN")
                rec["dy_g"], rec["dy_exc"] = spectrum(jic["states"], "E_icPT2_Dyall")
                rec["states"] = jic["states"]
            if jdk:
                rec["dk_g"], rec["dk_exc"] = spectrum(jdk["states"], "E_DK")
            bench["%s/%s" % (sysname, basis)] = rec
            print("done %s/%s" % (sysname, basis))
    json.dump(bench, open(os.path.join(HERE, "benchmark_results.json"), "w"), indent=2)
    emit_tables(bench)


def emit_tables(bench):
    md, tex = [], []
    md.append("# QMRSF benchmark results\n")
    # Table 1: ground-state totals (Ha) + dynamic correlation recovered
    md.append("## Table 1. Ground-state total energy (Hartree)\n")
    md.append("| system/basis | ref ROHF | CAS=DK | icPT2-EN | icPT2-Dyall | FCI | dyn.corr (Dyall) |")
    md.append("|---|---|---|---|---|---|---|")
    for k, r in bench.items():
        fci = "%.6f" % r["fci"] if r.get("fci") is not None else "--"
        dyn = "%.5f" % (r["dy_g"] - r["cas_g"]) if "dy_g" in r else "--"
        md.append("| %s | %.6f | %.6f | %.6f | %.6f | %s | %s |" %
                  (k, r.get("ref", float("nan")), r["cas_g"], r["en_g"], r["dy_g"], fci, dyn))
    # Table 2: lowest excitation energies (eV)
    md.append("\n## Table 2. Lowest vertical excitation energies S0->Sn (eV)\n")
    md.append("| system/basis | state | CAS=DK | icPT2-EN | icPT2-Dyall |")
    md.append("|---|---|---|---|---|")
    for k, r in bench.items():
        for n in (1, 2, 3):
            if "cas_exc" in r and len(r["cas_exc"]) > n:
                md.append("| %s | S0->S%d | %.3f | %.3f | %.3f |" %
                          (k, n, r["cas_exc"][n], r["en_exc"][n], r["dy_exc"][n]))
    # Table 3: spin-resolved lowest excitations (needs S^2 labels)
    def spin_exc(states, ek, mk):
        ss = sorted(states, key=lambda s: s[ek]); g = ss[0][ek]; gm = ss[0].get(mk)
        lo = {}
        for s in ss[1:]:
            m = s.get(mk)
            if m in (1, 3, 5) and m not in lo:
                lo[m] = (s[ek] - g) * EV
        return gm, lo
    SPINMETH = [('E_CAS', 'mult_cas', 'CAS=DK'),
                ('E_icPT2_EN', 'mult_en', 'icPT2-EN'),
                ('E_icPT2_Dyall', 'mult_dy', 'icPT2-Dyall')]
    md.append("\n## Table 3. Spin-resolved lowest excitations (eV): lowest excited singlet / triplet\n")
    md.append("| system/basis | method | ground 2S+1 | S0->S(singlet) | S0->T(triplet) |")
    md.append("|---|---|---|---|---|")
    for k, r in bench.items():
        if 'states' not in r or not r['states'] or r['states'][0].get('mult_cas') is None:
            continue
        for ek, mk, lab in SPINMETH:
            gm, lo = spin_exc(r['states'], ek, mk)
            s1 = ('%.3f' % lo[1]) if 1 in lo else '--'
            t1 = ('%.3f' % lo[3]) if 3 in lo else '--'
            md.append("| %s | %s | %s | %s | %s |" % (k, lab, gm, s1, t1))
    open(os.path.join(HERE, "benchmark_tables.md"), "w").write("\n".join(md) + "\n")

    # LaTeX (Table 1 + Table 2)
    tex.append(r"% QMRSF benchmark tables (auto-generated by run_benchmarks.py)")
    tex.append(r"\begin{table}[t]\centering\footnotesize")
    tex.append(r"\caption{QMRSF ground-state total energies (Hartree). CAS(4,4) equals the "
               r"QMRSF-DK dressed-kernel spectrum (DK==CAS on HF integrals); icPT2 adds the "
               r"external-$Q$ correlation (EN with a $0.1$~Eh imaginary shift, and Dyall). "
               r"For H$_4$ the exact FCI in the same window is shown.}")
    tex.append(r"\begin{tabular}{lrrrrr}\hline")
    tex.append(r"system/basis & ref ROHF & CAS$=$DK & icPT2-EN & icPT2-Dyall & FCI \\ \hline")
    for k, r in bench.items():
        fci = "%.6f" % r["fci"] if r.get("fci") is not None else r"--"
        tex.append("%s & %.6f & %.6f & %.6f & %.6f & %s \\\\" %
                   (k.replace("_", r"\_"), r.get("ref", float("nan")), r["cas_g"], r["en_g"], r["dy_g"], fci))
    tex.append(r"\hline\end{tabular}\label{tab:bench-ground}\end{table}")
    tex.append("")
    tex.append(r"\begin{table}[t]\centering\footnotesize")
    tex.append(r"\caption{QMRSF lowest vertical excitation energies $S_0\!\to\!S_n$ (eV).}")
    tex.append(r"\begin{tabular}{llrrr}\hline")
    tex.append(r"system/basis & state & CAS$=$DK & icPT2-EN & icPT2-Dyall \\ \hline")
    for k, r in bench.items():
        for n in (1, 2, 3):
            if "cas_exc" in r and len(r["cas_exc"]) > n:
                tex.append(r"%s & $S_0\!\to\!S_%d$ & %.3f & %.3f & %.3f \\" %
                           (k.replace("_", r"\_"), n, r["cas_exc"][n], r["en_exc"][n], r["dy_exc"][n]))
    tex.append(r"\hline\end{tabular}\label{tab:bench-exc}\end{table}")
    tex.append("")
    tex.append(r"\begin{table}[t]\centering\footnotesize")
    tex.append(r"\caption{QMRSF spin-resolved lowest excitations (eV): the lowest excited singlet "
               r"and lowest triplet above the ground state, from per-state $\langle\hat S^2\rangle$ "
               r"labels. Ground state is a singlet in every case.}")
    tex.append(r"\begin{tabular}{llrr}\hline")
    tex.append(r"system/basis & method & $S_0\!\to\!S$ (singlet) & $S_0\!\to\!T$ (triplet) \\ \hline")
    for k, r in bench.items():
        if 'states' not in r or not r['states'] or r['states'][0].get('mult_cas') is None:
            continue
        for ek, mk, lab in SPINMETH:
            gm, lo = spin_exc(r['states'], ek, mk)
            s1 = ('%.3f' % lo[1]) if 1 in lo else r'--'
            t1 = ('%.3f' % lo[3]) if 3 in lo else r'--'
            tex.append(r"%s & %s & %s & %s \\" % (k.replace('_', r'\_'), lab, s1, t1))
    tex.append(r"\hline\end{tabular}\label{tab:bench-spin}\end{table}")
    open(os.path.join(HERE, "benchmark_tables.tex"), "w").write("\n".join(tex) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
