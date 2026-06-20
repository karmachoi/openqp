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
{func}method=tdhf

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


def write_inp(tag, geom, basis, td, shift, functional=None):
    extra = "qmrsf_icpt2_shift=%g" % shift if td == "qmrsf_icpt2" else ""
    # functional=bhhlyp -> ROKS reference + the genuine DFT-dressed grid kernel
    func = "functional=%s\n" % functional if functional else ""
    path = os.path.join(HERE, "bench_%s.inp" % tag)
    open(path, "w").write(TEMPLATE.format(geom=geom, basis=basis, td=td, extra=extra, func=func))
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
            # genuine grid-derived kernel (DFT-dressed: adiabatic f_xc + transverse f^{+-}).
            # Requires a populated beta channel (nelec_B>0); the quintet H4 reference is
            # fully spin-polarized (nelec_B=0, M_s=+2) so the spin-flip-down kernels are
            # ill-defined there -- run only for the cored diradical (CBD).
            jdkg = None
            if sysname == "CBD":
                jdkg = run(write_inp("%s_%s_dkg" % (sysname, bkey), geom, basis, "qmrsf_dk", 0.0,
                                     functional="bhhlyp"))
            rec = {"system": sysname, "basis": basis, "fci": fci}
            if jic:
                rec["ref"] = jic["reference_energy"]
                rec["cas_g"], rec["cas_exc"] = spectrum(jic["states"], "E_CAS")
                rec["en_g"], rec["en_exc"] = spectrum(jic["states"], "E_icPT2_EN")
                rec["dy_g"], rec["dy_exc"] = spectrum(jic["states"], "E_icPT2_Dyall")
                rec["states"] = jic["states"]
            if jdk:
                rec["dk_g"], rec["dk_exc"] = spectrum(jdk["states"], "E_DK")
            # genuine-kernel DFT-DK states carry per-state mult + E_DK_DFT (grid kernel)
            if jdkg and jdkg.get("is_dft_dressed"):
                rec["dkg_states"] = jdkg["states"]
            bench["%s/%s" % (sysname, basis)] = rec
            print("done %s/%s" % (sysname, basis))
    json.dump(bench, open(os.path.join(HERE, "benchmark_results.json"), "w"), indent=2)
    emit_tables(bench)


SPINMETH = [('E_CAS', 'mult_cas', 'CAS=DK'),
            ('E_icPT2_EN', 'mult_en', 'icPT2-EN'),
            ('E_icPT2_Dyall', 'mult_dy', 'icPT2-Dyall')]


def singlet_ladder(states, ek, mk, nmax=3):
    """Excitations (eV) to the lowest nmax SINGLET states above the singlet ground.

    Spin-matched: filters to <S^2>-labelled singlets (mk==1) BEFORE ranking, so the
    n-th entry is genuinely the n-th singlet -- not whatever root happens to sit n-th
    in the all-multiplicity energy order (which is the state-mismatch artifact)."""
    ss = sorted(states, key=lambda s: s[ek])
    sing = [s for s in ss if s.get(mk) == 1]
    if not sing:
        return []
    g = sing[0][ek]
    return [(s[ek] - g) * EV for s in sing[1:1 + nmax]]


def lowest_mult(states, ek, mk, m):
    """Lowest excitation (eV) of multiplicity m, measured from the singlet ground."""
    ss = sorted(states, key=lambda s: s[ek])
    sing = [s for s in ss if s.get(mk) == 1]
    if not sing:
        return None
    g = sing[0][ek]
    for s in ss:
        if s.get(mk) == m and s[ek] > g + 1e-9:
            return (s[ek] - g) * EV
    return None


def has_spin(r):
    return ('states' in r and r['states'] and r['states'][0].get('mult_cas') is not None)


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
    # Table 2: SINGLET-ONLY excitations (spin-matched via <S^2>; no all-multiplicity mixing)
    md.append("\n## Table 2. Lowest vertical SINGLET excitation energies S0->Sn (eV)\n")
    md.append("_Spin-matched: each S_n is the n-th <S^2>-labelled singlet. "
              "DK==CAS on HF integrals, including the spin labels (GATE 1)._\n")
    md.append("_DK-DFT(grid) = the genuine grid-derived kernel (adiabatic f_xc + transverse "
              "f^{+-}) on a BHHLYP/ROKS reference; bare DK==CAS on HF integrals._\n")
    md.append("| system/basis | state | CAS=DK | DK-DFT(grid) | icPT2-EN | icPT2-Dyall |")
    md.append("|---|---|---|---|---|---|")
    for k, r in bench.items():
        if not has_spin(r):
            continue
        lc = singlet_ladder(r['states'], 'E_CAS', 'mult_cas')
        le = singlet_ladder(r['states'], 'E_icPT2_EN', 'mult_en')
        ld = singlet_ladder(r['states'], 'E_icPT2_Dyall', 'mult_dy')
        lg = singlet_ladder(r['dkg_states'], 'E_DK_DFT', 'mult') if r.get('dkg_states') else []
        for n in range(1, 4):
            if len(lc) >= n and len(le) >= n and len(ld) >= n:
                gstr = '%.3f' % lg[n - 1] if len(lg) >= n else '--'
                md.append("| %s | S0->S%d | %.3f | %s | %.3f | %.3f |" %
                          (k, n, lc[n - 1], gstr, le[n - 1], ld[n - 1]))
    # Table 3: the artifact made explicit -- lowest TRIPLET T1 vs lowest SINGLET S1
    md.append("\n## Table 3. Spin-resolved: lowest triplet T1 vs lowest singlet S1 (eV)\n")
    md.append("_The naive 'state 1' (lowest root above ground) is the TRIPLET; comparing it "
              "to a singlet-only column overstates the DK--icPT2 gap. T1 < S1 in every row._\n")
    md.append("| system/basis | method | ground 2S+1 | S0->T1 (triplet) | S0->S1 (singlet) |")
    md.append("|---|---|---|---|---|")
    for k, r in bench.items():
        if not has_spin(r):
            continue
        for ek, mk, lab in SPINMETH:
            gm = sorted(r['states'], key=lambda s: s[ek])[0].get(mk)
            t1 = lowest_mult(r['states'], ek, mk, 3)
            s1l = singlet_ladder(r['states'], ek, mk, 1)
            t1s = '%.3f' % t1 if t1 is not None else '--'
            s1s = '%.3f' % s1l[0] if s1l else '--'
            md.append("| %s | %s | %s | %s | %s |" % (k, lab, gm, t1s, s1s))
    open(os.path.join(HERE, "benchmark_tables.md"), "w").write("\n".join(md) + "\n")

    # LaTeX (Table 1 + Table 2 singlet-only + Table 3 spin-resolved)
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
    tex.append(r"\caption{QMRSF lowest vertical \emph{singlet} excitation energies "
               r"$S_0\!\to\!S_n$ (eV). Spin-matched: each $S_n$ is the $n$-th "
               r"$\langle\hat S^2\rangle$-labelled singlet, so the columns compare like with "
               r"like. DK$=$CAS on HF integrals (including the spin labels).}")
    tex.append(r"\begin{tabular}{llrrrr}\hline")
    tex.append(r"system/basis & state & CAS$=$DK & DK-DFT(grid) & icPT2-EN & icPT2-Dyall \\ \hline")
    for k, r in bench.items():
        if not has_spin(r):
            continue
        lc = singlet_ladder(r['states'], 'E_CAS', 'mult_cas')
        le = singlet_ladder(r['states'], 'E_icPT2_EN', 'mult_en')
        ld = singlet_ladder(r['states'], 'E_icPT2_Dyall', 'mult_dy')
        lg = singlet_ladder(r['dkg_states'], 'E_DK_DFT', 'mult') if r.get('dkg_states') else []
        for n in range(1, 4):
            if len(lc) >= n and len(le) >= n and len(ld) >= n:
                gstr = '%.3f' % lg[n - 1] if len(lg) >= n else r'--'
                tex.append(r"%s & $S_0\!\to\!S_%d$ & %.3f & %s & %.3f & %.3f \\" %
                           (k.replace("_", r"\_"), n, lc[n - 1], gstr, le[n - 1], ld[n - 1]))
    tex.append(r"\hline\end{tabular}\label{tab:bench-exc}\end{table}")
    tex.append("")
    tex.append(r"\begin{table}[t]\centering\footnotesize")
    tex.append(r"\caption{Spin-resolved lowest excitations (eV): the lowest triplet $T_1$ and "
               r"lowest excited singlet $S_1$ above the (singlet) ground state, from per-state "
               r"$\langle\hat S^2\rangle$ labels. $T_1<S_1$ in every case, so a naive "
               r"energy-ordered ``state~1'' is the triplet, not $S_1$.}")
    tex.append(r"\begin{tabular}{llrr}\hline")
    tex.append(r"system/basis & method & $S_0\!\to\!T_1$ & $S_0\!\to\!S_1$ \\ \hline")
    for k, r in bench.items():
        if not has_spin(r):
            continue
        for ek, mk, lab in SPINMETH:
            t1 = lowest_mult(r['states'], ek, mk, 3)
            s1l = singlet_ladder(r['states'], ek, mk, 1)
            t1s = '%.3f' % t1 if t1 is not None else r'--'
            s1s = '%.3f' % s1l[0] if s1l else r'--'
            tex.append(r"%s & %s & %s & %s \\" % (k.replace('_', r'\_'), lab, t1s, s1s))
    tex.append(r"\hline\end{tabular}\label{tab:bench-spin}\end{table}")
    open(os.path.join(HERE, "benchmark_tables.tex"), "w").write("\n".join(tex) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
