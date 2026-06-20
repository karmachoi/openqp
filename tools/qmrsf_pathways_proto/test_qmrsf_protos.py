#!/usr/bin/env python3
"""
Regression harness for the QMRSF pathway prototypes (pure NumPy, NO pytest / NO scipy).

Guards the load-bearing VALIDATED NUMBERS of four prototypes so future edits cannot
silently break the science. Each module is imported (its main() is NOT run -- every
proto guards main() under `if __name__ == "__main__"`), and the key quantities are
recomputed from the LOW-LEVEL public functions and asserted within tolerances chosen
tight enough to catch regressions but loose enough to pass the current validated values.

Modules guarded:
  qmrsf_icpt2_ppp_proto.py  -- H Hermitian; U=0 FCI == sum of lowest orbital energies;
                               external-Q EN downfold recovers a sizable fraction of the
                               CAS->FCI gap for the default model.
  qmrsf_dk_proto.py         -- adiabatic MISSES the double (error > 0.1); dressed pole
                               search matches exact to < 1e-8.
  qmrsf_dk_block_proto.py   -- P1 recovers the injected 0OS double pole to high accuracy
                               in the well-separated case; P2 anchor matches exact < 1e-8.
  qmrsf_icpt2_multistate.py -- H_eff Hermitian (< 1e-10); avoided-crossing smoothness
                               (max|2nd diff|) small / comparable to FCI.

Run:  python3 test_qmrsf_protos.py   (exit code 0 = all pass)
"""
import os
import sys
import traceback

import numpy as np

# Make sure the protos (same directory) are importable regardless of CWD.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Import the modules WITHOUT triggering their main() (all guarded by __name__ check).
import qmrsf_icpt2_ppp_proto as ppp
import qmrsf_dk_proto as dk
import qmrsf_dk_block_proto as dkb
import qmrsf_icpt2_multistate as ms


# ----------------------------------------------------------------------
# Tiny test runner (plain asserts, no pytest)
# ----------------------------------------------------------------------
_RESULTS = []


def check(name, fn):
    """Run one check fn(); record PASS/FAIL; print one line; never abort the suite."""
    try:
        detail = fn()
        _RESULTS.append((name, True, detail))
        print(f"PASS: {name}" + (f"   [{detail}]" if detail else ""))
    except AssertionError as exc:
        _RESULTS.append((name, False, str(exc)))
        print(f"FAIL: {name}   -- {exc}")
    except Exception as exc:  # unexpected error: also a failure, with traceback
        _RESULTS.append((name, False, f"{type(exc).__name__}: {exc}"))
        print(f"FAIL: {name}   -- UNEXPECTED {type(exc).__name__}: {exc}")
        traceback.print_exc()


def approx(name, got, want, tol):
    assert abs(got - want) < tol, (
        f"{name}: |{got:.6g} - {want:.6g}| = {abs(got - want):.3e} >= tol {tol:.1e}")


def le(name, got, bound):
    assert got <= bound, f"{name}: {got:.3e} > bound {bound:.3e}"


def ge(name, got, bound):
    assert got >= bound, f"{name}: {got:.3e} < bound {bound:.3e}"


def in_band(name, got, lo, hi):
    assert lo <= got <= hi, f"{name}: {got:.4f} not in [{lo}, {hi}]"


# ======================================================================
# 1. qmrsf_icpt2_ppp_proto.py
# ======================================================================
def t_ppp_hermitian():
    # Default model: hexatriene PPP, CAS(4,4), thop=1.0 (proto's run_case).
    h_mo, eri_mo, eps, dets, Hfull, Pidx, part = ppp.run_case(6, 6, 1.0)
    res = np.abs(Hfull - Hfull.T).max()
    le("max|H-H^T|", res, 1e-10)
    return f"max|H-H^T|={res:.2e}, dets={len(dets)}, P={len(Pidx)}"


def t_ppp_u0_fci():
    # U=0 limit: zero the 2e tensor -> FCI ground == sum of lowest nelec spin-orbital
    # energies. Validates the 1e Slater-Condon signs / assembly.
    h_mo, eri_mo, eps, dets, Hfull, Pidx, part = ppp.run_case(6, 6, 1.0)
    H1z, gz, _ = ppp.spinorb(h_mo, np.zeros_like(eri_mo))
    Hz = ppp.build_H(dets, H1z, gz)
    e0_ci = np.linalg.eigvalsh(Hz)[0]
    so_eps = np.sort(np.concatenate([eps, eps]))
    e0_ref = so_eps[:6].sum()
    approx("U=0 FCI vs sum-eps", e0_ci, e0_ref, 1e-10)
    return f"FCI={e0_ci:.6f}, sum-eps={e0_ref:.6f}, diff={abs(e0_ci - e0_ref):.2e}"


def t_ppp_recovery_band():
    # External-Q internally-contracted EN downfold recovers a sizable fraction of the
    # CAS->FCI dynamic-correlation gap for the DEFAULT model (thop=1.0).
    # Current validated value ~0.855 (85.5%); guard the sane band [0.55, 0.98].
    h_mo, eri_mo, eps, dets, Hfull, Pidx, part = ppp.run_case(6, 6, 1.0)
    e_fci = np.linalg.eigvalsh(Hfull)[0]
    e_cas, e_pt2, sig = ppp.icpt2(Hfull, dets, Pidx, nroots=1)[0]
    gap = e_cas - e_fci
    assert abs(gap) > 1e-9, "degenerate CAS==FCI gap; nothing to recover"
    recov = (e_cas - e_pt2) / gap
    in_band("EN recovery fraction", recov, 0.55, 0.98)
    return f"recovery={recov:.4f} (CAS={e_cas:.4f}, PT2={e_pt2:.4f}, FCI={e_fci:.4f})"


# ======================================================================
# 2. qmrsf_dk_proto.py
# ======================================================================
def t_dk_adiabatic_misses_double():
    # Canonical MCZB model (Case 1): A=5.0, omega_d=5.6, V=0.8.
    # ADIABATIC sits exactly at A and MUST mis-place the single (and miss the double):
    # error vs nearest exact root > 0.1. Current validated value ~0.554.
    A, wd, V = 5.0, 5.6, 0.8
    exact = dk.exact_spectrum(A, wd, V)
    adiab = dk.adiabatic_energy(A, wd, V)
    # error of adiabatic vs the NEAREST exact root (it still misses the other root).
    err = np.abs(exact - adiab).min()
    ge("adiabatic error vs nearest exact", err, 0.1)
    # also confirm a double is structurally absent: adiabatic gives 1 root, exact gives 2.
    assert exact.size == 2, f"expected 2 exact roots, got {exact.size}"
    return f"adiabatic err={err:.4f} (>0.1), exact roots={exact.size}"


def t_dk_dressed_matches_exact():
    # DRESSED pole search must recover the FULL single+double spectrum to < 1e-8.
    # Guard every model the proto validates (Cases 1-3 + the two 0OS-like channels).
    cases = [
        ("Case 1", 5.0, np.array([5.6]), np.array([0.8])),
        ("Case 2", 4.0, np.array([4.05]), np.array([0.5])),
        ("Case 3", 6.0, np.array([4.5, 6.4, 8.1]), np.array([0.45, 0.7, 0.3])),
        ("0OS #1", 5.2, np.array([4.9, 6.1]), np.array([0.6, 0.4])),
        ("0OS #2", 7.0, np.array([6.7, 8.3]), np.array([0.5, 0.55])),
    ]
    worst = 0.0
    for nm, A, wd, V in cases:
        exact = dk.exact_spectrum(A, wd, V)
        dressed = dk.dressed_roots(A, wd, V)
        assert dressed.size == exact.size, (
            f"{nm}: root count {dressed.size} != exact {exact.size}")
        err = np.abs(np.sort(dressed) - np.sort(exact)).max()
        le(f"{nm} dressed vs exact", err, 1e-8)
        worst = max(worst, err)
    return f"max|dressed-exact| over 5 models = {worst:.2e} (<1e-8)"


# ======================================================================
# 3. qmrsf_dk_block_proto.py
# ======================================================================
def _injected_double(p1_spec, Omega_ad):
    """The injected double-like root = the P1 root farthest from any adiabatic single."""
    return p1_spec[int(np.argmax([min(abs(r - Omega_ad)) for r in p1_spec]))]


def t_dkb_p1_recovers_double_abstract():
    # Well-separated abstract block (proto case i.a): 3 coupled singles + 1 isolated
    # double at 7.6. P1 must recover the injected double pole to high accuracy.
    # Current validated: injected double matches exact to ~1e-13.
    H_ss = np.array([[5.0, 0.25, 0.00],
                     [0.25, 5.8, 0.20],
                     [0.00, 0.20, 6.7]])
    omega_d = np.array([7.6])
    W = np.array([[0.5], [0.4], [0.3]])
    exact = dkb.exact_block_spectrum(H_ss, omega_d, W)
    Omega_ad, _ = dkb.adiabatic_backbone(H_ss)
    p1, _, _ = dkb.p1_full_spectrum(H_ss, omega_d, W)
    extra = _injected_double(p1, Omega_ad)
    dbl_exact = exact[int(np.argmin(np.abs(exact - extra)))]
    err = abs(extra - dbl_exact)
    le("i.a injected double err", err, 1e-6)
    return f"injected double P1={extra:.6f} vs exact={dbl_exact:.6f}, err={err:.2e}"


def t_dkb_p1_recovers_double_ppp():
    # Real PPP CAS(4,4) ii-B: the lowest, well-separated 0OS-type double dressing the
    # singles backbone. P1 must INJECT that double pole the adiabatic backbone cannot
    # produce. Current validated: P1=71.5183 vs exact=71.5183 (~1e-13).
    H_cas, dets, ref, singles, doubles, higher = dkb.ppp_cas44_block()
    Ns = len(singles)
    sub = singles + doubles
    Hsub = H_cas[np.ix_(sub, sub)]
    sl = list(range(Ns))
    dl = list(range(Ns, len(sub)))
    H_ss = Hsub[np.ix_(sl, sl)]
    W = Hsub[np.ix_(sl, dl)]
    omega_d = np.diag(Hsub[np.ix_(dl, dl)])
    dcol = int(np.argmin(omega_d))                       # lowest, cleanest satellite
    Wf = W[:, [dcol]]
    wdf = omega_d[[dcol]]
    exact_f = dkb.exact_block_spectrum(H_ss, wdf, Wf)
    Om_f, _ = dkb.adiabatic_backbone(H_ss)
    p1_f, _, _ = dkb.p1_full_spectrum(H_ss, wdf, Wf)
    extra = _injected_double(p1_f, Om_f)
    dbl_exact = exact_f[int(np.argmin(np.abs(exact_f - extra)))]
    err = abs(extra - dbl_exact)
    # also confirm the adiabatic backbone is one root short (the double is absent).
    assert exact_f.size == Om_f.size + 1, "adiabatic should miss exactly one double root"
    le("ii-B injected 0OS double err", err, 1e-6)
    return f"injected 0OS double P1={extra:.4f} vs exact={dbl_exact:.4f}, err={err:.2e}"


def t_dkb_p2_anchor_exact():
    # P2 (full B-matrix downfold) is the EXACT downfold onto the singles sector: in the
    # well-separated case (i.a) its single-like roots must match exact to < 1e-8.
    H_ss = np.array([[5.0, 0.25, 0.00],
                     [0.25, 5.8, 0.20],
                     [0.00, 0.20, 6.7]])
    omega_d = np.array([7.6])
    W = np.array([[0.5], [0.4], [0.3]])
    exact = dkb.exact_block_spectrum(H_ss, omega_d, W)
    p2 = dkb.p2_block_spectrum(H_ss, omega_d, W)
    p2_err, _ = dkb.adiabatic_match_error(p2, exact)   # Ns roots vs nearest exact
    le("i.a P2 anchor vs exact", p2_err, 1e-8)
    return f"P2 single-like max|err|={p2_err:.2e} (<1e-8)"


# ======================================================================
# 4. qmrsf_icpt2_multistate.py
# ======================================================================
def t_ms_heff_hermitian():
    # Multistate des Cloizeaux symmetric downfold: H_eff must be Hermitian (< 1e-10)
    # for BOTH denominator choices. Current validated residual ~0.
    case0 = ms.build_case(6, 6, thop=1.0, delta=0.0)
    worst = 0.0
    for denom in ("EN", "Dyall"):
        Ed, Heff, eP, hres = ms.icpt2_multistate(case0, nP=4, denom=denom)
        le(f"H_eff Hermiticity ({denom})", hres, 1e-10)
        # belt-and-suspenders: recompute residual directly from the returned matrix.
        direct = float(np.abs(Heff - Heff.T).max())
        le(f"H_eff sym ({denom}) direct", direct, 1e-10)
        worst = max(worst, hres, direct)
    return f"max Hermiticity residual (EN & Dyall) = {worst:.2e}"


def t_ms_avoided_crossing_smoothness():
    # Avoided-crossing continuity: scan the terminal donor/acceptor bias delta through
    # the avoided crossing of the two lowest states and compare the smoothness
    # (max|2nd diff|) of the multistate icPT2 curves to the smooth FCI reference. A kink
    # or root-flip would spike the metric. Current validated: FCI~3.2e-4, MS~1.2e-3
    # (~3.8x FCI). Guard: MS smooth (small absolute) AND comparable to FCI (< 8x).
    nP_ac = 2
    ls_ac = 0.4
    deltas = np.linspace(2.4, 4.4, 49)
    E_fci = np.zeros((len(deltas), 2))
    E_ms = np.zeros((len(deltas), 2))
    for i, dl in enumerate(deltas):
        c = ms.build_case(6, 6, thop=1.0, delta=dl)
        E_fci[i] = np.linalg.eigvalsh(c['Hfull'])[:2]
        Ed, _, eP, _ = ms.icpt2_multistate(c, nP=nP_ac, denom='EN', level_shift=ls_ac)
        E_ms[i] = Ed[:2]

    def maxd2(arr):
        return float(np.abs(np.diff(arr, 2, axis=0)).max())

    d2_fci = maxd2(E_fci)
    d2_ms = maxd2(E_ms)
    # absolute smoothness: no kink (current ~1.2e-3; guard generously at 1e-2).
    le("multistate max|2nd diff| (absolute)", d2_ms, 1e-2)
    # comparable to FCI: ratio must stay near the validated ~3.8x (guard < 8x).
    ratio = d2_ms / d2_fci
    le("multistate/FCI 2nd-diff ratio", ratio, 8.0)
    # no root flip: the two lowest states never cross (gap stays positive).
    gap_ms = E_ms[:, 1] - E_ms[:, 0]
    assert gap_ms.min() > 1e-6, f"multistate gap collapsed to {gap_ms.min():.2e} (root flip)"
    return (f"d2_FCI={d2_fci:.2e}, d2_MS={d2_ms:.2e} (ratio {ratio:.2f}x), "
            f"min gap={gap_ms.min():.3f}")


# ======================================================================
# Driver
# ======================================================================
def main():
    print("=" * 78)
    print("QMRSF prototype regression harness (pure NumPy)")
    print("=" * 78)

    print("\n-- qmrsf_icpt2_ppp_proto.py --")
    check("ppp: H Hermitian (max|H-H^T| < 1e-10)", t_ppp_hermitian)
    check("ppp: U=0 FCI == sum lowest orbital energies (< 1e-10)", t_ppp_u0_fci)
    check("ppp: EN downfold recovers CAS->FCI gap fraction in [0.55, 0.98]",
          t_ppp_recovery_band)

    print("\n-- qmrsf_dk_proto.py --")
    check("dk: adiabatic misses the double (error > 0.1)", t_dk_adiabatic_misses_double)
    check("dk: dressed pole search == exact (< 1e-8)", t_dk_dressed_matches_exact)

    print("\n-- qmrsf_dk_block_proto.py --")
    check("dkb: P1 recovers injected double, abstract well-separated (< 1e-6)",
          t_dkb_p1_recovers_double_abstract)
    check("dkb: P1 recovers injected 0OS double, real PPP ii-B (< 1e-6)",
          t_dkb_p1_recovers_double_ppp)
    check("dkb: P2 anchor == exact, well-separated (< 1e-8)", t_dkb_p2_anchor_exact)

    print("\n-- qmrsf_icpt2_multistate.py --")
    check("ms: H_eff Hermitian, EN & Dyall (< 1e-10)", t_ms_heff_hermitian)
    check("ms: avoided-crossing smoothness comparable to FCI", t_ms_avoided_crossing_smoothness)

    # ---- summary ----
    npass = sum(1 for _, ok, _ in _RESULTS if ok)
    ntot = len(_RESULTS)
    print("\n" + "=" * 78)
    print(f"SUMMARY: {npass}/{ntot} checks passed")
    if npass != ntot:
        print("FAILED checks:")
        for nm, ok, detail in _RESULTS:
            if not ok:
                print(f"  - {nm}: {detail}")
    print("=" * 78)
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
