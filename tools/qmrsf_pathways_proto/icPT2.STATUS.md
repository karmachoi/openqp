# QMRSF-icPT2 — Status & Handoff

_Last updated: 2026-06-21. Single source of truth for the QMRSF-icPT2 manuscript,
benchmarks, reference calculations, code, and open issues._

---

## 0. TL;DR

QMRSF-icPT2 = a single high-spin **quintet ROHF** reference → spin-pure **CAS(4,4)**
double-spin-flip backbone (M_s=0, 36 dets = 20 singlet + 15 triplet + 1 quintet) →
**internally-contracted external-Q second-order downfold** for dynamic correlation
(EN and Dyall partitions; Dyall is the production default), Hermitised by a multistate
des Cloizeaux average. One SCF + one-shot O(N^5).

The proper JCTC-style paper is **`main.v1.tex`** in the `overleaf-icpt2` project.
It is **pushed and current** on Overleaf `main` (commit `764404b`).

**Headline results:** H4 reproduces in-window FCI; CBD automerization barrier at the
reference geometry = **icPT2-Dyall +7.37 kcal/mol** (in the literature
CASPT2(4,4)–NEVPT2(4,4) bracket, near TBE 8.93); spin-resolved T1 excitations match the
TBE; scalability streams 4.06M perturbers (CBD/cc-pVDZ).

**Known limits (honest):** the active-only downfold under-correlates **doubly-excited
dark states** (butadiene 2¹Ag, CBD ionic singlets); the contracted engine **fails at
cc-pVTZ for CBD** (~330M perturbers). See §6.

---

## 1. Repositories, build & run

### Manuscript
- **Overleaf project `overleaf-icpt2`**: `git.overleaf.com/6a3618d3b7f6eca6d4c9e24b`
  - local clone: `~/Documents/Manuscripts/overleaf-icpt2`
  - **`main.v1.tex`** = the paper (compile with `latexmk -pdf main.v1.tex`, ≥2 passes).
  - `refs.bib` = bibliography.
  - `Figures/` = the figure PDFs (generated from the session folder, see §7).
  - Policy: **branch first, merge to `main`, then push** (a co-author, Alireza
    Lashkaripour, shares the project). Never edit `main` directly / rebase.
  - Working branch used this session: `claude/icpt2-caspt2-reference` (merged to main).

### Engine (OpenQP)
- Worktree: **`~/Documents/openqp-private-qmrsf-pathways`**, branch **`feat/qmrsf-dual-pathways`**.
  Implements `type=qmrsf_icpt2` (and `type=qmrsf_dk`).
- Build:
  ```
  cmake -S . -B build -G Ninja -DCMAKE_{C,CXX,Fortran}_COMPILER=g{cc,++,fortran}-15 \
    -DUSE_LIBINT=OFF -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF \
    -DENABLE_OPENTRAH=OFF -DENABLE_PYTHON=OFF -DBUILD_TESTING=OFF
  ninja -C build oqp
  ```
  Stage a private root at `/tmp/qmrsf_root` (include/oqp.h → worktree, lib/liboqp.dylib →
  build/source, share → installed share).
- Run:
  ```
  OPENQP_ROOT=/tmp/qmrsf_root PYTHONPATH=<wt>/pyoqp \
    python3 <wt>/pyoqp/oqp/pyoqp.py input.inp
  ```
- Input keys: `[scf] multiplicity=5 type=rohf converger_type=diis`,
  `[tdhf] type=qmrsf_icpt2 qmrsf_icpt2_shift=0.1`. Use `save_molden=False save_mol=False`.
  Basis names route to BSE (`cc-pvdz`, `6-31g`, `cc-pvtz` OK; `6-31g(d)` FAILS → use `6-31g`).

### Reference calculations (external tools — allowed as references, NOT inside OpenQP)
- **ORCA 6.1**: `/Users/cheolhochoi/Library/orca_6_1_0/orca`.
  - `PTMethod FIC_CASPT2` (canonical CASPT2), `FIC_NEVPT2` (= PC-NEVPT2), `FIC_CASPT2K`.
  - CASPT2 uses **no IPEA / no level shift** here (standard H0).
  - **Gotcha:** FIC_CASPT2 converges the PT2 energy then **aborts (exit 11) in a
    post-energy step**. Parse the energy anyway: `CASPT2 total = (Final CASSCF energy)
    + (converged EPT2 from the iteration table)`. NEVPT2 exits cleanly (FINAL SINGLE
    POINT ENERGY). Escape `%` as `%%` in Python format strings for ORCA inputs.
- **pyscf**: allowed as an external reference tool only (NOT in the OpenQP engine, build,
  or deps). Used for FCI and NEVPT2 cross-checks.

---

## 2. Manuscript structure (current, commit 764404b, 16 pp)

§1 Introduction · §2 Theory · §2.x Test systems & basis sets ·
**§3 Results** = three tiers: **Correctness** (in-space exactness/spin purity; dynamic
correlation recovery), **Physics** (geometry scans; CBD barrier; spin-resolved
excitations; butadiene internal conversion), **Scalability** · §4 Conclusions.

**Tables (current numbering):**
1. `tab:ground` — ground-state totals (H4 STO-3G/6-31G/cc-pVDZ + CBD rect/square),
   columns ROHF/CAS/EN/Dyall/CASPT2/PC-NEVPT2/FCI. CBD CASPT2 shown **active-only**
   (4e, matched to icPT2) with **all-valence** values in the caption.
2. `tab:barrier` — CBD automerization barrier at the **reference (Loos) geometry**,
   cc-pVDZ (this work) vs aug-cc-pVTZ literature reference.
3. `tab:spin` — spin-resolved excitations (T1/S1/S2), CBD D2h+D4h + H4, cc-pVDZ,
   CAS/EN/Dyall/CASPT2(4,4)/PC-NEVPT2(4,4)/Ref(TBE or FCI).
4. `tab:scal` — scalability (perturber count + RAM).

**Figures:** 1 mixed-reference schematic (tikz) · 2 P/Q partition (tikz) · 3 test-systems
molecules · 4 H4 scans (with CASPT2/NEVPT2) · 5 butadiene PES · 6 scalability.

**Removed this session:** validation-gate table, correlation-recovery figure,
MRSF-comparison table+subsection, spin-manifold figure, TMM/TME polyradical table.
Dev-flavoured wording (gate/oracle/pyscf-free) scrubbed.

---

## 3. Key validated numbers

**Ground states (E_h):**
- H4/STO-3G: CAS = EN = Dyall = CASPT2 = NEVPT2 = FCI = **−2.102608** (no external space).
- H4/6-31G: FCI −2.181861 | CASPT2 −2.180552 | NEVPT2 −2.180219 | EN −2.187565 |
  Dyall −2.200845. CASPT2/NEVPT2 hug FCI from above (1–2 mEh); icPT2 overshoots below.
- H4/cc-pVDZ: FCI −2.207431 | CASPT2 −2.202648 | NEVPT2 −2.199945 | Dyall −2.247510.
- CBD/cc-pVDZ rect: ROHF −153.3963 | CAS −153.7013 | EN −153.7340 | Dyall −153.7415.
  (6-31G rect Dyall −153.6571.)
  - CASPT2 active-only (4e) −153.6949; **all-valence** −154.1852 (correlates all 20
    valence e⁻ → not directly comparable to active-only icPT2 totals; compare differences).

**CBD automerization barrier (kcal/mol), at the Loos D2h/D4h reference geometries:**
- icPT2 cc-pVDZ: CAS **+2.05**, EN **+6.40**, Dyall **+7.37**.
- Literature (Monino/Loos 2022, aug-cc-pVTZ): CASPT2(4,4) 7.77, PC-NEVPT2(4,4) 9.24,
  **TBE 8.93**. icPT2-Dyall lands inside the bracket — magnitude, not just sign.
- (Idealized, unrelaxed geometry gave a wrong-sign CAS barrier −0.96; geometry-sensitive.)

**Spin-resolved excitations (eV), Loos geom, cc-pVDZ (T1/S1/S2):**
- CBD D2h Dyall 1.43/3.99/4.35; NEVPT2 1.39/3.03/4.00; **TBE 1.46/3.13/4.04**.
- CBD D4h Dyall 0.19/1.93/2.72; NEVPT2 0.08/1.26/1.54; **TBE 0.14/1.50/1.85**.
- H4 Dyall 3.54/7.95/8.06; **FCI 3.30/7.49/8.03**.
- T1 matches the TBE well everywhere; ionic/doubly-excited singlets are overestimated by
  ALL (4,4) methods (icPT2 and CASPT2(4,4) alike) — frontier-only-active-space signature.

**Butadiene FC→1¹Bu→2¹Ag internal conversion (cc-pVDZ, Fig 5):**
- bright 1¹Bu (single excitation): icPT2 within ~0.1 eV of SA-CASSCF(4,4)+PC-NEVPT2.
- dark 2¹Ag (double excitation): CAS(4,4) ~8.9 eV → icPT2 ~7.6 → NEVPT2 ~6.2→5.0
  (funnel below 1¹Bu). icPT2 improves on CAS but leaves 2¹Ag ~1.5–2 eV high.

**Scalability:** CBD/cc-pVDZ = **4,064,220 streamed perturbers**, 21 s, 4.2 GB peak RAM.
Validated closed-form perturber count `count_model(nvirt)`.

---

## 4. Cross-validation (all consistent)
- FCI: project NumPy det-CI == pyscf (exact).
- CASSCF: ORCA == pyscf (exact).
- NEVPT2: ORCA ≈ pyscf (≈0.07 mEh = SC- vs FIC-NEVPT2 variant difference).
- icPT2 vs independent NumPy/closed-form oracle: 1e-8…1e-15 (the old validation gates).

---

## 5. The single most important interpretive point

icPT2's external-Q downfold correlates **only the four active electrons**. This is why:
- it reproduces **single excitations** (bright 1¹Bu, all triplets) at PC-NEVPT2 quality;
- it **under-correlates doubly-excited dark states** (butadiene 2¹Ag, CBD ionic singlets)
  that draw their stabilization from **inactive-valence** dynamic correlation;
- CBD **total** energies are not directly comparable to standard (all-valence) CASPT2;
  only **differences** (barrier, excitations) are.

The clear improvement direction: **extend the downfold to inactive-shell correlation.**

---

## 6. Open issues / decisions pending

1. **[PENDING USER DECISION] cc-pVTZ barrier.** icPT2 cc-pVTZ for CBD **fails** — the
   contracted downfold returns the bare CAS (EN=Dyall=CAS) at ~330M perturbers; engine
   limit beyond the 4M cc-pVDZ case. Options offered: (a) keep icPT2/cc-pVDZ vs the
   aug-cc-pVTZ literature reference [recommended]; (b) engine-level work to push icPT2 to
   cc-pVTZ; (c) redo our CASPT2/NEVPT2 with **state-averaged** CASSCF for a reliable
   in-house cc-pVTZ barrier column. **Currently option (a) is in the paper.**
2. **Our ORCA CASPT2/NEVPT2 barrier is unreliable.** State-specific CASSCF(4,4)
   over-stabilizes the D4h diradical (CASSCF barrier −25 kcal/mol vs Loos +7.38), so the
   paper uses the **published** Loos CASPT2/NEVPT2, not ours. SA-CASSCF would fix this
   (= option 1c). Nice corollary worth a sentence: the quintet reference is a more
   *balanced* backbone than state-specific CASSCF for this barrier.
3. **Active-only correlation limit** (§5) — the headline methodological caveat; motivates
   the inactive-correlation extension.
4. **refs.bib author annotations** flagged `CHECK:` in the original paper still need a
   pass before submission (park2025quintet, casanova2015rasnsfpt2, etc.).
5. Quantitative benchmark on optimized geometries + polarized basis vs XMS-CASPT2(4,4) /
   MR-ADC(2) remains the target to beat.

---

## 7. Code & data inventory

### Engine harness — `~/Documents/openqp-private-qmrsf-pathways/tools/qmrsf_pathways_proto/`
- `qmrsf_icpt2_ppp_proto.py` — pure-NumPy proto (det-CI, spinorb, icPT2 downfold, FCI).
- `qmrsf_icpt2_contracted_proto.py` — contracted external-Q engine + `count_model` (perturber count).
- `stageB/run_benchmarks.py` — main driver (H4/CBD, STO-3G/6-31G; emits tables).
- `stageB/scan_h4.py` — H4 square→rect + linear geometry scans (Fig 4 data).
- `stageB/run_polyrad.py` — TMM/TME (programmatic geometries). [table removed from paper]
- `stageB/run_scalability.py` — perturber count + wall/RAM (Fig 6 / Table 4).
- `stageB/run_tablefill.py` — H4/cc-pVDZ + CBD rect/square cc-pVDZ icPT2 (Table 1).
- `stageB/run_loos_barrier.py` — CBD icPT2 barrier at the **reconstructed Loos geometries**,
  cc-pVDZ + cc-pVTZ (cc-pVTZ fails, see §6).
- `stageB/run_icpt2_tests.sh` + `gate_*.py`, `route_a_oracle.py` — the validation suite (7/7).

### Reference-calc scripts — session folder `caspt2_probe/` and `butadiene/`
Session folder: **`/Volumes/External_Storage/claude/sessions/20260620_223000_icpt2_bench_figs/`**
- `caspt2_probe/run_caspt2_ref.py` — ORCA CASPT2/NEVPT2 for H4 scans + ground.
- `caspt2_probe/cbd_caspt2.py` — standard (all-valence) CBD CASPT2 barrier.
- `caspt2_probe/cbd_activeonly.py` — active-only (frozen-12) CBD CASPT2 (Table 1).
- `caspt2_probe/loos_barrier_caspt2.py` — our CASPT2/NEVPT2 barrier at Loos geom (unreliable, §6).
- `caspt2_probe/run_spin_caspt2.py` — SA-CASSCF+NEVPT2/CASPT2 excitations (Table 3).
- `butadiene/scan_butadiene.py` — QMRSF-icPT2 PES along the FC→1Bu→2Ag path.
- `butadiene/run_butadiene_nevpt2.py` — SA-CASSCF+PC-NEVPT2 overlay.
- `butadiene/make_butadiene_fig.py` — Fig 5 (energies rel. to FC S0; CAS/QMRSF/NEVPT2 tiers).
- `scratch/make_molecules.py`, `make_realdata_figs.py`, `make_scan_fig_v2.py`,
  `make_scalability_and_tables.py` — figure generators.

### Data (session folder, JSON)
`caspt2_results.json` (H4 scans+ground CASPT2/NEVPT2), `loos_barrier.json` (icPT2 barrier,
cc-pVTZ broken), `loos_barrier_caspt2.json` (our CASPT2 barrier — unreliable),
`spin_caspt2.json` (SA NEVPT2/CASPT2 excitations), `spin_icpt2_loos_D{2,4}h_ccpvdz.json`
(icPT2 excitations at Loos geom), `cbd_activeonly_caspt2.json`, `cbd_caspt2_barrier.json`,
`scan_results.json`, `scalability_results.json`, `table_fill.json`,
`butadiene_scan_ccpvdz.json`, `butadiene_nevpt2_ccpvdz.json`, `polyrad_results.json`.

### Figures: session `figs/` → copied into `overleaf-icpt2/Figures/`
In the paper: `test_systems.pdf`, `fig_h4_scan.pdf`, `fig_butadiene.pdf`, `fig_scalability.pdf`.
Not used: `fig_corr_recovery.pdf`, `fig_spin_manifold.pdf`, `mol_*.pdf`.

### Reference geometries (reconstructed from Monino/Loos 2022, Table 1)
CBD ring is a rectangle/square (internal angle 90° ⇒ H–C–C = 135°, H along the external
bisector at C–H = 1.076 Å):
- D2h (1¹Ag, min): C–C 1.354 / 1.566 Å.  → C(±0.677,±0.783), H(±1.438,±1.544).
- D4h (1¹B1g, TS): C–C 1.449 Å.          → C(±0.7245,±0.7245), H(±1.485,±1.485).
Butadiene FC→1Bu→2Ag geodesic path: SI of jz1c02707 (`jz1c02707_si_002.pdf`,
parsed by `scan_butadiene.py`).

---

## 8. References
- **Monino, Boggio-Pasqua, Scemama, Jacquemin, Loos**, "Reference Energies for
  Cyclobutadiene," *J. Phys. Chem. A* 2022, 126, 4664. DOI 10.1021/acs.jpca.2c02480.
  (CBD barrier TBE 8.93; D2h/D4h geometries; CASPT2/NEVPT2(4,4) and (12,12).)
- **Park, Komarov, Filatov, Choi**, "Internal Conversion between Bright (1¹Bu) and Dark
  (2¹Ag) States…," *J. Phys. Chem. Lett.* 2021, 12. DOI 10.1021/acs.jpclett.1c02707.
  (Butadiene path geometries.)
- ORCA: Neese et al., *JCP* 2020, 152, 224108. pyscf: Sun et al., *JCP* 2020, 153, 024109.

---

## 9. Suggested next steps
1. Resolve the cc-pVTZ barrier question (§6.1) — likely option (a) or (c).
2. Inactive-correlation extension of the downfold (the real method advance; would fix the
   butadiene 2¹Ag and CBD ionic-singlet under-correlation).
3. Engine: lift the contracted-downfold perturber wall to reach cc-pVTZ for CBD-sized π systems.
4. Optimized-geometry / polarized-basis benchmark vs XMS-CASPT2(4,4) and MR-ADC(2).
5. Pre-submission: refs.bib `CHECK:` pass; final honesty/overclaim read of §3.
