# QMRSF-DK: the dressed (frequency-dependent) xc kernel for the 0OS double-spin-flip channel

Pathway II of `DESIGN_QMRSF_DUAL_PATHWAYS.md`. Companion prototype:
`qmrsf_dk_proto.py` (pure-NumPy mechanism demo). Created 2026-06-20.

This note derives the explicit frequency-dependent ("dressed") diagonal element that the
QMRSF-DK pathway adds to the **closed-shell 0OS double-spin-flip** configurations in the
QMRSF response matrix, shows why an adiabatic kernel cannot do it, and proposes a
prescription for the coupled CAS(4,4) block.

---

## 0. Setup and notation

QMRSF builds the Ms=0 response space of a quintet (S=2) ROHF/ROKS reference as a CAS(4,4)
of single-spin-flip (SSF) excitations from the four singly-occupied orbitals. Three config
classes (cf. design note):

- **2OS / 4OS** — open-shell configs reached by a *single* spin-flip (α→β in one of the
  four open shells). These are ordinary one-electron MRSF response operators: an adiabatic
  TDDFT kernel handles them correctly.
- **0OS** — six *closed-shell* configs reached by a **double** spin-flip (two α→β flips
  simultaneously). Relative to the closed-shell determinants that dominate the singlet
  manifold, a 0OS config differs by a **two-orbital (double) excitation**. This is the
  channel where the physics breaks.

Throughout, response amplitudes are labeled by configuration index; `ω` is the excitation
energy (eigenvalue of the response problem). We use the Casida/TDDFT response matrix in the
Tamm–Dancoff-friendly form where the "backbone" supplies an effective Hamiltonian-like
matrix `H` acting on configuration amplitudes, and `ω` are its eigenvalues. We write the
adiabatic backbone matrix as `A` (real, symmetric, frequency-independent) and add the
dressing as a frequency-dependent correction.

Key physical fact (Maitra, Cave, Zhang, Burke, *JCP* **120**, 5932 (2004); hereafter MCZB):

> An **adiabatic** (frequency-independent) xc kernel `f_xc` produces a response matrix whose
> dimension equals the number of *single* excitations. It therefore has **exactly** as many
> poles as there are single excitations — it **cannot** produce a state of predominantly
> double-excitation character. The missing double-excitation pole is restored **only** by a
> **frequency-dependent** kernel.

The 0OS class of QMRSF is precisely a set of *double* excitations relative to the singlet
closed-shell reference. So the adiabatic MRSF kernel that QMRSF inherits is structurally
blind to the dynamic correlation that the 0OS channel is supposed to carry — unless we
dress it. That is Pathway II.

---

## 1. The MCZB construction (single coupled to double)

Consider the canonical model: one single-excitation configuration `|s⟩` of energy `A` (the
adiabatic-TDDFT value) coupled, through the *bare* (frequency-independent, non-adiabatic)
two-electron interaction, to one double-excitation configuration `|d⟩` of energy `ω_d` with
coupling `V ≡ ⟨s|Ĥ|d⟩`. In the small *exact* (single + double) Hamiltonian basis,

```
        ┌            ┐
   H =  │   A     V  │                                                        (1)
        │   V*   ω_d │
        └            ┘
```

The exact spectrum is the two roots of `det(H − ω 1) = 0`:

```
   (A − ω)(ω_d − ω) − |V|^2 = 0.                                             (2)
```

This is a 2×2 problem with **two** roots. One root is "single-like", the other
"double-like". An adiabatic TDDFT, by contrast, only has the *single* sector `|s⟩` in its
basis — it never sees `|d⟩` and so produces at most **one** root near `A`. That is the
missing pole.

**Downfolding the double exactly (the dressing).** Eliminate the double amplitude `c_d` from
the eigenvalue equations `A c_s + V c_d = ω c_s`, `V* c_s + ω_d c_d = ω c_d`. The second
equation gives `c_d = V* c_s / (ω − ω_d)`. Substituting into the first yields a *scalar*,
*frequency-dependent* effective equation in the single sector alone:

```
   ω = A + B(ω),     B(ω) = |V|^2 / (ω − ω_d).                              (3)
```

Equation (3) is the **dressed-TDDFT eigenvalue condition**. The single-sector "diagonal" is
no longer the constant `A` but the *dressed diagonal*

```
   D(ω) = A + |V|^2 / (ω − ω_d).                                            (4)
```

`B(ω)` is exactly the frequency-dependent ("dynamical") part of the kernel — the piece an
adiabatic kernel sets to zero. Solving the **self-consistent** condition `ω = D(ω)`
(equivalently `[ω − A − B(ω)] = 0`) reproduces *both* roots of (2):
`(ω − A)(ω − ω_d) = |V|^2`, identical to (2). The single eigenvalue equation `ω = A + B(ω)`,
because `B` has a pole at `ω = ω_d`, supports **two** crossings of the line `y = ω` and
hence **two** solutions — the single-like and the double-like state. The pole the adiabatic
kernel lacked is injected by the `1/(ω − ω_d)` residue.

This is the MCZB/Casida result: the dynamical kernel matrix element in the dressed-TDDFT of
a single `s` coupled to a double `d` is

```
   f_xc^{dressed}(ω)_{ss}  ∝  |⟨s|Ĥ|d⟩|^2 / (ω − ω_d).                       (5)
```

(In the literature `A` is the adiabatic Casida matrix element `ω_s + 2⟨s|f_xc^{ad}|s⟩`-type
quantity; the dressing (5) is *added on top*. See also Loos & Blase, dynamical kernels,
*JCP* **150**, 124117 (2019); Romaniello *et al.*, *JCP* **130**, 044108 (2009) for the BSE
analog with the same `|V|^2/(ω − ω_d)` structure.)

---

## 2. The QMRSF-DK dressed 0OS diagonal

Specialize to QMRSF. Index the closed-shell 0OS double-spin-flip configuration as `0OS`,
with adiabatic backbone diagonal `A^{0OS}_adiabatic` (the value the existing MRSF backbone
already puts on that diagonal, including the adiabatic `v_xc`/`f_xc` contributions). The 0OS
config is a *double* excitation; the "satellite(s)" it must mix with are the genuine
doubly-excited (2SF) configurations `d` that share its space, with bare double energies `ω_d`
and couplings `V_{0OS,d} ≡ ⟨0OS|Ĥ|d⟩` through the (frequency-independent) electron–electron
interaction. Then the **dressed 0OS diagonal** is

```
   D(ω) = A^{0OS}_adiabatic  +  Σ_d  |V_{0OS,d}|^2 / (ω − ω_d).             (6)
```

and the dressed-response eigenvalue condition for that channel is the self-consistent

```
   [ ω − A^{0OS}_adiabatic − B(ω) ] = 0,    B(ω) = Σ_d |V_{0OS,d}|^2/(ω − ω_d).   (7)
```

Equation (6)/(7) is the central QMRSF-DK result. Properties:

1. **It injects the missing double-excitation pole.** Each satellite `d` contributes a
   pole in `B(ω)` at `ω = ω_d`; the self-consistent solution of (7) therefore has one extra
   root per satellite (the double-like states) beyond the single-like root near
   `A^{0OS}_adiabatic`. With `N_d` satellites the channel supports `1 + N_d` roots, exactly
   matching exact diagonalization of the `(1 + N_d)` model block (Section 4 / prototype).

2. **Correlation comes from the functional, not from an external Q space.** The satellite
   energies `ω_d` and couplings `V_{0OS,d}` are evaluated on the KS orbitals with the same
   two-electron interaction already in the backbone. No internally-contracted external-Q
   self-energy (that is Pathway I, icPT2) is introduced. The dressing reorganizes the
   *response* of the existing functional; it does not add a new correlation operator.

### 2.1 Avoiding double counting (the load-bearing subtlety)

The adiabatic part `A^{0OS}_adiabatic` — including the adiabatic `v_xc` and `f_xc` already
in the backbone diagonal — is **kept as is**. We add **only** the frequency-dependent
*residue* `B(ω)`, which is **identically zero** in the static (`ω`-independent) limit one
would obtain from an adiabatic kernel. Concretely:

```
   lim_{frequency-independent}  B(ω)  →  0 ,                                 (8)
```

i.e. an adiabatic kernel sets the dynamical residue to zero (MCZB). Therefore (6) does
**not** re-count the adiabatic contribution: the constant `A^{0OS}_adiabatic` is the full
adiabatic value, and `B(ω)` is a *pure addition* that vanishes in the adiabatic limit. There
is exactly **one** correlation source on that diagonal — the functional — split into

- the static piece `A^{0OS}_adiabatic` (already present in the backbone), and
- the dynamical piece `B(ω)` (the new dressing),

with no overlap between them. This is the formal statement of "no double counting" promised
in the design note (§4): unlike Pathway I, DK must **not** also add a Q-space self-energy on
the same diagonal, or the dynamical correlation would be counted twice.

A practical consistency requirement follows: the satellites `d` summed in (6) must be the
configurations that are **excluded** from the backbone response space (the genuine doubles
that the adiabatic single-particle response cannot reach), *not* configurations already
diagonalized inside the backbone. If a satellite is already in the backbone P-space, its
coupling must be removed from `B(ω)` to avoid re-folding a state that is treated exactly.
(In the prototype the "backbone" is just the bare single sector, so this reduces to: the
satellites are exactly the doubles left out of that sector.)

---

## 3. Solving the frequency-dependent eigenvalue problem

The condition `ω = A + B(ω)` (Eq. 3/7) is nonlinear in `ω`. Two equivalent solution routes:

**(a) Linearized / pole-expanded form (exact for the model).** Writing `B(ω) =
Σ_d |V_d|^2/(ω − ω_d)`, multiplying (7) through by `Π_d (ω − ω_d)` turns it into a
polynomial of degree `1 + N_d` whose roots are *exactly* the eigenvalues of the augmented
matrix

```
        ┌                              ┐
        │  A      V_1   V_2  ...  V_Nd │
        │  V_1*   ω_1    0   ...   0   │
   Ĥ =  │  V_2*    0    ω_2  ...   0   │ .                                   (9)
        │  ...                        │
        │  V_Nd*   0     0   ...  ω_Nd │
        └                              ┘
```

So **the dressed scalar equation (7) is the exact downfold of the augmented Hamiltonian (9)
onto its single sector**: solving (7) self-consistently and diagonalizing (9) give the same
spectrum. This is the cleanest validation target (prototype gate). It also shows the
dressing is *Hermitian/variational* — it is an exact Löwdin/Feshbach partitioning, not an
ad-hoc level shift.

**(b) Root search on `g(ω) = ω − A − B(ω)`.** Because `B(ω)` has simple poles at each
`ω_d`, `g(ω)` is monotone increasing on each interval between consecutive poles and crosses
zero exactly once per interval (plus once above the top pole). A bracketed bisection/secant
per interval finds all `1 + N_d` roots robustly without forming (9). This is the route that
generalizes when the satellites are too numerous to augment explicitly, or when `V_d`, `ω_d`
are only available as on-the-fly contractions.

---

## 4. Well-definedness beyond one single + one double: the CAS(4,4) block

MCZB derive the dressing **rigorously only for the isolated one-single/one-double mixing
limit** (Eqs. 1–5): a single `|s⟩` resonant with a single `|d⟩`, well separated from the
rest of the spectrum. QMRSF's response space is a **coupled 20-singlet CAS(4,4) block**
(plus the triplet/quintet blocks), in which several singles and several 0OS doubles mix
simultaneously. The clean scalar downfold (3) is no longer uniquely defined, because the
"single sector" is now multidimensional and the partition single↔double is basis-dependent.
Two admissible prescriptions:

### Prescription P1 — state-by-state pole search (recommended default)

Diagonalize the **adiabatic** backbone block first to obtain the singlet states and their
energies `{Ω_k^{ad}}` and vectors `{C^{(k)}}`. For each target state `k`, project the
dynamical dressing onto that state: define the *state-resolved* coupling to satellite `d`

```
   Ṽ_{k,d} = Σ_{0OS configs c}  C^{(k)}_c  V_{c,d} ,                        (10)
```

and solve the scalar dressed equation **per state**

```
   ω = Ω_k^{ad} + Σ_d |Ṽ_{k,d}|^2 / (ω − ω_d).                             (11)
```

This recovers the MCZB scalar form (3) with `A → Ω_k^{ad}`, `V → Ṽ_{k,d}`, and the root
search of §3(b) returns the dressed state energy (and, near a resonance, the extra
double-like root that becomes a *new* state). It is exactly the limit in which MCZB is
rigorous, applied one state at a time.

- **Assumptions / limits.** (i) The adiabatic states are a good zeroth order — the dressing
  is a perturbative reorganization, valid when `|Ṽ_{k,d}|` is small compared to the spacing
  to neighboring adiabatic states. (ii) Off-diagonal *re-coupling* between two adiabatic
  states *through* a shared satellite is neglected (each state is dressed independently). Near
  a true single/single avoided crossing this can mis-assign the extra root; a degeneracy
  guard is needed.

### Prescription P2 — block dressed kernel (frequency-dependent matrix)

Keep the dressing as a **matrix** in the 0OS sub-block and solve the frequency-dependent
matrix eigenproblem

```
   [ A_block + B_block(ω) ] C(ω) = ω C(ω),                                 (12)

   [B_block(ω)]_{c c'} = Σ_d  V_{c,d} V_{c',d}* / (ω − ω_d).               (13)
```

This is again the exact downfold of the augmented matrix
`[[A_block, W],[W†, diag(ω_d)]]` (the multi-channel generalization of (9), `W_{c,d}=V_{c,d}`)
onto the 0OS sub-block, so it is well-defined and Hermitian/variational. It preserves
off-diagonal re-coupling that P1 drops, and it cannot split degenerate multiplets because it
is a similarity-faithful downfold of a symmetric augmented matrix (symmetry of `A_block` and
`diag(ω_d)` is inherited).

- **Assumptions / limits.** (i) Requires the full satellite list and couplings as a matrix
  (more expensive than P1's per-state contraction). (ii) The frequency-dependent matrix
  eigenproblem is solved by iterating `ω` to self-consistency (or by one-shot diagonalization
  of the augmented matrix when `N_d` is small); intruder satellites with `ω_d` inside the
  target window need a regularizer, as in any partitioned PT.

**Recommendation.** Default to **P1** (cheap, matches the rigorous MCZB limit, easy to gate
for degeneracy preservation) and use **P2** as the reference/validation when off-diagonal
re-coupling matters (near-degenerate 0OS configs, e.g. the CBD square point). Both reduce to
the exact scalar/matrix downfold of an augmented Hamiltonian, which is the property the
prototype validates against exact diagonalization.

---

## 5. Summary of the key equations

| Eq.  | Statement |
|------|-----------|
| (1)  | exact (single+double) 2×2 model Hamiltonian |
| (2)  | exact two-root secular equation `(A−ω)(ω_d−ω) − |V|^2 = 0` |
| (3)  | **dressed scalar eigenvalue condition** `ω = A + |V|^2/(ω − ω_d)` |
| (4)  | dressed diagonal `D(ω) = A + |V|^2/(ω − ω_d)` |
| (5)  | dynamical kernel element `∝ |⟨s|Ĥ|d⟩|^2/(ω − ω_d)` |
| (6)  | **QMRSF-DK dressed 0OS diagonal** `D(ω) = A^{0OS}_ad + Σ_d |V_{0OS,d}|^2/(ω − ω_d)` |
| (7)  | self-consistent condition `[ω − A^{0OS}_ad − B(ω)] = 0` |
| (8)  | adiabatic limit `B(ω) → 0` ⇒ **no double counting** |
| (9)  | augmented Hamiltonian whose exact spectrum = roots of (7) |
| (10)–(11) | **P1**: state-resolved coupling + per-state scalar dressing |
| (12)–(13) | **P2**: block frequency-dependent matrix dressing |

The prototype `qmrsf_dk_proto.py` validates (1)–(9) numerically: it shows the **adiabatic**
treatment (drop `B`, project out the double) **misses the doubly-excited state** and gets the
single-state energy wrong, while solving the **dressed** equation (3)/(7) recovers the exact
two-root (and multi-root) spectrum to machine precision, validated against exact
diagonalization of (1)/(9). It also demonstrates the P1 pole-search on a small "0OS-like"
3-state block.
