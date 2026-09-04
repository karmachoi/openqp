# MECI search on the exact branching plane (analytic MRSF derivative coupling)

`meci_exact_plane.py` minimizes the energy of a degenerate pair on its seam with the
Bearpark–Robb composite gradient

    F = 2 dE g_hat + (1 - g_hat g_hat^T - h_hat h_hat^T) s ,

where g is half the gradient difference, h the analytic MRSF derivative-coupling vector
(`oqp.library.nac_analytic.analytic_nac`, energy-weighted), s the mean gradient and dE the gap.
The plane is exact at every step; no update-branching-plane approximation and no penalty
function are used. Steps are a trust-limited BFGS on F (initial inverse Hessian 2 bohr^2/Eh);
the trust radius grows by 1.2 when the mean energy falls with the gap below 2e-4 Eh and is
reduced (and the BFGS memory reset) otherwise.

Usage (inside a venv built from this branch):

    python meci_exact_plane.py input.inp MAXIT TRUST_BOHR [start_geometry_bohr.txt]

`input.inp` is an ordinary MECI input with `[optimize] lib=scipy meci_search=ubp istate=I jstate=J`
(1-based response roots; for MRSF singlets S0 is root 1) and `[tdhf] conv <= 1e-8` (1e-10
recommended, 1e-14 stalls the Davidson). The driver only borrows the energy/gradient machinery of
`MECIOpt`; the optimization loop is its own. Each step appends to `sd_trace.txt` the pair energies,
gap, |g|, |h|, the projected mean gradient, the composite gradient, and the Yarkony invariants
(delta_gh, Delta_gh, P, B) computed after the rotation that makes g and h orthogonal, and saves
`point_N.npz` (geometry, g, h, mean gradient, state gradients, energies). Convergence is declared
at max|F| < 3e-5 Eh/bohr and |dE| < 5e-6 Eh; `converged_geom_bohr.txt` is written then.

Results obtained with it on 2026-09-04 (MRSF-BH&HLYP/6-31G*):
- trans-butadiene S2/S1: the penalty-search point already has max|F| = 3.7e-5 (a seam minimum).
- thymine S2/S1: the penalty-search point has max|F| = 3.7e-2; the search lowers the pair by
  0.57 eV to E(S1) = -453.6763550 Eh with max|F| = 3.5e-4, rms 8.7e-5, gap 0.011 meV,
  P = 1.04 (boundary), B = 1.05, Delta_gh = +0.987.
