"""Analytic nuclear derivative of the MRSF-TDDFT transition dipole.

This is the missing ingredient for analytic Herzberg-Teller vibronic
intensities: ``d mu_IJ / dR`` between two MRSF roots, with no finite
differences anywhere.

Theory
------
Within the MRSF response space the state-to-state transition dipole is the
interstate one-particle transition density contracted with the MO dipole
matrix,

    mu^k_IJ = - sum_pq  gamma^IJ_pq (b_I, b_J)  m^k_pq ,
    m^k     =  C^T R^k_AO C ,

with ``gamma^IJ`` exactly the object ``get_mrsf_transition_density`` builds
(occupied-occupied and virtual-virtual blocks plus the SOMO sqrt(2) cross
terms).  It is bilinear in the two packed amplitude vectors and traceless, so
``mu_IJ`` and its derivative are independent of the dipole origin.

Differentiating gives four contributions:

    T1  explicit AO dipole-integral derivative at frozen MOs and amplitudes
            -sum_uv gamma^IJ_AO,uv  dR^k_uv/dx
    T2  orbital response of the MO dipole matrix
            sum_pq Xdip_pq U^x_pq ,   Xdip = -m^k (gamma + gamma^T)
    T3  amplitude response of state I
            ytil_I^T (dA/dx) b_I
    T4  the same for state J

where ``A`` is the MRSF TDA response matrix and the eigenvector response is
eliminated by one projected resolvent solve per state and Cartesian dipole
component,

    (omega_I - A) ytil_I = Q_I w_I ,   w_I = d mu^k_IJ / d b_I ,
    Q_I = 1 - b_I b_I^T  (also projecting out the redundant OO coordinate).

For a nonadiabatic coupling the source ``w`` is the other eigenvector and the
solve collapses to ``b_I/(omega_J - omega_I)``; that closed form is what
``nac_analytic`` uses.  A transition property has no such collapse, so the
resolvent is solved explicitly here.  Everything else -- the frozen-MO
skeleton engines, the orbital source of the bilinear form, the ROHF/ROKS
adjoint Z solve and the derivative-integral contractions -- is the certified
analytic-NAC machinery, driven with a different property.

Each ``U^x`` dependence is collected into a single orbital source matrix
``Xtot = Xdip + Xamp(I) + Xamp(J)`` and handed to
``mrsf_nac_rohf_pair_overlap`` -> ``mrsf_nac_rohf_zvector`` -> HF/XC adjoints,
which fix the dependent (symmetric) rotations by ``U_sym = -S^x/2`` and set the
same-space antisymmetric rotations to zero.  That gauge is legitimate exactly
because the same-space antisymmetric part of ``Xtot`` vanishes identically --
a property that follows from the invariance of ``mu_IJ`` under block-diagonal
MO rotations and that :func:`analytic_transition_dipole_derivative` reports as
a self-check.

Validity
--------
* singlet fold only (``[tdhf] multiplicity = 1``);
* ROHF/ROKS reference (``scftype == 3``), no UMRSF;
* Cartesian AO basis -- the derivative-integral drivers inherited from
  ``der_overlap_matrix_ket`` perform no spherical reduction;
* linear response on a (near-)degenerate manifold is ill conditioned; the
  resolvent denominator is the diagnostic.
"""
import os

import numpy as np

RSQ2 = 1.0/np.sqrt(2.0)
SQ2M1 = np.sqrt(2.0) - 1.0

__all__ = ["MRSFTransitionDipole", "analytic_transition_dipole_derivative"]


def _fortran_view(raw, shape):
    """Re-read a TagArray buffer with its true Fortran shape."""
    return np.asarray(raw).ravel(order="C").reshape(shape, order="F")


def _unpack_lt(packed, n):
    m = np.zeros((n, n))
    r, c = np.tril_indices(n)
    m[r, c] = packed
    m[c, r] = packed
    return m


def _flat_fortran(matrix):
    return np.asfortranarray(matrix).reshape(-1, order="F").copy()


class MRSFTransitionDipole:
    """Fold bookkeeping and the transition-dipole bilinear algebra."""

    def __init__(self, mol):
        nbf = int(round(np.asarray(mol.data["OQP::VEC_MO_A"]).size ** 0.5))
        noca = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
        nocb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])
        self.mol = mol
        self.nbf, self.noca, self.nocb = nbf, noca, nocb
        self.nvirb = nbf - nocb
        self.nij = noca*self.nvirb
        self.mult = int(mol.config["tdhf"]["multiplicity"])
        self.C = np.asarray(mol.data["OQP::VEC_MO_A"]).reshape(nbf, nbf).T
        dip = _fortran_view(mol.data["OQP::dip_ao_org0"],
                            (nbf*(nbf + 1)//2, 3))
        self.R = [_unpack_lt(dip[:, k], nbf) for k in range(3)]
        self.m = [self.C.T @ r @ self.C for r in self.R]
        self.energies = np.asarray(mol.data["OQP::td_energies"]).ravel().copy()
        self.nstate = self.energies.size
        self.bvec = _fortran_view(mol.data["OQP::td_bvec_mo"],
                                  (self.nij, self.nstate))
        # packed slots of the folded open-open pair (SI Fig. S2 L and R)
        self.ilr1 = (noca - 1 - nocb - 1)*noca + noca - 2
        self.ilr2 = (noca - nocb - 1)*noca + noca - 1
        self.redundant = self.ilr2
        somo = (noca - 2, noca - 1)
        self.mask_vv = np.zeros((noca, self.nvirb))
        self.mask_oo = np.zeros((noca, self.nvirb))
        for k in somo:
            self.mask_vv[k, 0:2] = 1.0
            self.mask_oo[k, 0:2] = 1.0

    # -- packed <-> unfolded ------------------------------------------------
    def unfold(self, b):
        """Packed amplitude vector -> (noca, nvirb) rectangle."""
        noca = self.noca
        flat = np.array(b, dtype=float, copy=True)
        lr = flat[self.ilr1]
        flat[self.ilr1] = 0.0
        flat[self.ilr2] = 0.0
        x = flat.reshape((self.nvirb, noca)).T.copy()
        x[noca - 2, 0] = lr*RSQ2
        x[noca - 1, 1] = (-lr*RSQ2) if self.mult == 1 else (lr*RSQ2)
        return x

    def fold_adjoint(self, g):
        """Adjoint of :meth:`unfold`: rectangle gradient -> packed gradient."""
        noca = self.noca
        gg = np.array(g, dtype=float, copy=True)
        lr = gg[noca - 2, 0]*RSQ2
        lr += (-gg[noca - 1, 1]*RSQ2) if self.mult == 1 else (gg[noca - 1, 1]*RSQ2)
        gg[noca - 2, 0] = 0.0
        gg[noca - 1, 1] = 0.0
        out = gg.T.reshape(-1).copy()
        out[self.ilr1] = lr
        out[self.ilr2] = 0.0
        return out

    # -- transition density / dipole ---------------------------------------
    def gamma_mo(self, xi, xj):
        """Interstate 1-TDM in the alpha-MO basis (== get_mrsf_transition_density)."""
        nocb, noca = self.nocb, self.noca
        mvv, moo = self.mask_vv, self.mask_oo
        p = xj*mvv
        q = xi*mvv
        tvv = xj.T @ xi + SQ2M1*(p.T @ xi + xj.T @ q - 2.0*(p.T @ q))
        u = xi*moo
        v = xj*moo
        too = -(xi @ xj.T + SQ2M1*(u @ xj.T + xi @ v.T - 2.0*(u @ v.T)))
        trd = np.zeros((self.nbf, self.nbf))
        trd[nocb:, nocb:] += tvv
        trd[:noca, :noca] += too
        return trd

    def transition_dipole(self, xi, xj):
        trd = self.gamma_mo(xi, xj)
        return np.array([-float(np.sum(trd*mk)) for mk in self.m])

    def dipole_source(self, xi, xj, k):
        """d mu_k(xi, xj) / d xi on the unfolded rectangle (exact adjoint)."""
        m = self.m[k]
        mvv = m[self.nocb:, self.nocb:]
        moo = m[:self.noca, :self.noca]
        mask_v, mask_o = self.mask_vv, self.mask_oo
        p = xj*mask_v
        a1 = xj @ mvv
        a2 = p @ mvv
        wvv = -(a1 + SQ2M1*(a2 + a1*mask_v - 2.0*(a2*mask_v)))
        v = xj*mask_o
        b1 = moo @ xj
        b2 = moo @ v
        woo = b1 + SQ2M1*(b1*mask_o + b2 - 2.0*(b2*mask_o))
        return wvv + woo


# --------------------------------------------------------------------------
# resident engines
# --------------------------------------------------------------------------
def _dense_response_matrix(mol, ctx):
    """A of the MRSF TDA eigenproblem, one column per unit trial vector."""
    import oqp
    n = ctx.nij
    saved = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
    a = np.zeros((n, n))
    for col in range(n):
        work = np.array(saved, copy=True)
        flat = work.reshape(-1)
        flat[:n] = 0.0
        flat[col] = 1.0
        mol.data["OQP::td_bvec_mo"] = work
        oqp.mrsf_matvec_apply(mol)
        a[:, col] = np.array(mol.data["OQP::nac_mvax"], copy=True).ravel()[:n]
    mol.data["OQP::td_bvec_mo"] = saved
    red = ctx.redundant
    a[red, :] = 0.0
    a[:, red] = 0.0
    return 0.5*(a + a.T)


def _resolvent(a, omega, b, w, redundant):
    """Solve (omega - A) y = Q w with y orthogonal to b and to the redundant slot."""
    n = a.shape[0]
    keep = [i for i in range(n) if i != redundant]
    bb = b[keep]
    ww = w[keep] - bb*float(bb @ w[keep])/float(bb @ bb)
    shifted = omega*np.eye(len(keep)) - a[np.ix_(keep, keep)]
    proj = np.eye(len(keep)) - np.outer(bb, bb)/float(bb @ bb)
    regular = proj @ shifted @ proj + np.outer(bb, bb)/float(bb @ bb)
    yy = np.linalg.solve(regular, ww)
    yy -= bb*float(bb @ yy)/float(bb @ bb)
    y = np.zeros(n)
    y[keep] = yy
    residual = np.abs((omega*np.eye(n) - a) @ y - (w - b*float(b @ w)))
    residual[redundant] = 0.0
    return y, float(residual.max())


def _restore_bvec(mol, ctx, bvec):
    flat = np.zeros(ctx.nstate*ctx.nij)
    for s in range(ctx.nstate):
        flat[s*ctx.nij:(s + 1)*ctx.nij] = bvec[:, s]
    shape = np.array(mol.data["OQP::td_bvec_mo"]).shape
    mol.data["OQP::td_bvec_mo"] = flat.reshape(shape)


def _amplitude_term(mol, ctx, y, target, slot, bvec):
    """Skeleton vector and orbital source of  y^T (dA/dx) b_target."""
    import oqp
    natom = int(mol.data["natom"])
    mol.data["OQP::nac_ytil"] = np.ascontiguousarray(y, dtype=float)
    mol.data["OQP::nac_xstate"] = np.ascontiguousarray(bvec[:, target],
                                                       dtype=float)
    oqp.mrsf_nac_wpair(mol, slot + 1, target + 1)
    frozen = _fortran_view(np.array(mol.data["OQP::nac_mt_frozen"], copy=True),
                           (ctx.nbf, ctx.nbf))

    work = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
    flat = work.reshape(-1)
    flat[slot*ctx.nij:(slot + 1)*ctx.nij] = y
    mol.data["OQP::td_bvec_mo"] = work
    oqp.mrsf_nac_amp_pair(mol, slot + 1, target + 1)
    amp = np.array(mol.data["OQP::nac_amp"], copy=True).reshape(
        ctx.nstate, ctx.nstate, natom, 3)
    skeleton = amp[target, slot].reshape(-1).copy()
    oqp.mrsf_nac_esum(mol, slot + 1, target + 1)
    skeleton += np.array(mol.data["OQP::nac_esum"], copy=True).reshape(-1)
    _restore_bvec(mol, ctx, bvec)

    oqp.mrsf_nac_response(mol)
    response = _fortran_view(
        np.array(mol.data["OQP::nac_mt_response"], copy=True),
        (ctx.nbf, ctx.nbf))
    return skeleton, frozen + response


def _orbital_response(mol, ctx, xtot):
    """Dependent-MO term and Z-vector term of  sum_pq Xtot_pq U^x_pq."""
    import oqp
    mol.data["OQP::nac_mt_frozen"] = _flat_fortran(xtot)
    mol.data["OQP::nac_mt_response"] = np.zeros(ctx.nbf*ctx.nbf)
    mol.data["OQP::nac_gamma_pair"] = np.zeros(ctx.nbf*ctx.nbf)
    oqp.mrsf_nac_rohf_pair_overlap(mol)
    dependent = np.array(mol.data["OQP::nac_pair_vmask"], copy=True).reshape(-1)
    metric = np.array(mol.data["OQP::nac_pair_gsk"], copy=True).reshape(-1)
    if np.abs(metric).max() > 1.0e-12:
        raise RuntimeError("the NAC metric channel must be inactive here")
    oqp.mrsf_nac_rohf_zvector(mol)
    mol.data["OQP::nac_rohf_z"] = np.array(
        mol.data["OQP::nac_rohf_solution"], copy=True).reshape(-1)
    oqp.mrsf_nac_rohf_hf_adjoint(mol)
    zterm = np.array(mol.data["OQP::nac_rohf_hf_adjoint"], copy=True).reshape(-1)
    oqp.mrsf_nac_xc_adjoint(mol)
    zterm = zterm + np.array(mol.data["OQP::nac_rohf_xc_adjoint"],
                             copy=True).reshape(-1)
    return dependent, zterm


def _orbital_space(index, nocb, noca):
    if index < nocb:
        return 1
    return 2 if index < noca else 3


def analytic_transition_dipole_derivative(mol, pairs=None, response_matrix=None):
    """Analytic ``d mu_IJ / dR`` for the requested MRSF state pairs.

    Parameters
    ----------
    mol
        A molecule whose MRSF energy run has completed in this process.
    pairs
        Iterable of zero-based ``(I, J)`` state pairs.  Defaults to every
        unordered pair.
    response_matrix
        Pre-built dense ``A``; rebuilt when omitted.

    Returns
    -------
    dict
        ``(I, J) -> (3, 3*natom)`` derivative in a.u. per bohr, plus the key
        ``"diagnostics"`` carrying, per pair, the resolvent residuals and the
        same-space antisymmetric gauge residual of the orbital source.
    """
    import oqp

    if mol.config["tdhf"]["multiplicity"] != 1:
        raise NotImplementedError(
            "the analytic MRSF transition-dipole derivative implements the "
            "singlet fold only")
    if float(mol.config["scf"]["conv"]) > 1.0e-8:
        raise RuntimeError(
            "the analytic transition-dipole derivative requires [scf] conv "
            "<= 1e-8; 1e-10 is recommended")
    if float(mol.config["tdhf"]["conv"]) > 1.0e-8:
        raise RuntimeError(
            "the analytic transition-dipole derivative requires [tdhf] conv "
            "<= 1e-8; 1e-10 is recommended")

    oqp.electric_dipole_der_bra(mol)
    ctx = MRSFTransitionDipole(mol)
    natom = int(mol.data["natom"])
    ncoord = 3*natom
    nbf, noca, nocb = ctx.nbf, ctx.noca, ctx.nocb

    dbra = _fortran_view(np.array(mol.data["OQP::dip_dbra"], copy=True),
                         (nbf, nbf, 3, ncoord))
    dip_deriv = dbra + dbra.transpose(1, 0, 2, 3)

    control = mol.data._data.control
    cutoff = control.int2e_cutoff
    control.int2e_cutoff = 1.0e-20
    try:
        bvec = ctx.bvec.copy()
        amplitudes = [ctx.unfold(bvec[:, s]) for s in range(ctx.nstate)]
        amat = (_dense_response_matrix(mol, ctx) if response_matrix is None
                else response_matrix)

        same_space = np.zeros((nbf, nbf))
        for p in range(nbf):
            for q in range(nbf):
                if (_orbital_space(p, nocb, noca)
                        == _orbital_space(q, nocb, noca)):
                    same_space[p, q] = 1.0

        if pairs is None:
            pairs = [(i, j) for i in range(ctx.nstate)
                     for j in range(i + 1, ctx.nstate)]

        out = {}
        diagnostics = {}
        for (istate, jstate) in pairs:
            slot_i = jstate          # any slot other than the target
            slot_j = istate
            gamma = ctx.gamma_mo(amplitudes[istate], amplitudes[jstate])
            gamma_ao = ctx.C @ gamma @ ctx.C.T
            result = np.zeros((3, ncoord))
            info = []
            for k in range(3):
                explicit = np.array(
                    [-float(np.sum(gamma_ao*dip_deriv[:, :, k, c]))
                     for c in range(ncoord)])
                xdip = -ctx.m[k] @ (gamma + gamma.T)
                w_i = ctx.fold_adjoint(ctx.dipole_source(
                    amplitudes[istate], amplitudes[jstate], k))
                w_j = ctx.fold_adjoint(ctx.dipole_source(
                    amplitudes[jstate], amplitudes[istate], k))
                y_i, r_i = _resolvent(amat, ctx.energies[istate],
                                      bvec[:, istate], w_i, ctx.redundant)
                y_j, r_j = _resolvent(amat, ctx.energies[jstate],
                                      bvec[:, jstate], w_j, ctx.redundant)
                skel_i, x_i = _amplitude_term(mol, ctx, y_i, istate, slot_i,
                                              bvec)
                skel_j, x_j = _amplitude_term(mol, ctx, y_j, jstate, slot_j,
                                              bvec)
                xtot = xdip + x_i + x_j
                gauge = np.abs(0.5*(xtot - xtot.T)*same_space).max()
                dependent, zterm = _orbital_response(mol, ctx, xtot)
                result[k] = explicit + skel_i + skel_j + dependent + zterm
                info.append(dict(component=k, resolvent_residual=(r_i, r_j),
                                 gauge_residual=float(gauge),
                                 source_scale=float(np.abs(xtot).max())))
            out[(istate, jstate)] = result
            diagnostics[(istate, jstate)] = info
        out["diagnostics"] = diagnostics
        return out
    finally:
        control.int2e_cutoff = cutoff
        _restore_bvec(mol, ctx, ctx.bvec)
