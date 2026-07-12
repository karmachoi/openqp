"""Shared dispatch helpers for the tight-binding (TB) backends.

OpenQP currently ships two decoupled TB backends behind the same PyOQP
plumbing:

* ``method=dftb`` -> the [dftb] section and the openqp-dftb library
  (``oqp.library.openqp_dftb.OpenQPDFTBAdapter``);
* ``method=xtb``  -> the [xtb] section and the openqp-xtb library
  (``oqp.library.openqp_xtb.OpenQPXTBAdapter``).

Every ``method == 'dftb'`` dispatch site in single_point.py / namd.py /
qmmm_driver.py / runfunc.py goes through the helpers below so both backends
share one code path.

INTENTIONALLY SHARED between the two backends (do not fork per backend):

* the wavefunction/response tag names the adapters publish and the drivers
  consume -- ``OQP::dftb_wf_dims``, ``OQP::VEC_MO_A``, ``OQP::VEC_MO_B``,
  ``OQP::E_MO_A``, ``OQP::E_MO_B``, ``OQP::td_bvec_mo``, ``OQP::td_energies``,
  ``OQP::partial_charges``, the ``OQP::soc_*`` family, ...;
* the ``mol.dftb_external_potential`` attribute that carries the per-atom
  QM/MM embedding potential into the SCC Hamiltonian.

Both libraries implement the same C API family, so the drivers stay agnostic:
only the section name, the shared-library/symbol prefix, and a few model
options differ (see the adapter classes).
"""

from __future__ import annotations


TB_METHODS = {"dftb", "xtb"}


def is_tb_method(method_or_mol) -> bool:
    """True when the input method is one of the tight-binding backends.

    Accepts either the method name itself (any case) or an object carrying a
    ``config`` dict (a Molecule), so both ``is_tb_method(self.method)`` and
    ``is_tb_method(mol)`` read naturally at the dispatch sites.
    """
    if isinstance(method_or_mol, str):
        return method_or_mol.strip().lower() in TB_METHODS
    config = getattr(method_or_mol, "config", None)
    if config is None:
        return False
    return str(config.get("input", {}).get("method", "")).strip().lower() in TB_METHODS


def tb_section_name(config) -> str:
    """Name of the TB options section for this config: 'dftb' or 'xtb'.

    The section name equals the input method by construction ([dftb] for
    method=dftb, [xtb] for method=xtb).
    """
    method = str(config.get("input", {}).get("method", "")).strip().lower()
    if method not in TB_METHODS:
        raise ValueError(
            f"tb_section_name called for non-TB method {method!r}; "
            f"expected one of {sorted(TB_METHODS)}"
        )
    return method


def tb_config(config) -> dict:
    """The TB options section ([dftb] or [xtb]) of this config, {} by default."""
    return config.get(tb_section_name(config), {})


def make_tb_adapter(mol):
    """Instantiate the adapter matching mol's input method (dftb or xtb)."""
    section = tb_section_name(mol.config)
    if section == "xtb":
        from oqp.library.openqp_xtb import OpenQPXTBAdapter  # noqa: PLC0415
        return OpenQPXTBAdapter(mol)
    from oqp.library.openqp_dftb import OpenQPDFTBAdapter  # noqa: PLC0415
    return OpenQPDFTBAdapter(mol)


# --------------------------------------------------------------------------
# tb_operator: user-facing operator presets for the MRSF response exchange.
# --------------------------------------------------------------------------
#
# A single [dftb]/[xtb] key, ``tb_operator = <name>``, selects a named bundle
# of the individual response-exchange flags below so a user does not have to
# hand-set the ~10 knobs that define a given MRSF-TB response operator. The
# preset dict is backend-agnostic: both the openqp-dftb and openqp-xtb C ABIs
# now carry the same three DTCAM-TB response knobs (c_mrsf,
# response_global_hybrid, onsite_exchange_scale) in the same position, so one
# dict serves both backends.
#
# The four presets:
#   * ``yukawa``    -- default LR-MRSF-DFTB (Yukawa-Slater gamma^lr), the
#                      historical baseline (n->pi* singlet-triplet split == 0);
#   * ``erf``       -- same but with the erf(omega R)/R long-range gamma;
#   * ``erf-tuned`` -- the retuned erf kernel (omega=0.25, cam_beta=1.2, LC
#                      ground state);
#   * ``dtcam-tb``  -- the doubly-tuned DTCAM-TB operator that generates the
#                      paper numbers (erf long-range gamma, response global
#                      hybrid, Garcia beyond-monopole on-site exchange scaled
#                      1.60, and c_MRSF=1.02). NOTE the spc_*=-1.0 entries are
#                      REQUIRED: a negative per-channel spin-pairing coupling
#                      routes the Fortran merge(spc, exchange_fraction, spc>=0)
#                      onto response_exchange_fraction == c_mrsf (1.02), i.e.
#                      it ties the spin-pairing scale to c_MRSF instead of the
#                      0.5 baseline used by the other presets.
TB_OPERATOR_PRESETS = {
    'yukawa': {
        'lc_gamma': 'yukawa', 'omega': 0.3, 'cam_alpha': 0.0, 'cam_beta': 1.0,
        'spc_coco': 0.5, 'spc_ovov': 0.5, 'spc_coov': 0.5,
        'c_mrsf': -1.0, 'response_global_hybrid': 0, 'onsite_exchange_scale': 0.0,
    },
    'erf': {
        'lc_gamma': 'erf', 'omega': 0.3, 'cam_alpha': 0.0, 'cam_beta': 1.0,
        'spc_coco': 0.5, 'spc_ovov': 0.5, 'spc_coov': 0.5,
        'c_mrsf': -1.0, 'response_global_hybrid': 0, 'onsite_exchange_scale': 0.0,
    },
    'erf-tuned': {
        'lc_gamma': 'erf', 'omega': 0.25, 'cam_alpha': 0.0, 'cam_beta': 1.2,
        'lc_ground_state': True,
        'spc_coco': 0.5, 'spc_ovov': 0.5, 'spc_coov': 0.5,
        'c_mrsf': -1.0, 'response_global_hybrid': 0, 'onsite_exchange_scale': 0.0,
    },
    'dtcam-tb': {
        'lc_gamma': 'erf', 'omega': 0.3, 'cam_alpha': 0.0, 'cam_beta': 1.0,
        'lc_ground_state': False,
        'spc_coco': -1.0, 'spc_ovov': -1.0, 'spc_coov': -1.0,
        'c_mrsf': 1.02, 'response_global_hybrid': 1, 'onsite_exchange_scale': 1.60,
    },
}


def apply_tb_operator_preset(cfg_dict) -> None:
    """Expand ``tb_operator`` in a [dftb]/[xtb] section dict, in place.

    If ``cfg_dict['tb_operator']`` is a non-empty known preset name, every key
    of that preset is written into ``cfg_dict``. The preset WINS: it OVERRIDES
    any individually-set operator flag (lc_gamma, omega, cam_alpha, cam_beta,
    spc_*, lc_ground_state, c_mrsf, response_global_hybrid,
    onsite_exchange_scale) so the named operator is reproduced exactly. Set
    ``tb_operator`` to '' (the default) to keep hand-set individual flags.

    An unknown preset name raises ValueError. A no-op when tb_operator is unset
    or empty.
    """
    if not isinstance(cfg_dict, dict):
        return
    name = str(cfg_dict.get('tb_operator', '') or '').strip().lower()
    if not name:
        return
    if name not in TB_OPERATOR_PRESETS:
        raise ValueError(
            f"Unknown tb_operator {name!r}; expected one of "
            f"{sorted(TB_OPERATOR_PRESETS)} (or '' to use individual flags)."
        )
    for key, value in TB_OPERATOR_PRESETS[name].items():
        cfg_dict[key] = value
