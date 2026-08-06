"""Zhu--Nakamura global-switching production and numerical gates."""

import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _native_zn(active=2, random_value=0.0, center_gap=0.05,
               threshold=0.1, same_force_sign=False, velocity_x=0.1):
    import oqp

    natom, nstate = 1, 2
    left = np.ascontiguousarray([[-1.0, 0.0, 0.0]])
    center = np.ascontiguousarray([[0.0, 0.0, 0.0]])
    right = np.ascontiguousarray([[1.0, 0.0, 0.0]])
    energy_left = np.ascontiguousarray([0.0, 0.2])
    energy_center = np.ascontiguousarray([0.0, center_gap])
    energy_right = np.ascontiguousarray([0.0, 0.2])
    # Cross interpolation gives F1=-0.1 and F2=+0.1 at the centre.
    gradient_left = np.ascontiguousarray([
        [[-0.1, 0.0, 0.0]],
        [[+0.1, 0.0, 0.0]],
    ])
    if same_force_sign:
        # Cross interpolation gives F1=-0.1 and F2=-0.2 (F1.F2 > 0).
        gradient_left = np.ascontiguousarray([
            [[+0.2, 0.0, 0.0]],
            [[+0.1, 0.0, 0.0]],
        ])
        gradient_right = np.ascontiguousarray([
            [[+0.1, 0.0, 0.0]],
            [[+0.2, 0.0, 0.0]],
        ])
    else:
        gradient_right = np.ascontiguousarray([
            [[+0.1, 0.0, 0.0]],
            [[-0.1, 0.0, 0.0]],
        ])
    mass = np.ascontiguousarray([1.0])
    velocity = np.ascontiguousarray([[velocity_x, 0.2, 0.0]])
    probabilities = np.zeros(nstate)
    a2 = np.zeros(nstate)
    b2 = np.zeros(nstate)
    hopped = np.zeros(1, dtype=np.int32)
    blocked = np.zeros(1, dtype=np.int32)
    target = np.array([active], dtype=np.int64)
    status = int(oqp.oqp_namd_zhu_nakamura_step(
        natom, nstate, active, threshold, random_value,
        oqp.ffi.cast('double *', left.ctypes.data),
        oqp.ffi.cast('double *', center.ctypes.data),
        oqp.ffi.cast('double *', right.ctypes.data),
        oqp.ffi.cast('double *', energy_left.ctypes.data),
        oqp.ffi.cast('double *', energy_center.ctypes.data),
        oqp.ffi.cast('double *', energy_right.ctypes.data),
        oqp.ffi.cast('double *', gradient_left.ctypes.data),
        oqp.ffi.cast('double *', gradient_right.ctypes.data),
        oqp.ffi.cast('double *', mass.ctypes.data),
        oqp.ffi.cast('double *', velocity.ctypes.data),
        oqp.ffi.cast('double *', probabilities.ctypes.data),
        oqp.ffi.cast('double *', a2.ctypes.data),
        oqp.ffi.cast('double *', b2.ctypes.data),
        oqp.ffi.cast('int *', hopped.ctypes.data),
        oqp.ffi.cast('int *', blocked.ctypes.data),
        oqp.ffi.cast('int64_t *', target.ctypes.data),
    ))
    return {
        'status': status, 'probabilities': probabilities, 'a2': a2,
        'b2': b2, 'hopped': int(hopped[0]), 'blocked': int(blocked[0]),
        'target': int(target[0]), 'velocity': velocity,
    }


def test_zhu_nakamura_path_is_fortran_resident_and_restart_aware():
    source = (ROOT / 'source/modules/namd.F90').read_text()
    header = (ROOT / 'include/oqp.h').read_text()
    wrapper = (ROOT / 'pyoqp/oqp/__init__.py').read_text()
    driver = (ROOT / 'pyoqp/oqp/library/namd.py').read_text()
    schema = (ROOT / 'pyoqp/oqp/molecule/oqpdata.py').read_text()

    assert 'subroutine namd_zhu_nakamura_pair(' in source
    assert 'bind(C, name="oqp_namd_zhu_nakamura_step")' in source
    assert 'int oqp_namd_zhu_nakamura_step(' in header
    assert repr('oqp_namd_zhu_nakamura_step') in wrapper
    assert 'def _run_zhu_nakamura(self):' in driver
    assert 'trial ``k`` point is discarded' in driver
    assert 'zn_history_gradients' in driver
    assert "'hop_method': {'type': str, 'default': 'fssh'}" in schema
    assert 'NAMD_RESTART_SCHEMA_VERSION = 9' in driver


def test_zhu_nakamura_preflight_rejects_unsupported_representations():
    from oqp.utils.input_checker import CheckReport, _check_runtype

    def diagnostics(md, qmmm=False):
        config = {
            'input': {
                'method': 'tdhf', 'runtype': 'namd', 'qmmm_flag': qmmm,
            },
            'tdhf': {'type': 'mrsf', 'nstate': 2},
            'md': {'active': 1, **md},
        }
        report = CheckReport()
        _check_runtype(config, report)
        return [item for item in report.errors if item.path == 'md.hop_method']

    assert not diagnostics({'hop_method': 'zhu_nakamura', 'soc': False})
    assert diagnostics({'hop_method': 'unpublished_guess', 'soc': False})
    assert diagnostics({'hop_method': 'zhu_nakamura', 'soc': True})
    assert diagnostics(
        {'hop_method': 'zhu_nakamura', 'soc': False}, qmmm=True)


def test_native_zhu_nakamura_probability_and_energy_conservation():
    result = _native_zn()
    assert result['status'] == 0
    assert result['hopped'] == 1
    assert result['blocked'] == 0
    assert result['target'] == 1

    expected_a2 = 0.5 * 0.2 * 0.1 / 0.05**3
    expected_b2 = (0.05 + 0.5*0.1**2 - 0.025) * 0.2 / (0.1*0.05)
    # F1.F2 < 0 selects the |b^4 - 1| branch of ZN eq. 3.
    radial = expected_b2 + math.sqrt(abs(expected_b2**2 - 1.0))
    expected_probability = math.exp(
        -math.pi/(4.0*math.sqrt(expected_a2)*math.sqrt(radial)))
    assert np.isclose(result['a2'][0], expected_a2, rtol=1.0e-13)
    assert np.isclose(result['b2'][0], expected_b2, rtol=1.0e-13)
    assert np.isclose(
        result['probabilities'][0], expected_probability, rtol=1.0e-13)

    before = 0.5*(0.1**2 + 0.2**2) + 0.05
    after = 0.5*np.sum(result['velocity']**2) + 0.0
    assert np.isclose(after, before, atol=1.0e-14)
    # Perpendicular y velocity is untouched by the self-consistent direction.
    assert result['velocity'][0, 1] == 0.2


def test_native_zhu_nakamura_probability_uses_positive_force_product_branch():
    result = _native_zn(active=1, random_value=np.nextafter(1.0, 0.0),
                        same_force_sign=True)
    expected_gamma = 0.1
    expected_omega = math.sqrt(0.02)
    expected_a2 = 0.5*expected_gamma*expected_omega/0.05**3
    expected_b2 = (0.0 + 0.5*0.1**2 - 0.025) * (
        expected_gamma/(expected_omega*0.05))
    radial = expected_b2 + math.sqrt(expected_b2**2 + 1.0)
    expected_probability = math.exp(
        -math.pi/(4.0*math.sqrt(expected_a2)*math.sqrt(radial)))

    assert result['status'] == 0
    assert result['hopped'] == 0
    assert np.isclose(result['a2'][1], expected_a2, rtol=1.0e-13)
    assert np.isclose(result['b2'][1], expected_b2, rtol=1.0e-13)
    assert np.isclose(
        result['probabilities'][1], expected_probability, rtol=1.0e-13)


def test_native_zhu_nakamura_requires_a_local_gap_minimum_inside_threshold():
    no_minimum = _native_zn(center_gap=0.25, threshold=1.0)
    assert no_minimum['status'] == 0
    assert no_minimum['hopped'] == 0
    assert np.all(no_minimum['probabilities'] == 0.0)

    gated = _native_zn(center_gap=0.05, threshold=0.01)
    assert gated['status'] == 0
    assert gated['hopped'] == 0
    assert np.all(gated['probabilities'] == 0.0)


def test_native_zhu_nakamura_blocks_an_energy_forbidden_upward_switch():
    result = _native_zn(active=1, random_value=0.0, velocity_x=0.25)
    assert result['status'] == 0
    assert result['hopped'] == 0
    assert result['blocked'] == 1
    assert result['target'] == 1
    assert np.array_equal(result['velocity'], [[0.25, 0.2, 0.0]])


def test_zhu_nakamura_restart_round_trip_preserves_three_point_history():
    from oqp.library.namd import NAMD

    driver = NAMD.__new__(NAMD)
    driver.hop_method = 'zhu_nakamura'
    driver.natom = 2
    driver.nstate = 3
    driver._zn_history = [
        {
            'step': step,
            'coordinates': np.full((2, 3), step + 0.125),
            'energies': np.arange(3, dtype=float) + step,
            'gradients': np.full((3, 2, 3), step + 0.25),
        }
        for step in (4, 5)
    ]

    saved = driver._restart_extra_payload()
    restored = driver._load_restart_extra(saved)
    driver._zn_history = []
    driver._restore_restart_extra(restored)

    assert [point['step'] for point in driver._zn_history] == [4, 5]
    for original, recovered in zip(
            ({
                'coordinates': np.full((2, 3), step + 0.125),
                'energies': np.arange(3, dtype=float) + step,
                'gradients': np.full((3, 2, 3), step + 0.25),
            } for step in (4, 5)), driver._zn_history):
        for key in ('coordinates', 'energies', 'gradients'):
            assert np.array_equal(recovered[key], original[key])
