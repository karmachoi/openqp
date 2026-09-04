"""Structural gates for exact batched MRSF XC-kernel derivatives."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/dftlib/dft_gridint_mrsf_xc_kernel_derivative.F90"
ENGINE = ROOT / "source/dftlib/dft_gridint_mrsf_xc_slice_gemm.F90"


def _body(text, name):
    start = text.index("  subroutine " + name)
    stop = text.index("  end subroutine " + name, start)
    return text[start:stop]


def test_point_workspace_is_reused_across_grid_points():
    text = SOURCE.read_text()
    point = _body(text, "kernel_point_potentials")
    assert "ws%" in point
    assert "allocate(" not in point
    assert "deallocate(" not in point


def test_all_density_and_response_fields_share_slice_level_dgemm_passes():
    text = SOURCE.read_text()
    slice_body = _body(text, "kernel_slice")
    assert slice_body.count("call slice_stack_values") == 2
    assert slice_body.count("call slice_stack_fixed") == 1
    assert slice_body.count("call slice_fock_derivative_accumulate") == 1
    point = _body(text, "kernel_point_potentials")
    assert "call total_density_derivative" not in point
    assert "moving_ao_pair_derivative" not in text
    assert "spin_fock_point_derivative" not in text


def test_engine_accumulates_with_stacked_dgemm_and_symmetrizes_once():
    text = ENGINE.read_text()
    accumulate = _body(text, "slice_fock_derivative_accumulate")
    assert accumulate.count("call dgemm(") == 2
    assert "n*ncart*nprobe" in accumulate
    values = _body(text, "slice_stack_values")
    assert values.count("call dgemm(") == 1
    fixed = _body(text, "slice_stack_fixed")
    assert fixed.count("call dgemm(") == 1
    source = SOURCE.read_text()
    assert source.count("call symmetrize_half_accumulator") == 1
