from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def _subroutine_body(text, name):
    match = re.search(
        rf"subroutine {name}\b(?P<body>.*?)end subroutine {name}",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    assert match, f"{name} not found"
    return match.group("body")


def test_shellquartet_libint_normalizes_with_libint_shell_order():
    text = (ROOT / "source/integrals/int2.F90").read_text()
    body = _subroutine_body(text, "shellquartet")

    libint_branch = body.split("else if (libint) then", 1)[1].split("else if (rys) then", 1)[0]
    assert "call normalize_ints(nbf, eri_data%am(eri_data%flips), eri_data%pints)" in libint_branch
    assert "eri_data%gdat%am" not in libint_branch


def test_ints_exchange_libint_uses_cartesian_counts_before_spherical_schwarz():
    text = (ROOT / "source/integrals/int2.F90").read_text()
    body = _subroutine_body(text, "ints_exchange")
    libint_branch = body.split("else if (libint) then", 1)[1].split("else if (rys) then", 1)[0]

    assert "nbf = NUM_CART_BF(am(flips))" in libint_branch
    assert "nbf = am(flips)" not in libint_branch
    assert "call normalize_ints(nbf, am(flips), pints)" in libint_branch


def test_ints_exchange_libint_schwarz_matches_spherical_runtime_path():
    text = (ROOT / "source/integrals/int2.F90").read_text()
    body = _subroutine_body(text, "ints_exchange")

    assert "use constants, only: HARMONIC_ACTIVE, NUM_CART_BF" in body
    assert "use cart2sph, only: cart2sph_eri" in body

    libint_branch = body.split("else if (libint) then", 1)[1].split("else if (rys) then", 1)[0]
    assert "if (zero_shq) then" in libint_branch
    assert "if (HARMONIC_ACTIVE) then" in libint_branch
    assert "pure_s(s_) = basis%harmonic(shell_ids(orig_))" in libint_branch
    assert "call cart2sph_eri(ints, am_s, pure_s, nbf_s, nbf_out_s)" in libint_branch
    assert "vmax = maxval(abs(ints(1:product(nbf_out_s))))" in libint_branch
