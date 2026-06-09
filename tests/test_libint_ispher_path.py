from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def _compact_ws(text):
    return re.sub(r"\s+", " ", text)


def _no_ws(text):
    return re.sub(r"\s+", "", text)


def _compact_code(text):
    return re.sub(r"[\s&]", "", text)


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


def test_shellquartet_prefers_libint_for_pure_spherical_quartets():
    text = (ROOT / "source/integrals/int2.F90").read_text()
    body = _subroutine_body(text, "shellquartet")
    selector = _compact_ws(body.split("rotspd = max_am <= 2", 1)[0])
    assert _compact_code("libint = int2_quartet_prefers_libint(basis, eri_data%ids, max_am, eri_data%attenuated_ints, eri_data%rys_only)") in _compact_code(selector)
    assert re.search(r"rotspd\s*=\s*max_am\s*<=\s*2\s*\.and\.\s*\.not\.libint", body, re.IGNORECASE) is not None


def test_ints_exchange_prefers_libint_for_pure_spherical_quartets():
    text = (ROOT / "source/integrals/int2.F90").read_text()
    body = _subroutine_body(text, "ints_exchange")
    selector = _compact_ws(body.split("rotspd = max_am <= 2", 1)[0])
    assert _compact_code("libint = int2_quartet_prefers_libint(basis, shell_ids, max_am, attenuated, rys_only_)") in _compact_code(selector)
    assert re.search(r"rotspd\s*=\s*max_am\s*<=\s*2\s*\.and\.\s*\.not\.libint", body, re.IGNORECASE) is not None


def test_prefers_libint_helper_tracks_pure_5d_7f_shells():
    text = (ROOT / "source/integrals/int2.F90").read_text()
    start = text.find("logical function int2_quartet_prefers_libint")
    assert start >= 0, "helper not found"
    end = text.find("end function int2_quartet_prefers_libint", start)
    assert end >= 0, "helper end not found"
    body = _compact_ws(text[start:end])
    assert "if (HARMONIC_ACTIVE .and. allocated(basis%harmonic)) then" in body
    assert _compact_code("has_pure_shell = any(basis%harmonic(shell_ids) == 1 .and. basis%am(shell_ids) >= 2)") in _compact_code(body)
    assert _compact_code("(max_am > 2 .or. has_pure_shell)") in _compact_code(body)
    assert _compact_code(".not.attenuated") in _compact_code(body)


def test_libint_cxx_pure_bridge_is_optional_and_energy_only():
    cmake = (ROOT / "source/CMakeLists.txt").read_text()
    wrapper = (ROOT / "source/wrapper/libint_engine_wrapper.cpp").read_text()
    f90 = (ROOT / "source/integrals/int_libint.F90").read_text()
    int2 = (ROOT / "source/integrals/int2.F90").read_text()

    assert "find_path(EIGEN3_INCLUDE_DIR Eigen/Core" in cmake
    assert "OQP_LIBINT_CXX_ENGINE" in cmake
    assert "wrapper/libint_engine_wrapper.cpp" in cmake
    assert "libint2::Engine engine(libint2::Operator::coulomb" in wrapper
    assert "false};" in wrapper
    assert 'std::getenv("OQP_LIBINT_CXX_PURE_STATS")' in wrapper
    assert "std::atexit(print_stats)" in wrapper
    assert "Libint C++ pure ERI stats: attempts=" in wrapper
    assert "oqp_libint_engine_eri0" in f90
    assert 'get_environment_variable("OQP_LIBINT_CXX_PURE"' in f90
    assert 'trim(env_value) /= "1"' in f90
    assert 'get_environment_variable("OQP_LIBINT_CXX_PURE_STATS"' in f90
    assert "Libint C++ pure ERI stats: attempts=" in f90
    assert "libint_engine_attempts = libint_engine_attempts + 1_8" in f90
    assert "libint_engine_compute_eri(basis, eri_data%ids" in int2
    assert ".not. direct_libint" in int2
