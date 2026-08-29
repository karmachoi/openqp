from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_CMAKE = ROOT / "external" / "CMakeLists.txt"


def test_tagarray_uses_one_promoted_integer_interface():
    text = EXTERNAL_CMAKE.read_text(encoding="utf-8")
    assert 'set(_OQP_TAGARRAY_VERSION "1.0.0")' in text
    assert 'string(APPEND _TAGARRAY_Fortran_FLAGS " -fdefault-integer-8")' in text
    assert "libtagarray-v0.0.6-default-integer-shapes.patch" not in text
    assert "fix_tagarray_linelen.py" in text
