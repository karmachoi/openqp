"""Runtime helpers for experimental OpenQP GPU backends.

The GPU branches are intentionally split by target.  This helper stays pure
Python so input parsing and branch-level unit tests do not require CUDA
libraries or a built OpenQP shared library.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class GpuConfig:
    """Normalized GPU runtime configuration."""

    enabled: bool = False
    backend: str = "cuda"
    target: str = "metc"
    device: int = 0
    precision: str = "float64"
    fallback: str = "cpu"

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> "GpuConfig":
        """Create a normalized GPU config from an OpenQP config dictionary."""

        gpu = config.get("gpu", {}) or {}
        return cls(
            enabled=bool(gpu.get("enabled", False)),
            backend=str(gpu.get("backend", "cuda")).lower(),
            target=str(gpu.get("target", "metc")).lower(),
            device=int(gpu.get("device", 0)),
            precision=str(gpu.get("precision", "float64")).lower(),
            fallback=str(gpu.get("fallback", "cpu")).lower(),
        )

    @property
    def targets_metc(self) -> bool:
        """Return True when the requested target is METC contractions."""

        return self.target == "metc"

    @property
    def targets_xc_response(self) -> bool:
        """Return True when the requested target is TDHF/TDDFT XC response."""

        return self.target == "xc_response"

    def supports_xc_response(self, config: dict[str, Any]) -> bool:
        """Return True for the first scoped XC-response GPU validation path."""

        input_section = config.get("input", {}) or {}
        tdhf_section = config.get("tdhf", {}) or {}
        method = str(input_section.get("method", "hf")).lower()
        functional = str(input_section.get("functional", "")).strip()
        td_type = str(tdhf_section.get("type", "rpa")).lower()
        return (
            self.enabled
            and self.backend == "cuda"
            and self.target == "xc_response"
            and self.precision == "float64"
            and method == "tdhf"
            and td_type in {"rpa", "tda"}
            and bool(functional)
        )
