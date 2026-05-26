"""Pure-Python persistent METC allocation runtime scaffolding.

This module deliberately does not allocate CUDA memory.  It is a source-level
contract for future Fortran/CUDA wiring: a caller must provide the allocation
table generated from :mod:`gpu_metc_buffers`, and the registry validates that
ABI table before remembering the reusable plan metadata.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


AllocationTable = tuple[tuple[int, str, int, str], ...]
ReuseKey = tuple[int, int, int, int, int]


@dataclass(frozen=True)
class PersistentMetcAllocationRecord:
    """Validated metadata for one reusable METC allocation shape."""

    reuse_key: ReuseKey
    total_bytes: int
    table: AllocationTable


class PersistentMetcAllocationRegistry:
    """Record validated persistent-buffer plans without touching CUDA."""

    def __init__(self) -> None:
        self._records: dict[ReuseKey, PersistentMetcAllocationRecord] = {}

    def register_plan(self, plan: Any, table: AllocationTable) -> PersistentMetcAllocationRecord:
        """Validate and remember a persistent METC allocation plan.

        ``plan`` is intentionally duck-typed so source-level tests can load this
        module without importing the full OpenQP package.  The real runtime path
        is expected to pass a ``PersistentMetcBufferPlan`` instance.
        """

        validated_table = plan.validate_fortran_allocation_table(table)
        record = PersistentMetcAllocationRecord(
            reuse_key=plan.reuse_key,
            total_bytes=plan.total_bytes,
            table=validated_table,
        )
        self._records[record.reuse_key] = record
        return record

    def lookup(self, reuse_key: ReuseKey) -> PersistentMetcAllocationRecord | None:
        """Return the validated record for ``reuse_key`` if one is registered."""

        return self._records.get(reuse_key)
