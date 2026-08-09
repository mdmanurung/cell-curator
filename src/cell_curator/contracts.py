#!/usr/bin/env python3
"""Shared deterministic contracts for the bundled skill utilities."""

from __future__ import annotations

import base64
import hashlib
import re
from pathlib import Path


class ContractError(RuntimeError):
    """Raised when an artifact violates a reusable workflow contract."""


def safe_path_component(value: object) -> str:
    """Return a collision-free path component while preserving already-safe values."""

    rendered = str(value)
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", rendered):
        return rendered
    encoded = base64.urlsafe_b64encode(rendered.encode("utf-8")).decode("ascii").rstrip("=")
    # Safe literal identifiers cannot start with ``_``, keeping the encoded
    # namespace disjoint from every value preserved by the branch above.
    return f"_u-{encoded}"


def require_safe_path_component(value: object, *, field: str) -> str:
    """Reject a semantic identifier that would escape or alias a path component."""

    rendered = str(value)
    if safe_path_component(rendered) != rendered:
        raise ContractError(f"{field} is not a safe path component: {rendered!r}")
    return rendered


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def values_sha256(values: list[str], *, sort_values: bool = False) -> str:
    normalized = sorted(map(str, values)) if sort_values else list(map(str, values))
    digest = hashlib.sha256()
    for value in normalized:
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no", ""}:
        return False
    raise ContractError(f"cannot interpret Boolean value {value!r}")


def require_columns(frame, columns: set[str], *, artifact: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ContractError(f"{artifact} misses columns: {sorted(missing)}")
