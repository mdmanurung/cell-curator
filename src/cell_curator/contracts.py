#!/usr/bin/env python3
"""Shared deterministic contracts for the bundled skill utilities."""

from __future__ import annotations

import base64
import hashlib
import math
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


def identifier_text(value: object) -> str:
    """Normalize an identifier cell, mapping the stringified null to absent.

    A missing identifier arrives as ``None`` or ``NaN``, and the obvious
    `str(value)` turns either into a short string that reads as data: the
    resulting text is truthy, so it survives an `if value` guard, is collected
    as a real identifier, and is written into durable, hashed artifacts as an
    accession the run never had. Anything consuming those artifacts later then
    fails on a term named "None", with nothing pointing at the stage that
    invented it.

    An accession is never the literal text `none`, `nan`, or `null`, so treating
    those as absent is unambiguous. Everything else is returned stripped, and a
    value that is already absent stays absent rather than becoming its own
    representation of absence.
    """

    text = str(value).strip()
    return "" if text.lower() in {"", "none", "nan", "null"} else text


def measured_float(value: object) -> float | None:
    """Parse a measurement cell, treating an unmeasured value as absent.

    A quantity that was never measured is `NaN` in memory, which pandas writes to
    a tab-separated file as an empty field, and a reader using
    `keep_default_na=False` gets that field back as the empty string. A bare
    `float()` then raises `could not convert string to float: ''`, so a run dies
    the first time an evidence table contains an unmeasured cell -- which is
    exactly the targeted-panel case, where a marker program's expected-negative
    violation cannot be measured because the panel does not carry the markers.

    Returning `None` instead lets a caller state what absence means. A gate that
    asks whether a value is within a bound should refuse an absent one, because
    "not measured" is not "measured as acceptable". The round trip is the only
    place the two spellings of absent disagree, so it is the only place this
    matters, and a test that builds its frame in memory cannot reach it.
    """

    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null", "na"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


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
