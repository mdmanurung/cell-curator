"""Content-addressed artifact and receipt helpers."""

from __future__ import annotations

import json
import os
import platform
import sys
from pathlib import Path
from typing import Any, Literal

from .contracts import ContractError, sha256, values_sha256
from .models import BackendReceipt, Capability


def runtime_fingerprint() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "executable": sys.executable,
    }


def atomic_publish(path: Path, payload: bytes, *, resume: bool = True) -> str:
    """Publish bytes atomically, accepting an identical completed output on resume."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        observed = sha256(path)
        import hashlib

        expected = hashlib.sha256(payload).hexdigest()
        if resume and observed == expected:
            return observed
        raise ContractError(f"refusing to replace non-identical immutable artifact: {path}")
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        partial.write_bytes(payload)
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)
    return sha256(path)


def atomic_replace(path: Path, payload: bytes) -> str:
    """Atomically replace explicitly mutable audit state."""

    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.mutable-partial.{os.getpid()}")
    try:
        partial.write_bytes(payload)
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)
    return sha256(path)


def publish_json(path: Path, value: Any, *, resume: bool = True) -> str:
    return atomic_publish(
        path,
        (json.dumps(value, indent=2, sort_keys=True, default=str) + "\n").encode(),
        resume=resume,
    )


def replace_json(path: Path, value: Any) -> str:
    """Atomically replace a declared mutable JSON ledger or state file."""

    return atomic_replace(
        path,
        (json.dumps(value, indent=2, sort_keys=True, default=str) + "\n").encode(),
    )


def make_receipt(
    *,
    run_id: str,
    lane_key: str,
    capability: Capability,
    input_path: Path,
    config_path: Path,
    output_path: Path | None,
    artifacts: dict[str, Path] | None = None,
    feature_ids: list[str],
    groups: list[str],
    status: Literal["complete", "unavailable", "failed"] = "complete",
    limitations: tuple[str, ...] = (),
    runtime: dict[str, Any] | None = None,
    marker_methods: tuple[Literal["wilcoxon", "logreg"], ...] = (),
    convergence_validated: bool | None = None,
) -> BackendReceipt:
    return BackendReceipt(
        run_id=run_id,
        lane_key=lane_key,
        backend=capability.backend,
        capability_level=capability.level,
        status=status,
        input_sha256=sha256(input_path),
        config_sha256=sha256(config_path),
        output_sha256=sha256(output_path) if output_path and output_path.is_file() else None,
        artifact_sha256={
            name: sha256(path) for name, path in (artifacts or {}).items() if path.is_file()
        },
        feature_order_sha256=values_sha256(feature_ids),
        group_membership_sha256=values_sha256(groups),
        fallback_used=False,
        marker_methods=marker_methods,
        cpu_fallback_available=False,
        cpu_fallback_used=False,
        convergence_validated=convergence_validated,
        runtime=runtime or runtime_fingerprint(),
        limitations=limitations,
    )
