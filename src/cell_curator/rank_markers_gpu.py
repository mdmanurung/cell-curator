#!/usr/bin/env python3
"""Validate the native RAPIDS runtime contract with no CPU fallback."""

from __future__ import annotations

import inspect
import os
import socket
from importlib.metadata import distribution, version
from typing import Any

from .contracts import ContractError


class GpuContractError(ContractError):
    """Raised before evidence creation when the GPU contract is invalid."""


def validate_runtime_contract_values(
    *,
    observed_version: str,
    expected_version: str,
    observed_distribution: str,
    expected_distribution: str,
    device_count: int,
    device_name: str,
    compute_capability: str,
    allowed_gpu_contracts: list[dict],
    observed_cuda_major: int,
    expected_cuda_major: int,
) -> dict:
    if not expected_version or observed_version != expected_version:
        raise GpuContractError(
            f"rapids-singlecell version mismatch: {observed_version!r} != {expected_version!r}"
        )
    if not expected_distribution or observed_distribution != expected_distribution:
        raise GpuContractError("RAPIDS distribution does not match configuration")
    if device_count != 1:
        raise GpuContractError(f"exactly one visible GPU is required, observed {device_count}")
    if observed_cuda_major != expected_cuda_major:
        raise GpuContractError("CUDA runtime major does not match configuration")
    if not allowed_gpu_contracts:
        raise GpuContractError("allowed_gpu_contracts must not be empty")
    for contract in allowed_gpu_contracts:
        name = str(contract["name_contains"])
        capability = str(contract["compute_capability"])
        if name.lower() in device_name.lower() and capability == compute_capability:
            return {
                "name_contains": name,
                "compute_capability": capability,
            }
    raise GpuContractError(
        f"GPU model/compute capability is not allowed: {device_name!r} {compute_capability}"
    )


def installed_rapids_version() -> str:
    """Report the installed rapids-singlecell version for diagnostics.

    Used only in messages, so a missing distribution must not mask the error it
    is describing.
    """

    from importlib.metadata import PackageNotFoundError

    try:
        return version("rapids-singlecell")
    except PackageNotFoundError:
        return "not installed"


def native_ranking_callable(rsc: Any) -> Any:
    """Return ``rapids_singlecell.tl.rank_genes_groups``, or explain its absence.

    Native GPU Wilcoxon and logreg both come from this one entry point. It was
    introduced in rapids-singlecell 0.14; releases before that expose only
    ``rank_genes_groups_logreg`` and contain no Wilcoxon implementation, so the
    lane cannot satisfy its receipt contract on them. Fail with the installed
    version named rather than with a bare ``AttributeError``.
    """

    native = getattr(rsc.tl, "rank_genes_groups", None)
    if native is None:
        installed = installed_rapids_version()
        raise GpuContractError(
            "rapids_singlecell.tl.rank_genes_groups is required for native GPU "
            f"Wilcoxon and logreg, but rapids-singlecell {installed} does not "
            "provide it (that entry point arrived in 0.14). Install "
            "rapids-singlecell>=0.14, which needs Python 3.12 or newer."
        )
    return native


def supported_ranking_options(native: Any) -> frozenset[str]:
    """Return the keyword parameters the installed ranking entry point declares.

    Only declared parameters count. The entry point also takes ``**kwds``, which
    would swallow an unknown name without applying it — exactly the silent
    behaviour this contract exists to prevent.
    """

    parameters = inspect.signature(native).parameters
    return frozenset(
        name
        for name, item in parameters.items()
        if item.kind
        in {inspect.Parameter.KEYWORD_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD}
    )


def optional_ranking_requests(marker: dict[str, Any]) -> tuple[tuple[str, bool], ...]:
    """The optional ranking switches this configuration asks for.

    One definition, used both by the pre-flight contract check and by the call
    that builds the keyword arguments, so the two cannot disagree.
    """

    return (
        ("use_continuity", bool(marker["wilcoxon"]["continuity_correct"])),
        ("multi_gpu", False),
    )


def resolve_optional_ranking_option(
    *,
    name: str,
    requested: bool,
    supported: frozenset[str],
    installed_version: str,
) -> tuple[bool, str]:
    """Decide whether an optional statistical switch can be honoured.

    Returns ``(pass_it, note)``. A switch the installed version does not declare
    is only acceptable when it was requested off, because absent means "not
    applied" — that is the same result, and the note records it. Requesting it on
    against a version that cannot apply it fails closed instead of silently
    producing statistics the configuration did not ask for.
    """

    if name in supported:
        return True, f"{name}={requested} applied natively"
    if requested:
        raise GpuContractError(
            f"markers requests {name}=True, but rapids-singlecell "
            f"{installed_version} does not expose that parameter, so the option "
            "cannot be applied. Install a release that supports it, or set the "
            "option to false."
        )
    return False, f"{name} unavailable in {installed_version}; requested off, so not applied"


def format_compute_capability(raw: Any) -> str:
    """Render CuPy's packed compute capability as ``major.minor``.

    ``cupy.cuda.Device.compute_capability`` is a digit string in which only the
    final character is the minor version, so a two-digit major such as Blackwell
    ``sm_120`` must not be split positionally. Indexing ``[0]`` and ``[1]`` is
    correct only for single-digit majors and silently yields ``1.2`` for
    ``"120"``, which then compares against ``allowed_gpu_contracts`` as a
    different device.
    """

    if isinstance(raw, (tuple, list)):
        if len(raw) != 2:
            raise GpuContractError(f"unreadable compute capability: {raw!r}")
        return f"{int(raw[0])}.{int(raw[1])}"
    text = str(raw).strip()
    if "." in text:
        major, _, minor = text.partition(".")
        text = f"{major}{minor}"
    if not text.isdigit() or len(text) < 2:
        raise GpuContractError(f"unreadable compute capability: {raw!r}")
    return f"{int(text[:-1])}.{int(text[-1])}"


def inspect_runtime(config: dict[str, Any]) -> tuple[Any, Any, dict[str, Any]]:
    marker = config["markers"]
    if marker.get("cpu_fallback") is not False:
        raise GpuContractError("markers.cpu_fallback must be false")
    import cupy as cp
    import rapids_singlecell as rsc

    observed_version = version("rapids-singlecell")
    expected_distribution = str(marker["rapids_distribution"])
    package = distribution("rapids-singlecell")
    observed_distribution = str(package.metadata["Name"])
    device_count = int(cp.cuda.runtime.getDeviceCount())
    if device_count == 1:
        device = cp.cuda.Device(0)
        properties = cp.cuda.runtime.getDeviceProperties(0)
        raw_name = properties["name"]
        device_name = raw_name.decode() if isinstance(raw_name, bytes) else str(raw_name)
        capability = format_compute_capability(device.compute_capability)
    else:
        device_name, capability = "", ""
    runtime_version = int(cp.cuda.runtime.runtimeGetVersion())
    cuda_major = runtime_version // 1000
    validated = validate_runtime_contract_values(
        observed_version=observed_version,
        expected_version=str(marker["rapids_singlecell_version"]),
        observed_distribution=observed_distribution,
        expected_distribution=expected_distribution,
        device_count=device_count,
        device_name=device_name,
        compute_capability=capability,
        allowed_gpu_contracts=list(marker["allowed_gpu_contracts"]),
        observed_cuda_major=cuda_major,
        expected_cuda_major=int(marker["cuda_major"]),
    )
    native = native_ranking_callable(rsc)
    # Resolve the optional statistical switches here, before the potentially
    # large input is opened, so an unhonourable option fails closed early and the
    # decision is recorded in the receipt rather than inferred later.
    supported = supported_ranking_options(native)
    ranking_options = [
        resolve_optional_ranking_option(
            name=name,
            requested=requested,
            supported=supported,
            installed_version=observed_version,
        )[1]
        for name, requested in optional_ranking_requests(marker)
    ]
    return (
        cp,
        rsc,
        {
            "device_count": device_count,
            "device_name": device_name,
            "compute_capability": capability,
            "validated_gpu_contract": validated,
            "rapids_singlecell": observed_version,
            "rapids_distribution": observed_distribution,
            "cuda_major": cuda_major,
            "cuda_runtime_version": runtime_version,
            "cuda_driver_version": int(cp.cuda.runtime.driverGetVersion()),
            "native_api": f"{native.__module__}.{native.__name__}",
            "native_api_signature": str(inspect.signature(native)),
            "ranking_options": ranking_options,
            "hostname": socket.gethostname(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        },
    )
