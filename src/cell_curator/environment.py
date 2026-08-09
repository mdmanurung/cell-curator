"""Read-only compute-environment detection and fail-closed CUDA validation.

Detection never installs software or submits scheduler jobs.  A visible driver,
``CUDA_VISIBLE_DEVICES``, or ``nvidia-smi`` binary is only a CUDA *candidate*;
GPU-required work is permitted only after a small CuPy allocation and kernel have
completed on a real device.
"""

from __future__ import annotations

import importlib
import os
import platform
import re
import resource
import shutil
from collections.abc import Callable, Mapping
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal

from pydantic import Field

from .models import StrictModel


class ComputeProfile(StrEnum):
    """Executable compute profiles understood by routing and orchestration."""

    LOCAL_CPU = "local_cpu"
    LOCAL_GPU = "local_gpu"
    SCHEDULER = "scheduler"


class GPUValidationError(RuntimeError):
    """Raised when a CUDA method cannot prove that a real GPU is usable."""


class ResourceLimit(StrictModel):
    """One detected or configured resource ceiling."""

    source: str
    kind: Literal["memory", "cpu"]
    value: int


class SchedulerDetection(StrictModel):
    """Read-only evidence that a batch scheduler is present or active."""

    kind: Literal["slurm", "pbs", "lsf", "sge"]
    detected_by: tuple[Literal["environment", "path"], ...]
    active_job: bool
    job_id: str | None = None
    commands: tuple[str, ...] = ()
    environment_keys: tuple[str, ...] = ()


class GPUValidation(StrictModel):
    """Proof that a CUDA runtime executed work on a real visible device."""

    status: Literal["validated"] = "validated"
    runtime_backend: Literal["cupy"] = "cupy"
    runtime_version: int
    driver_version: int
    package_version: str
    device_count: int
    device_index: int
    device_name: str
    compute_capability: str
    total_memory_bytes: int | None = None
    cuda_visible_devices: str | None = None
    proof: str = "allocation-kernel-synchronization-host-transfer"
    cpu_fallback_allowed: Literal[False] = False


class EnvironmentReport(StrictModel):
    """Structured, auditable result of environment detection."""

    system: str
    machine: str
    apple_silicon: bool
    accelerator: Literal["none", "apple_mps", "cuda"]
    cpu_count: int
    usable_cpu_count: int
    total_ram_bytes: int
    usable_ram_bytes: int
    resource_limits: tuple[ResourceLimit, ...] = ()
    scheduler: SchedulerDetection | None = None
    cuda_candidate: bool = False
    cuda_candidate_reasons: tuple[str, ...] = ()
    gpu_validation: GPUValidation | None = None
    gpu_validation_error: str | None = None
    compute_profile: ComputeProfile
    profile_candidates: tuple[ComputeProfile, ...]
    requires_confirmation: bool = False
    ambiguity_reason: str | None = None
    blockers: tuple[str, ...] = ()
    next_action: str
    configured_cpu_resources: dict[str, Any] = Field(default_factory=dict)
    configured_gpu_resources: dict[str, Any] = Field(default_factory=dict)
    configured_gpu_requested: bool = False
    detection_only: Literal[True] = True

    @property
    def gpu_validated(self) -> bool:
        """Whether GPU-required methods may execute without a CPU substitution."""

        return self.gpu_validation is not None


_SCHEDULERS: tuple[tuple[str, tuple[str, ...], tuple[str, ...], tuple[str, ...]], ...] = (
    (
        "slurm",
        ("SLURM_JOB_ID",),
        ("SLURM_JOB_ID", "SLURM_CLUSTER_NAME", "SLURM_SUBMIT_HOST"),
        ("sbatch", "squeue"),
    ),
    (
        "pbs",
        ("PBS_JOBID",),
        ("PBS_JOBID", "PBS_ENVIRONMENT", "PBS_O_HOST"),
        ("qsub", "qstat"),
    ),
    (
        "lsf",
        ("LSB_JOBID",),
        ("LSB_JOBID", "LSF_ENVDIR"),
        ("bsub", "bjobs"),
    ),
    (
        "sge",
        ("JOB_ID",),
        ("JOB_ID", "SGE_ROOT"),
        ("qsub", "qstat"),
    ),
)


def _read_text(path: Path) -> str | None:
    try:
        return path.read_text().strip()
    except (OSError, UnicodeError):
        return None


def _positive_int(value: str | int | None) -> int | None:
    try:
        parsed = int(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


def _parse_memory(value: str | None, *, default_unit: int = 1024**2) -> int | None:
    """Parse scheduler memory values; unqualified values conventionally mean MiB."""

    if not value:
        return None
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)\s*([KMGTPE]?)(?:I?B)?\s*", value, re.I)
    if not match:
        return None
    amount = float(match.group(1))
    suffix = match.group(2).upper()
    multiplier = default_unit if not suffix else 1024 ** ("KMGTPE".index(suffix) + 1)
    parsed = int(amount * multiplier)
    return parsed if parsed > 0 else None


def _parse_cpu_set(value: str | None) -> int | None:
    """Count CPUs in Linux cpuset syntax such as ``0-3,8,10-11``."""

    if not value:
        return None
    members: set[int] = set()
    try:
        for segment in value.split(","):
            bounds = segment.strip().split("-", maxsplit=1)
            start = int(bounds[0])
            stop = int(bounds[-1])
            if start < 0 or stop < start:
                return None
            members.update(range(start, stop + 1))
    except ValueError:
        return None
    return len(members) or None


def _configured_limits(resources: Mapping[str, Any]) -> tuple[int | None, int | None]:
    """Read common Snakemake/HPC CPU and host-memory resource spellings."""

    cpu_values: list[int] = []
    memory_values: list[int] = []
    for key, raw_value in resources.items():
        normalized = key.strip().casefold().replace("-", "_")
        if normalized in {
            "cpu",
            "cpus",
            "cores",
            "threads",
            "cpu_count",
            "cpus_per_task",
        }:
            value = _positive_int(raw_value)
            if value is None:
                raise ValueError(f"configured CPU resource {key!r} must be a positive integer")
            cpu_values.append(value)
        elif normalized in {"memory_bytes", "mem_bytes"}:
            value = _positive_int(raw_value)
            if value is None:
                raise ValueError(f"configured memory resource {key!r} must be positive")
            memory_values.append(value)
        elif normalized in {"memory_mb", "mem_mb", "memory_mib", "mem_mib"}:
            value = _positive_int(raw_value)
            if value is None:
                raise ValueError(f"configured memory resource {key!r} must be positive")
            memory_values.append(value * 1024**2)
        elif normalized in {"memory_gb", "mem_gb", "memory_gib", "mem_gib"}:
            value = _positive_int(raw_value)
            if value is None:
                raise ValueError(f"configured memory resource {key!r} must be positive")
            memory_values.append(value * 1024**3)
    return (
        min(cpu_values) if cpu_values else None,
        min(memory_values) if memory_values else None,
    )


def _gpu_requested(resources: Mapping[str, Any]) -> bool:
    """Recognize a configured GPU request without treating it as hardware proof."""

    for key, raw_value in resources.items():
        normalized = key.strip().casefold().replace("-", "_")
        if normalized in {"gpu", "gpus", "gpu_count", "gpus_per_task"}:
            return (_positive_int(raw_value) or 0) > 0
        if normalized in {"gres", "accelerator", "device", "partition"}:
            text = str(raw_value).strip().casefold()
            if "gpu" in text or "cuda" in text or "nvidia" in text:
                return True
    return False


def _host_total_ram_bytes() -> int:
    try:
        pages = int(os.sysconf("SC_PHYS_PAGES"))
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        if pages > 0 and page_size > 0:
            return pages * page_size
    except (OSError, TypeError, ValueError):
        pass
    text = _read_text(Path("/proc/meminfo")) or ""
    match = re.search(r"^MemTotal:\s+(\d+)\s+kB$", text, re.MULTILINE)
    if match:
        return int(match.group(1)) * 1024
    raise RuntimeError("could not determine host RAM from sysconf or /proc/meminfo")


def _resource_limits(
    environ: Mapping[str, str],
    *,
    total_ram_bytes: int,
    configured_memory_bytes: int | None,
    configured_cpu_count: int | None,
) -> list[ResourceLimit]:
    limits: list[ResourceLimit] = []
    for source, path in (
        ("cgroup_v2", Path("/sys/fs/cgroup/memory.max")),
        ("cgroup_v1", Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")),
    ):
        raw = _read_text(path)
        value = _positive_int(raw) if raw not in {None, "max"} else None
        # Kernels often expose an enormous sentinel instead of infinity.
        if value is not None and value < 2**60:
            limits.append(ResourceLimit(source=source, kind="memory", value=value))

    try:
        affinity_count = len(os.sched_getaffinity(0))
        if affinity_count > 0:
            limits.append(
                ResourceLimit(source="process_affinity", kind="cpu", value=affinity_count)
            )
    except (AttributeError, OSError):
        pass

    for source, path in (
        ("cgroup_v2_cpuset", Path("/sys/fs/cgroup/cpuset.cpus.effective")),
        ("cgroup_v1_cpuset", Path("/sys/fs/cgroup/cpuset/cpuset.cpus")),
    ):
        value = _parse_cpu_set(_read_text(path))
        if value is not None:
            limits.append(ResourceLimit(source=source, kind="cpu", value=value))

    cpu_max = (_read_text(Path("/sys/fs/cgroup/cpu.max")) or "").split()
    if len(cpu_max) == 2 and cpu_max[0] != "max":
        quota = _positive_int(cpu_max[0])
        period = _positive_int(cpu_max[1])
        if quota is not None and period is not None:
            limits.append(
                ResourceLimit(
                    source="cgroup_v2_quota",
                    kind="cpu",
                    value=max(1, quota // period),
                )
            )
    quota = _positive_int(_read_text(Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")))
    period = _positive_int(_read_text(Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")))
    if quota is not None and period is not None:
        limits.append(
            ResourceLimit(source="cgroup_v1_quota", kind="cpu", value=max(1, quota // period))
        )

    try:
        for name, key in (
            ("rlimit_as", resource.RLIMIT_AS),
            ("rlimit_data", resource.RLIMIT_DATA),
        ):
            soft, _hard = resource.getrlimit(key)
            if soft not in {resource.RLIM_INFINITY, -1} and 0 < int(soft) < 2**60:
                limits.append(ResourceLimit(source=name, kind="memory", value=int(soft)))
    except (OSError, ValueError):
        pass

    slurm_node = _parse_memory(environ.get("SLURM_MEM_PER_NODE"))
    if slurm_node:
        limits.append(
            ResourceLimit(source="slurm_mem_per_node", kind="memory", value=slurm_node)
        )
    slurm_per_cpu = _parse_memory(environ.get("SLURM_MEM_PER_CPU"))
    slurm_cpus = _positive_int(
        environ.get("SLURM_CPUS_ON_NODE") or environ.get("SLURM_CPUS_PER_TASK")
    )
    if slurm_per_cpu and slurm_cpus:
        limits.append(
            ResourceLimit(
                source="slurm_mem_per_cpu", kind="memory", value=slurm_per_cpu * slurm_cpus
            )
        )
    for source, environment_key in (
        ("pbs_memory", "PBS_RESC_MEM"),
        ("configured_environment", "CELL_CURATOR_MEMORY_BYTES"),
    ):
        value = (
            _parse_memory(environ.get(environment_key))
            if environment_key == "PBS_RESC_MEM"
            else _positive_int(environ.get(environment_key))
        )
        if value:
            limits.append(ResourceLimit(source=source, kind="memory", value=value))
    if configured_memory_bytes is not None:
        if configured_memory_bytes <= 0:
            raise ValueError("configured_memory_bytes must be positive")
        limits.append(
            ResourceLimit(
                source="configured_resource_limit", kind="memory", value=configured_memory_bytes
            )
        )

    for source, environment_key in (
        ("slurm_cpus", "SLURM_CPUS_ON_NODE"),
        ("slurm_cpus_per_task", "SLURM_CPUS_PER_TASK"),
        ("pbs_cpus", "PBS_NCPUS"),
        ("lsf_cpus", "LSB_DJOB_NUMPROC"),
        ("sge_slots", "NSLOTS"),
    ):
        value = _positive_int(environ.get(environment_key))
        if value:
            limits.append(ResourceLimit(source=source, kind="cpu", value=value))
    if configured_cpu_count is not None:
        if configured_cpu_count <= 0:
            raise ValueError("configured_cpu_count must be positive")
        limits.append(
            ResourceLimit(
                source="configured_resource_limit", kind="cpu", value=configured_cpu_count
            )
        )
    # Keep the host total visible as the ultimate physical bound without calling it
    # a configured/cgroup limit.
    limits.append(ResourceLimit(source="host_physical", kind="memory", value=total_ram_bytes))
    return sorted(limits, key=lambda item: (item.kind, item.source, item.value))


def _detect_scheduler(
    environ: Mapping[str, str], which: Callable[[str], str | None]
) -> SchedulerDetection | None:
    detections: list[SchedulerDetection] = []
    for kind, job_keys, environment_keys, commands in _SCHEDULERS:
        present_env = tuple(key for key in environment_keys if environ.get(key))
        present_commands = tuple(command for command in commands if which(command))
        if not present_env and not present_commands:
            continue
        detected_by: list[Literal["environment", "path"]] = []
        if present_env:
            detected_by.append("environment")
        if present_commands:
            detected_by.append("path")
        job_id = next((environ[key] for key in job_keys if environ.get(key)), None)
        detections.append(
            SchedulerDetection(
                kind=kind,  # type: ignore[arg-type]
                detected_by=tuple(detected_by),
                active_job=job_id is not None,
                job_id=job_id,
                commands=present_commands,
                environment_keys=present_env,
            )
        )
    if not detections:
        return None
    # Active scheduler state is stronger than PATH alone; ties use a stable kind order.
    order = {item[0]: index for index, item in enumerate(_SCHEDULERS)}
    return sorted(detections, key=lambda item: (not item.active_job, order[item.kind]))[0]


def _cuda_candidate(
    environ: Mapping[str, str], which: Callable[[str], str | None]
) -> tuple[bool, tuple[str, ...]]:
    reasons: list[str] = []
    visible = environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None and visible.strip() not in {"", "-1", "none", "None"}:
        reasons.append("CUDA_VISIBLE_DEVICES names one or more devices")
    if Path("/dev/nvidiactl").exists() and any(
        Path(f"/dev/nvidia{index}").exists() for index in range(16)
    ):
        reasons.append("NVIDIA device nodes are present")
    gpu_root = Path("/proc/driver/nvidia/gpus")
    try:
        if gpu_root.is_dir() and any(gpu_root.iterdir()):
            reasons.append("the NVIDIA driver reports GPU entries")
    except OSError:
        pass
    if which("nvidia-smi"):
        reasons.append("nvidia-smi is on PATH")
    return bool(reasons), tuple(reasons)


def _property(properties: Mapping[Any, Any], name: str, default: Any = None) -> Any:
    return properties.get(name, properties.get(name.encode(), default))


def validate_real_gpu(
    *,
    cupy_module: Any | None = None,
    environ: Mapping[str, str] | None = None,
) -> GPUValidation:
    """Prove that a CUDA device can execute a kernel; never returns a CPU fallback."""

    environment = dict(os.environ if environ is None else environ)
    visible = environment.get("CUDA_VISIBLE_DEVICES")
    if visible is not None and visible.strip() in {"", "-1", "none", "None"}:
        raise GPUValidationError(
            "CUDA validation failed: CUDA_VISIBLE_DEVICES disables all devices. "
            "Run inside an authorized GPU allocation with a visible NVIDIA GPU; "
            "CPU substitution is disabled for GPU-required methods."
        )
    try:
        cp = cupy_module or importlib.import_module("cupy")
        if hasattr(cp.cuda, "is_available") and not bool(cp.cuda.is_available()):
            raise RuntimeError("CuPy reports that CUDA is unavailable")
        device_count = int(cp.cuda.runtime.getDeviceCount())
        if device_count < 1:
            raise RuntimeError("the CUDA runtime reports zero visible devices")
        device_index = 0
        device = cp.cuda.Device(device_index)
        properties = cp.cuda.runtime.getDeviceProperties(device_index)
        raw_name = _property(properties, "name", "unknown")
        device_name = raw_name.decode() if isinstance(raw_name, bytes) else str(raw_name)
        major = int(_property(properties, "major", device.compute_capability[0]))
        minor = int(_property(properties, "minor", device.compute_capability[1]))
        total_memory = _property(properties, "totalGlobalMem")
        probe = cp.asarray([1.0, 2.0], dtype=cp.float32)
        observed = float(cp.asnumpy((probe * probe).sum()))
        cp.cuda.Stream.null.synchronize()
        if abs(observed - 5.0) > 1e-6:
            raise RuntimeError("the CUDA validation kernel returned an invalid result")
        return GPUValidation(
            runtime_version=int(cp.cuda.runtime.runtimeGetVersion()),
            driver_version=int(cp.cuda.runtime.driverGetVersion()),
            package_version=str(getattr(cp, "__version__", "unknown")),
            device_count=device_count,
            device_index=device_index,
            device_name=device_name,
            compute_capability=f"{major}.{minor}",
            total_memory_bytes=int(total_memory) if total_memory is not None else None,
            cuda_visible_devices=visible,
        )
    except GPUValidationError:
        raise
    except Exception as exc:
        raise GPUValidationError(
            f"CUDA validation failed: {exc}. Run inside an authorized GPU allocation "
            "with a compatible CuPy/CUDA runtime and a visible NVIDIA GPU; CPU "
            "substitution is disabled for GPU-required methods."
        ) from exc


def require_real_gpu(
    *,
    cupy_module: Any | None = None,
    environ: Mapping[str, str] | None = None,
) -> GPUValidation:
    """Require runtime CUDA proof for a GPU-only operation, failing closed."""

    return validate_real_gpu(cupy_module=cupy_module, environ=environ)


def detect_environment(
    *,
    configured_memory_bytes: int | None = None,
    configured_cpu_count: int | None = None,
    configured_cpu_resources: Mapping[str, Any] | None = None,
    configured_gpu_resources: Mapping[str, Any] | None = None,
    validate_cuda: bool = True,
    environ: Mapping[str, str] | None = None,
    which: Callable[[str], str | None] | None = None,
    system: str | None = None,
    machine: str | None = None,
    cpu_count: int | None = None,
    total_ram_bytes: int | None = None,
    gpu_validator: Callable[[], GPUValidation] | None = None,
) -> EnvironmentReport:
    """Detect compute resources without mutation, installation, or job submission."""

    environment = dict(os.environ if environ is None else environ)
    cpu_resources = dict(configured_cpu_resources or {})
    gpu_resources = dict(configured_gpu_resources or {})
    mapped_cpu_count, mapped_memory_bytes = _configured_limits(cpu_resources)
    if mapped_cpu_count is not None:
        configured_cpu_count = min(
            value for value in (configured_cpu_count, mapped_cpu_count) if value is not None
        )
    if mapped_memory_bytes is not None:
        configured_memory_bytes = min(
            value
            for value in (configured_memory_bytes, mapped_memory_bytes)
            if value is not None
        )
    gpu_requested = _gpu_requested(gpu_resources)
    command_lookup = shutil.which if which is None else which
    detected_system = system or platform.system()
    detected_machine = machine or platform.machine()
    host_cpus = cpu_count if cpu_count is not None else (os.cpu_count() or 1)
    if host_cpus <= 0:
        raise ValueError("cpu_count must be positive")
    host_ram = total_ram_bytes if total_ram_bytes is not None else _host_total_ram_bytes()
    if host_ram <= 0:
        raise ValueError("total_ram_bytes must be positive")
    limits = _resource_limits(
        environment,
        total_ram_bytes=host_ram,
        configured_memory_bytes=configured_memory_bytes,
        configured_cpu_count=configured_cpu_count,
    )
    memory_limits = [item.value for item in limits if item.kind == "memory"]
    cpu_limits = [item.value for item in limits if item.kind == "cpu"]
    usable_ram = min(memory_limits)
    usable_cpus = min([host_cpus, *cpu_limits])

    scheduler = _detect_scheduler(environment, command_lookup)
    apple_silicon = detected_system == "Darwin" and detected_machine.lower() in {
        "arm64",
        "aarch64",
    }
    candidate, candidate_reasons = _cuda_candidate(environment, command_lookup)
    validation: GPUValidation | None = None
    validation_error: str | None = None
    blockers: list[str] = []
    if candidate and validate_cuda:
        try:
            validation = (gpu_validator or (lambda: validate_real_gpu(environ=environment)))()
        except GPUValidationError as exc:
            validation_error = str(exc)
            blockers.append(str(exc))
    if gpu_requested and validation is None and scheduler is None:
        blockers.append(
            "configured GPU resource request is not validated by a real CUDA GPU; "
            "CPU substitution is disabled for GPU-required methods"
        )

    candidates: tuple[ComputeProfile, ...]
    if scheduler is not None and scheduler.active_job:
        profile = ComputeProfile.SCHEDULER
        candidates = (ComputeProfile.SCHEDULER,)
        requires_confirmation = False
        ambiguity = None
    elif validation is not None:
        profile = ComputeProfile.LOCAL_GPU
        candidates = (
            (ComputeProfile.LOCAL_GPU, ComputeProfile.SCHEDULER)
            if scheduler is not None
            else (ComputeProfile.LOCAL_GPU,)
        )
        requires_confirmation = False
        ambiguity = None
    elif scheduler is not None:
        profile = ComputeProfile.SCHEDULER
        candidates = (ComputeProfile.SCHEDULER, ComputeProfile.LOCAL_CPU)
        requires_confirmation = True
        ambiguity = (
            f"{scheduler.kind} is visible but no active job or validated local CUDA device was "
            "found; this may be a login/submit node or an ordinary CPU host"
        )
        blockers.append(ambiguity)
    else:
        profile = ComputeProfile.LOCAL_CPU
        candidates = (ComputeProfile.LOCAL_CPU,)
        requires_confirmation = False
        ambiguity = None

    if requires_confirmation:
        next_action = (
            "Confirm scheduler versus local_cpu before routing compute; detection did not "
            "submit a job."
        )
    elif profile is ComputeProfile.SCHEDULER:
        next_action = (
            "Use the detected scheduler profile only through an explicitly authorized job "
            "submission workflow."
        )
    elif profile is ComputeProfile.LOCAL_GPU:
        next_action = "CUDA is validated; local GPU methods may be routed explicitly."
    elif candidate and validation is None:
        next_action = (
            "CUDA indicators were found but runtime validation did not pass; use local_cpu "
            "or obtain an authorized validated GPU environment."
        )
    else:
        next_action = "Use the local_cpu profile within the reported RAM and CPU limits."

    accelerator: Literal["none", "apple_mps", "cuda"]
    if validation is not None:
        accelerator = "cuda"
    elif apple_silicon:
        accelerator = "apple_mps"
    else:
        accelerator = "none"
    return EnvironmentReport(
        system=detected_system,
        machine=detected_machine,
        apple_silicon=apple_silicon,
        accelerator=accelerator,
        cpu_count=host_cpus,
        usable_cpu_count=usable_cpus,
        total_ram_bytes=host_ram,
        usable_ram_bytes=usable_ram,
        resource_limits=tuple(limits),
        scheduler=scheduler,
        cuda_candidate=candidate,
        cuda_candidate_reasons=candidate_reasons,
        gpu_validation=validation,
        gpu_validation_error=validation_error,
        compute_profile=profile,
        profile_candidates=candidates,
        requires_confirmation=requires_confirmation,
        ambiguity_reason=ambiguity,
        blockers=tuple(dict.fromkeys(blockers)),
        next_action=next_action,
        configured_cpu_resources=cpu_resources,
        configured_gpu_resources=gpu_resources,
        configured_gpu_requested=gpu_requested,
    )


__all__ = [
    "ComputeProfile",
    "EnvironmentReport",
    "GPUValidation",
    "GPUValidationError",
    "ResourceLimit",
    "SchedulerDetection",
    "detect_environment",
    "require_real_gpu",
    "validate_real_gpu",
]
