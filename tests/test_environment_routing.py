from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

import cell_curator.environment as environment_module
from cell_curator.environment import (
    ComputeProfile,
    GPUValidation,
    GPUValidationError,
    ResourceLimit,
    SchedulerDetection,
    detect_environment,
    require_real_gpu,
    validate_real_gpu,
)
from cell_curator.routing import (
    ALL_METHODS,
    CrosscheckMethod,
    HealthContext,
    MethodResource,
    RoutingRequest,
    Technology,
    route_crosscheck,
    route_hierarchy,
)

GIB = 1024**3


def _fixed_limits(
    *_args,
    total_ram_bytes: int,
    configured_memory_bytes=None,
    configured_cpu_count=None,
    **_kwargs,
):
    limits = [ResourceLimit(source="host_physical", kind="memory", value=total_ram_bytes)]
    if configured_memory_bytes:
        limits.append(
            ResourceLimit(
                source="configured_resource_limit", kind="memory", value=configured_memory_bytes
            )
        )
    if configured_cpu_count:
        limits.append(
            ResourceLimit(
                source="configured_resource_limit", kind="cpu", value=configured_cpu_count
            )
        )
    return limits


def _validated_gpu() -> GPUValidation:
    return GPUValidation(
        runtime_version=12040,
        driver_version=12060,
        package_version="13.0.0",
        device_count=1,
        device_index=0,
        device_name="Fixture GPU",
        compute_capability="8.0",
        total_memory_bytes=40 * GIB,
    )


def _isolate_detection(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(environment_module, "_resource_limits", _fixed_limits)
    monkeypatch.setattr(environment_module, "_cuda_candidate", lambda *_args: (False, ()))
    monkeypatch.setattr(environment_module, "_detect_scheduler", lambda *_args: None)


def test_local_cpu_and_configured_resource_limits(monkeypatch: pytest.MonkeyPatch) -> None:
    _isolate_detection(monkeypatch)
    report = detect_environment(
        environ={},
        system="Linux",
        machine="x86_64",
        cpu_count=32,
        total_ram_bytes=64 * GIB,
        configured_memory_bytes=12 * GIB,
        configured_cpu_count=6,
    )
    assert report.compute_profile is ComputeProfile.LOCAL_CPU
    assert report.usable_ram_bytes == 12 * GIB
    assert report.usable_cpu_count == 6
    assert report.accelerator == "none"
    assert report.detection_only is True


def test_guided_resource_mappings_are_applied_and_gpu_request_is_not_proof(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _isolate_detection(monkeypatch)
    report = detect_environment(
        environ={},
        system="Linux",
        machine="x86_64",
        cpu_count=32,
        total_ram_bytes=64 * GIB,
        configured_cpu_resources={"threads": 7, "mem_mb": 2048},
        configured_gpu_resources={"gres": "gpu:1"},
    )
    assert report.usable_cpu_count == 7
    assert report.usable_ram_bytes == 2 * GIB
    assert report.configured_gpu_requested is True
    assert report.gpu_validated is False
    assert any(
        "request" in blocker and "not validated" in blocker for blocker in report.blockers
    )


def test_cgroup_affinity_and_cpu_quota_are_reported(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    values = {
        "/sys/fs/cgroup/memory.max": str(8 * GIB),
        "/sys/fs/cgroup/cpuset.cpus.effective": "0-3,8-9",
        "/sys/fs/cgroup/cpu.max": "250000 100000",
    }
    monkeypatch.setattr(environment_module, "_read_text", lambda path: values.get(str(path)))
    monkeypatch.setattr(environment_module.os, "sched_getaffinity", lambda _pid: set(range(12)))
    monkeypatch.setattr(
        environment_module.resource,
        "getrlimit",
        lambda _key: (
            environment_module.resource.RLIM_INFINITY,
            environment_module.resource.RLIM_INFINITY,
        ),
    )
    limits = environment_module._resource_limits(
        {},
        total_ram_bytes=64 * GIB,
        configured_memory_bytes=None,
        configured_cpu_count=None,
    )
    by_source = {item.source: item.value for item in limits}
    assert by_source["cgroup_v2"] == 8 * GIB
    assert by_source["process_affinity"] == 12
    assert by_source["cgroup_v2_cpuset"] == 6
    assert by_source["cgroup_v2_quota"] == 2


def test_apple_silicon_is_visible_but_remains_local_cpu(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _isolate_detection(monkeypatch)
    report = detect_environment(
        environ={}, system="Darwin", machine="arm64", cpu_count=10, total_ram_bytes=32 * GIB
    )
    assert report.apple_silicon is True
    assert report.accelerator == "apple_mps"
    assert report.compute_profile is ComputeProfile.LOCAL_CPU
    assert report.gpu_validated is False


def test_scheduler_on_path_without_active_job_is_explicitly_ambiguous(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _isolate_detection(monkeypatch)
    monkeypatch.setattr(
        environment_module,
        "_detect_scheduler",
        lambda *_args: SchedulerDetection(
            kind="slurm", detected_by=("path",), active_job=False, commands=("sbatch",)
        ),
    )
    report = detect_environment(
        environ={}, system="Linux", machine="x86_64", cpu_count=8, total_ram_bytes=16 * GIB
    )
    assert report.compute_profile is ComputeProfile.SCHEDULER
    assert report.profile_candidates == (
        ComputeProfile.SCHEDULER,
        ComputeProfile.LOCAL_CPU,
    )
    assert report.requires_confirmation is True
    assert "login/submit node" in str(report.ambiguity_reason)
    assert report.blockers


def test_active_scheduler_job_is_not_the_login_node_ambiguity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _isolate_detection(monkeypatch)
    monkeypatch.setattr(
        environment_module,
        "_detect_scheduler",
        lambda *_args: SchedulerDetection(
            kind="slurm",
            detected_by=("environment", "path"),
            active_job=True,
            job_id="123",
        ),
    )
    report = detect_environment(
        environ={"SLURM_JOB_ID": "123"},
        system="Linux",
        machine="x86_64",
        cpu_count=8,
        total_ram_bytes=16 * GIB,
    )
    assert report.compute_profile is ComputeProfile.SCHEDULER
    assert report.requires_confirmation is False


def test_cuda_candidate_is_not_local_gpu_when_runtime_validation_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _isolate_detection(monkeypatch)
    monkeypatch.setattr(
        environment_module,
        "_cuda_candidate",
        lambda *_args: (True, ("nvidia-smi is on PATH",)),
    )

    def fail() -> GPUValidation:
        raise GPUValidationError("fixture CUDA validation failed; CPU substitution is disabled")

    report = detect_environment(
        environ={},
        system="Linux",
        machine="x86_64",
        cpu_count=8,
        total_ram_bytes=16 * GIB,
        gpu_validator=fail,
    )
    assert report.cuda_candidate is True
    assert report.compute_profile is ComputeProfile.LOCAL_CPU
    assert report.gpu_validated is False
    assert "CPU substitution is disabled" in str(report.gpu_validation_error)


def test_validated_cuda_selects_local_gpu(monkeypatch: pytest.MonkeyPatch) -> None:
    _isolate_detection(monkeypatch)
    monkeypatch.setattr(
        environment_module,
        "_cuda_candidate",
        lambda *_args: (True, ("NVIDIA device nodes are present",)),
    )
    report = detect_environment(
        environ={"CUDA_VISIBLE_DEVICES": "0"},
        system="Linux",
        machine="x86_64",
        cpu_count=8,
        total_ram_bytes=64 * GIB,
        gpu_validator=_validated_gpu,
    )
    assert report.compute_profile is ComputeProfile.LOCAL_GPU
    assert report.accelerator == "cuda"
    assert report.gpu_validated is True
    assert report.gpu_validation is not None
    assert report.gpu_validation.cpu_fallback_allowed is False


def test_real_gpu_validation_executes_allocation_kernel_and_sync() -> None:
    synchronized = {"value": False}

    class Runtime:
        @staticmethod
        def getDeviceCount() -> int:
            return 1

        @staticmethod
        def getDeviceProperties(_index: int):
            return {
                "name": b"Fake A100",
                "major": 8,
                "minor": 0,
                "totalGlobalMem": 40 * GIB,
            }

        @staticmethod
        def runtimeGetVersion() -> int:
            return 12040

        @staticmethod
        def driverGetVersion() -> int:
            return 12060

    class Device:
        compute_capability = (8, 0)

        def __init__(self, _index: int) -> None:
            pass

    class Stream:
        null = SimpleNamespace(synchronize=lambda: synchronized.__setitem__("value", True))

    fake_cupy = SimpleNamespace(
        __version__="13.0.0",
        float32=np.float32,
        asarray=np.asarray,
        asnumpy=np.asarray,
        cuda=SimpleNamespace(
            is_available=lambda: True,
            runtime=Runtime,
            Device=Device,
            Stream=Stream,
        ),
    )
    validation = validate_real_gpu(cupy_module=fake_cupy, environ={"CUDA_VISIBLE_DEVICES": "0"})
    assert validation.status == "validated"
    assert validation.device_name == "Fake A100"
    assert synchronized["value"] is True


def test_real_gpu_validation_fails_when_visibility_is_disabled() -> None:
    with pytest.raises(GPUValidationError, match="CPU substitution is disabled"):
        validate_real_gpu(cupy_module=object(), environ={"CUDA_VISIBLE_DEVICES": "-1"})
    with pytest.raises(GPUValidationError, match="CPU substitution is disabled"):
        require_real_gpu(cupy_module=object(), environ={"CUDA_VISIBLE_DEVICES": "-1"})


def _request(**updates) -> RoutingRequest:
    values = {
        "level": "L2",
        "parent_scope": "immune",
        "organism": "human",
        "tissue": "blood",
        "developmental_stage": "adult",
        "health_context": HealthContext.HEALTHY,
        "technology": Technology.DISSOCIATED,
        "compute_profile": ComputeProfile.LOCAL_CPU,
        "usable_ram_bytes": 32 * GIB,
    }
    values.update(updates)
    return RoutingRequest.model_validate(values)


def _assessment(decision, method: CrosscheckMethod):
    return next(item for item in decision.assessments if item.method is method)


def test_existing_prediction_is_primary_and_every_method_has_a_skip_reason() -> None:
    request = _request(
        existing_predictions=(
            MethodResource(
                name="atlas_pred",
                organisms=("human",),
                tissues=("blood",),
                levels=("L2",),
            ),
        )
    )
    decision = route_crosscheck(request)
    assert decision.primary_method is CrosscheckMethod.EXISTING_PREDICTION
    assert decision.primary_resource == "atlas_pred"
    assert {item.method for item in decision.assessments} == set(ALL_METHODS)
    assert all(item.reasons for item in decision.assessments)
    assert sum(item.selected for item in decision.assessments) == 1


def test_labeled_reference_routes_to_cpu_and_explains_gpu_skip() -> None:
    decision = route_crosscheck(
        _request(
            labeled_references=(
                MethodResource(
                    name="local-blood-reference",
                    organisms=("human",),
                    tissues=("blood",),
                    health_contexts=(HealthContext.HEALTHY,),
                ),
            )
        )
    )
    assert decision.primary_method is CrosscheckMethod.REFERENCE_KNN_CPU
    gpu = _assessment(decision, CrosscheckMethod.REFERENCE_KNN_GPU)
    assert gpu.executable is False
    assert any("validated real CUDA GPU" in reason for reason in gpu.reasons)


def test_scarches_requires_validated_gpu_and_then_has_priority() -> None:
    model = MethodResource(
        name="blood-scanvi",
        organisms=("human",),
        tissues=("blood",),
        health_contexts=(HealthContext.HEALTHY,),
        requires_gpu=True,
    )
    with pytest.raises(ValueError, match="validated real-GPU"):
        _request(compute_profile=ComputeProfile.LOCAL_GPU, reference_models=(model,))

    decision = route_crosscheck(
        _request(
            compute_profile=ComputeProfile.LOCAL_GPU,
            gpu_validated=True,
            reference_models=(model,),
        )
    )
    assert decision.primary_method is CrosscheckMethod.REFERENCE_MODEL_SCARCHES


def test_tumor_context_rejects_healthy_reference_and_selects_tumor_reference() -> None:
    healthy = MethodResource(
        name="healthy-atlas",
        organisms=("human",),
        health_contexts=(HealthContext.HEALTHY,),
    )
    tumor = MethodResource(
        name="tumor-atlas",
        organisms=("human",),
        health_contexts=(HealthContext.TUMOR,),
    )
    decision = route_crosscheck(
        _request(
            health_context=HealthContext.TUMOR,
            labeled_references=(healthy, tumor),
        )
    )
    assert decision.primary_method is CrosscheckMethod.REFERENCE_KNN_CPU
    assert decision.primary_resource == "tumor-atlas"


def test_targeted_panel_rejects_foundation_model_and_requires_panel_overlap() -> None:
    low_overlap = MethodResource(
        name="low-overlap-ref",
        gene_overlap_count=12,
        minimum_gene_overlap=30,
        technologies=(Technology.TARGETED_PANEL,),
    )
    callable_reference = MethodResource(
        name="panel-ref",
        gene_overlap_count=40,
        minimum_gene_overlap=30,
        technologies=(Technology.TARGETED_PANEL,),
    )
    foundation = MethodResource(name="scimilarity-v1", organisms=("human",))
    decision = route_crosscheck(
        _request(
            technology=Technology.TARGETED_PANEL,
            labeled_references=(low_overlap, callable_reference),
            scimilarity_models=(foundation,),
        )
    )
    assert decision.primary_method is CrosscheckMethod.REFERENCE_KNN_CPU
    assert decision.primary_resource == "panel-ref"
    cpu_knn = _assessment(decision, CrosscheckMethod.REFERENCE_KNN_CPU)
    assert any("low-overlap-ref" in reason for reason in cpu_knn.reasons)
    scimilarity = _assessment(decision, CrosscheckMethod.SCIMILARITY)
    assert scimilarity.suitable is False
    assert any("unsuitable for targeted-panel" in reason for reason in scimilarity.reasons)


def test_memory_limit_makes_method_unavailable_with_exact_reason() -> None:
    huge_model = MethodResource(
        name="huge-model",
        estimated_memory_bytes=64 * GIB,
        organisms=("human",),
    )
    decision = route_crosscheck(_request(scimilarity_models=(huge_model,)))
    assessment = _assessment(decision, CrosscheckMethod.SCIMILARITY)
    assert assessment.executable is False
    assert any(
        "only" in reason and "bytes are usable" in reason for reason in assessment.reasons
    )


def test_scheduler_ambiguity_blocks_routing_instead_of_assuming_cpu() -> None:
    reference = MethodResource(name="local-ref")
    decision = route_crosscheck(
        _request(
            compute_profile=ComputeProfile.SCHEDULER,
            scheduler_ambiguity=True,
            labeled_references=(reference,),
        )
    )
    assert decision.status == "blocked"
    assert decision.primary_method is None
    assert "scheduler/login-node ambiguity" in decision.blockers[0]


def test_no_resources_is_honestly_unavailable() -> None:
    decision = route_crosscheck(_request())
    assert decision.status == "unavailable"
    assert decision.primary_method is None
    assert len(decision.assessments) == len(ALL_METHODS)
    assert all(item.remediation for item in decision.assessments)


def test_census_requires_network_authorization_and_validated_cache() -> None:
    census = MethodResource(
        name="census-blood-v1",
        network_required=True,
        checksum_validated=False,
        organisms=("human",),
    )
    denied = route_crosscheck(_request(census_references=(census,)))
    assessment = _assessment(denied, CrosscheckMethod.CENSUS_REFERENCE_KNN)
    assert assessment.executable is False
    assert any("not authorized" in reason for reason in assessment.reasons)
    assert any("checksum" in reason for reason in assessment.reasons)


def test_hierarchy_routing_is_parent_scoped_and_rejects_duplicate_keys() -> None:
    reference = MethodResource(name="local-ref")
    requests = (
        _request(level="L2", parent_scope="immune", labeled_references=(reference,)),
        _request(level="L2", parent_scope="stromal", labeled_references=(reference,)),
    )
    decisions = route_hierarchy(requests)
    assert [item.parent_scope for item in decisions] == ["immune", "stromal"]
    with pytest.raises(ValueError, match="duplicate"):
        route_hierarchy((requests[0], requests[0]))
