from __future__ import annotations

import hashlib
import json
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cell_curator.crosschecks as crosschecks_module
import cell_curator.environment as environment_module
from cell_curator.contracts import ContractError
from cell_curator.crosschecks import (
    celltypist_hypothesis,
    discover_celltypist_models,
    fetch_census_reference,
    map_scarches_query,
    train_scanvi_reference,
)
from cell_curator.environment import GPUValidationError


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _install_fake_celltypist(
    monkeypatch: pytest.MonkeyPatch,
    cache_dir: Path,
) -> dict[str, Any]:
    calls: dict[str, Any] = {
        "annotate_models": [],
        "get_all_models": 0,
        "download_models": 0,
    }

    def get_all_models() -> list[str]:
        calls["get_all_models"] += 1
        pytest.fail("read-only discovery invoked the potentially remote model catalog")

    def download_models(*_args: Any, **_kwargs: Any) -> None:
        calls["download_models"] += 1
        pytest.fail("CellTypist automatic model download was invoked")

    models = SimpleNamespace(
        models_path=str(cache_dir),
        get_all_models=get_all_models,
        download_models=download_models,
    )

    def annotate(
        adata: ad.AnnData,
        *,
        model: str,
        majority_voting: bool,
        over_clustering: str | None,
    ) -> SimpleNamespace:
        resolved = Path(model)
        assert resolved.is_absolute()
        assert resolved.is_file()
        calls["annotate_models"].append(model)
        assert majority_voting is False
        assert over_clustering is None
        frame = pd.DataFrame(
            {
                "predicted_labels": ["T cell", "B cell"],
                "conf_score": [0.9, 0.75],
            },
            index=adata.obs_names,
        )
        return SimpleNamespace(predicted_labels=frame)

    fake = SimpleNamespace(models=models, annotate=annotate, __version__="1.7.fake")
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setitem(sys.modules, "celltypist", fake)
    return calls


def _celltypist_input() -> ad.AnnData:
    return ad.AnnData(
        X=np.ones((2, 2), dtype=np.float32),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )


def test_celltypist_discovery_reads_only_cache_and_marks_remote_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_dir = tmp_path / "celltypist-models"
    cache_dir.mkdir()
    cached = cache_dir / "cached.pkl"
    cached.write_bytes(b"cached CellTypist model")
    (cache_dir / "models.json").write_text(
        json.dumps(
            {
                "models": [
                    {"filename": "cached.pkl", "description": "local"},
                    {"filename": "remote.pkl", "description": "catalog only"},
                ]
            }
        )
    )
    calls = _install_fake_celltypist(monkeypatch, cache_dir)

    status = discover_celltypist_models(output_dir=tmp_path / "status")

    assert status["status"] == "complete"
    assert status["models"] == ["cached.pkl", "remote.pkl"]
    entries = {entry["identifier"]: entry for entry in status["model_entries"]}
    assert entries["cached.pkl"] == {
        "identifier": "cached.pkl",
        "cached": True,
        "availability": "cached",
        "path": str(cached.resolve()),
        "sha256": _digest(cached),
    }
    assert entries["remote.pkl"] == {
        "identifier": "remote.pkl",
        "cached": False,
        "availability": "remote_metadata_only",
        "path": "",
        "sha256": "",
    }
    assert status["discovery_mode"] == "local_cache_only"
    assert status["remote_catalog_consulted"] is False
    assert status["network_accessed"] is False
    assert status["downloaded"] is False
    assert status["celltypist_version"] == "1.7.fake"
    assert calls == {
        "annotate_models": [],
        "get_all_models": 0,
        "download_models": 0,
    }
    assert json.loads((tmp_path / "status" / "provider_status.json").read_text()) == status


@pytest.mark.parametrize("use_identifier", [False, True])
def test_celltypist_prediction_requires_and_hashes_local_model_without_download(
    use_identifier: bool,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_dir = tmp_path / "celltypist-models"
    cache_dir.mkdir()
    cached = cache_dir / "cached.pkl"
    cached.write_bytes(b"immutable cached CellTypist model")
    calls = _install_fake_celltypist(monkeypatch, cache_dir)
    requested = cached.name if use_identifier else str(cached)
    output_dir = tmp_path / "prediction"

    cells, status = celltypist_hypothesis(
        _celltypist_input(),
        model=requested,
        output_dir=output_dir,
    )

    expected_digest = _digest(cached)
    assert cells["predicted_label"].tolist() == ["T cell", "B cell"]
    assert cells["uncertainty"].tolist() == pytest.approx([0.1, 0.25])
    assert set(cells["provider"]) == {"celltypist"}
    assert set(cells["role"]) == {"hypothesis_not_verdict"}
    assert status["status"] == "complete"
    assert status["automated_output_role"] == "hypothesis_not_verdict"
    assert status["model_identifier"] == requested
    assert status["resolved_model_name"] == "cached.pkl"
    assert status["model_path"] == str(cached.resolve())
    assert status["model_source"] == (
        "celltypist_cache" if use_identifier else "explicit_local_path"
    )
    assert status["model_sha256"] == expected_digest
    assert status["model_version"] == f"sha256:{expected_digest}"
    assert status["celltypist_version"] == "1.7.fake"
    assert status["network_accessed"] is False
    assert status["automatic_downloads_allowed"] is False
    assert calls["annotate_models"] == [str(cached.resolve())]
    assert calls["get_all_models"] == 0
    assert calls["download_models"] == 0
    durable_cells = pd.read_csv(status["artifacts"]["cells"], sep="\t")
    assert set(durable_cells["role"]) == {"hypothesis_not_verdict"}
    assert json.loads((output_dir / "provider_status.json").read_text()) == status


def test_celltypist_missing_cached_model_is_actionable_and_never_downloads(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_dir = tmp_path / "celltypist-models"
    cache_dir.mkdir()
    calls = _install_fake_celltypist(monkeypatch, cache_dir)
    output_dir = tmp_path / "prediction"

    cells, status = celltypist_hypothesis(
        _celltypist_input(),
        model="remote-only.pkl",
        output_dir=output_dir,
    )

    assert cells.empty
    assert status["status"] == "unavailable"
    assert (
        status["reason"]
        == "CellTypist model 'remote-only.pkl' is not present as a local file or cached model"
    )
    assert status["remediation"] == (
        "Cache the model explicitly during an authorized preparation step, then pass its local "
        "file path or cached CellTypist identifier; automatic downloads are disabled."
    )
    assert status["model_identifier"] == "remote-only.pkl"
    assert status["model_path"] == ""
    assert status["model_sha256"] == ""
    assert status["network_accessed"] is False
    assert status["automatic_downloads_allowed"] is False
    assert status["automated_output_role"] == "hypothesis_not_verdict"
    assert calls == {
        "annotate_models": [],
        "get_all_models": 0,
        "download_models": 0,
    }
    assert json.loads((output_dir / "provider_status.json").read_text()) == status


@pytest.mark.parametrize("operation", ["train-scanvi", "map-scarches"])
def test_reference_models_fail_closed_before_dependency_or_cpu_routes(
    operation: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    message = "no validated real CUDA GPU; CPU substitution is prohibited"

    def reject_gpu() -> None:
        raise GPUValidationError(message)

    monkeypatch.setattr(environment_module, "require_real_gpu", reject_gpu)
    monkeypatch.setattr(
        crosschecks_module.importlib.util,
        "find_spec",
        lambda _name: pytest.fail("dependency discovery ran before the real-GPU gate"),
    )
    monkeypatch.setattr(
        crosschecks_module,
        "cKDTree",
        lambda *_args, **_kwargs: pytest.fail("a CPU reference mapper was substituted"),
    )
    output_dir = tmp_path / operation

    with pytest.raises(GPUValidationError) as error:
        if operation == "train-scanvi":
            train_scanvi_reference(None, labels_key="cell_type", output_dir=output_dir)
        else:
            map_scarches_query(None, model_dir=tmp_path / "model", output_dir=output_dir)

    assert str(error.value) == message
    assert not output_dir.exists()


def test_reference_models_report_exact_missing_extra_after_gpu_validation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(environment_module, "require_real_gpu", lambda: None)
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: None)

    training = train_scanvi_reference(
        None,
        labels_key="cell_type",
        output_dir=tmp_path / "training",
    )
    cells, mapping = map_scarches_query(
        None,
        model_dir=tmp_path / "model",
        output_dir=tmp_path / "mapping",
    )

    assert cells.empty
    for status, provider in (
        (training, "scanvi_training"),
        (mapping, "scarches_mapping"),
    ):
        assert status["provider"] == provider
        assert status["status"] == "unavailable"
        assert status["reason"] == "optional Python distribution 'scvi' is not installed"
        assert (
            status["remediation"]
            == "Install cell-curator[reference-models] in the execution environment and retry."
        )


def test_cached_census_receipt_records_version_timeout_filter_and_checksum(
    tmp_path: Path,
) -> None:
    cached = tmp_path / "reference.h5ad"
    cached.write_bytes(b"checksum-pinned cached Census fixture")
    expected = _digest(cached)

    status = fetch_census_reference(
        organism="Homo sapiens",
        tissue="lung",
        output_path=cached,
        census_version="2025-01-30",
        max_cells=123,
        label_columns=("cell_type", "subclass"),
        expected_sha256=expected,
        timeout_seconds=12.5,
    )

    assert status["status"] == "complete"
    assert status["cache_hit"] is True
    assert status["checksum_validated"] is True
    assert status["expected_sha256"] == expected
    assert status["artifacts"] == {"reference": str(cached), "sha256": expected}
    assert status["census_version"] == "2025-01-30"
    assert status["timeout_seconds"] == 12.5
    assert status["max_cells"] == 123
    assert status["label_columns"] == ["cell_type", "subclass"]
    assert status["tissue_filter"] == "tissue_general == 'lung'"
    assert json.loads((tmp_path / "provider_status.json").read_text()) == status


def test_cached_census_checksum_mismatch_is_explicit_and_preserves_cache(
    tmp_path: Path,
) -> None:
    cached = tmp_path / "reference.h5ad"
    payload = b"untrusted cached Census fixture"
    cached.write_bytes(payload)

    status = fetch_census_reference(
        organism="Homo sapiens",
        tissue="lung",
        output_path=cached,
        census_version="2025-01-30",
        expected_sha256="0" * 64,
        timeout_seconds=9.0,
    )

    assert status["status"] == "failed"
    assert status["cache_hit"] is True
    assert status["checksum_validated"] is False
    assert (
        status["reason"]
        == "cached Census reference checksum does not match the declared SHA-256"
    )
    assert (
        status["remediation"]
        == "Remove or replace the untrusted cache entry and retry with the pinned reference receipt."
    )
    assert status["artifacts"] == {
        "untrusted_reference": str(cached),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }
    assert cached.read_bytes() == payload
    assert status["census_version"] == "2025-01-30"
    assert status["timeout_seconds"] == 9.0


def test_uncached_census_missing_dependency_has_exact_remediation_and_receipt(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: None)
    output = tmp_path / "reference.h5ad"

    status = fetch_census_reference(
        organism="Homo sapiens",
        tissue="lung",
        output_path=output,
        census_version="2025-01-30",
        expected_sha256="a" * 64,
        timeout_seconds=7.5,
    )

    assert status["status"] == "unavailable"
    assert (
        status["reason"] == "optional Python distribution 'cellxgene-census' is not installed"
    )
    assert (
        status["remediation"]
        == "Install cell-curator[census] in the execution environment and retry."
    )
    assert status["cache_hit"] is False
    assert status["checksum_validated"] is False
    assert status["expected_sha256"] == "a" * 64
    assert status["census_version"] == "2025-01-30"
    assert status["timeout_seconds"] == 7.5
    assert not output.exists()
    assert json.loads((tmp_path / "provider_status.json").read_text()) == status


@pytest.mark.parametrize(
    ("override", "message"),
    [
        ({"census_version": "  "}, "Census census_version must be nonblank"),
        (
            {"timeout_seconds": float("nan")},
            "Census timeout_seconds must be finite and positive",
        ),
        ({"expected_sha256": "ABC"}, "Census expected_sha256 must be lowercase SHA-256 hex"),
        (
            {"label_columns": ("cell_type", " ")},
            "Census label_columns must contain at least one nonblank column name",
        ),
    ],
)
def test_census_request_contract_rejects_invalid_receipt_fields(
    override: dict[str, Any],
    message: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        crosschecks_module.importlib.util,
        "find_spec",
        lambda _name: pytest.fail("invalid request reached dependency discovery"),
    )
    arguments: dict[str, Any] = {
        "organism": "Homo sapiens",
        "tissue": "lung",
        "output_path": tmp_path / "reference.h5ad",
        "census_version": "2025-01-30",
        "timeout_seconds": 5.0,
    }
    arguments.update(override)

    with pytest.raises(ContractError) as error:
        fetch_census_reference(**arguments)
    assert str(error.value) == message


class _FakeQueue:
    def __init__(self, payload: dict[str, Any] | None = None) -> None:
        self.payload = payload
        self.closed = False

    def get(self, *, timeout: float) -> dict[str, Any]:
        assert timeout == 1.0
        if self.payload is None:
            raise AssertionError("the timeout path must not read the result queue")
        return self.payload

    def put(self, payload: dict[str, Any]) -> None:
        self.payload = payload

    def close(self) -> None:
        self.closed = True


class _FakeProcess:
    def __init__(
        self,
        context: _FakeContext,
        *,
        target: Any,
        kwargs: dict[str, Any],
    ) -> None:
        self.context = context
        self.target = target
        self.kwargs = kwargs
        self.started = False
        self.terminated = False
        self.join_calls: list[float | None] = []
        self.exitcode = context.exitcode

    def start(self) -> None:
        self.started = True
        if self.context.create_partial:
            Path(self.kwargs["destination"]).write_bytes(b"partial Census artifact")

    def join(self, timeout: float | None = None) -> None:
        self.join_calls.append(timeout)

    def is_alive(self) -> bool:
        return self.context.alive

    def terminate(self) -> None:
        self.terminated = True
        self.context.alive = False


class _FakeContext:
    def __init__(
        self,
        *,
        alive: bool,
        payload: dict[str, Any] | None,
        create_partial: bool = True,
        exitcode: int = 0,
    ) -> None:
        self.alive = alive
        self.create_partial = create_partial
        self.exitcode = exitcode
        self.queue = _FakeQueue(payload)
        self.process: _FakeProcess | None = None

    def Queue(self) -> _FakeQueue:
        return self.queue

    def Process(self, *, target: Any, kwargs: dict[str, Any]) -> _FakeProcess:
        self.process = _FakeProcess(self, target=target, kwargs=kwargs)
        return self.process


def test_census_worker_timeout_is_deterministic_atomic_and_actionable(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context = _FakeContext(alive=True, payload=None)
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        crosschecks_module.multiprocessing,
        "get_context",
        lambda _method: context,
    )
    output = tmp_path / "reference.h5ad"

    status = fetch_census_reference(
        organism="Homo sapiens",
        tissue="lung",
        output_path=output,
        census_version="2025-01-30",
        timeout_seconds=2.5,
    )

    assert status["status"] == "failed"
    assert status["reason"] == "bounded Census acquisition exceeded its 2.5-second timeout"
    assert (
        status["remediation"]
        == "Increase --timeout-seconds only after checking network access and the bounded query size."
    )
    assert status["census_version"] == "2025-01-30"
    assert status["timeout_seconds"] == 2.5
    assert status["cache_hit"] is False
    assert status["checksum_validated"] is False
    assert context.queue.closed is True
    assert context.process is not None
    assert context.process.started is True
    assert context.process.terminated is True
    assert context.process.join_calls == [2.5, None]
    assert not output.exists()
    assert not list(tmp_path.glob(".reference.h5ad.partial.*"))


def test_census_worker_failure_is_deterministic_atomic_and_actionable(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context = _FakeContext(
        alive=False,
        payload={"status": "failed", "reason": "synthetic provider failure"},
        exitcode=3,
    )
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        crosschecks_module.multiprocessing,
        "get_context",
        lambda _method: context,
    )
    output = tmp_path / "reference.h5ad"

    status = fetch_census_reference(
        organism="Homo sapiens",
        tissue="lung",
        output_path=output,
        census_version="2025-01-30",
        timeout_seconds=4.0,
    )

    assert status["status"] == "failed"
    assert status["reason"] == "bounded Census acquisition failed: synthetic provider failure"
    assert (
        status["remediation"]
        == "Check the pinned Census version, organism/tissue filter, network, and cache permissions."
    )
    assert status["census_version"] == "2025-01-30"
    assert status["timeout_seconds"] == 4.0
    assert context.queue.closed is True
    assert context.process is not None
    assert context.process.join_calls == [4.0]
    assert context.process.terminated is False
    assert not output.exists()
    assert not list(tmp_path.glob(".reference.h5ad.partial.*"))


def test_census_worker_uses_official_filter_and_excludes_blank_labels(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    observed: dict[str, Any] = {}

    @contextmanager
    def open_soma(*, census_version: str):
        observed["census_version"] = census_version
        yield object()

    def get_obs(
        _census: object,
        organism: str,
        *,
        value_filter: str,
        column_names: list[str],
    ) -> pd.DataFrame:
        observed.update(
            {
                "organism": organism,
                "value_filter": value_filter,
                "column_names": column_names,
            }
        )
        return pd.DataFrame(
            {
                "soma_joinid": [101, 102, 103],
                "cell_type": ["T cell", "", "B cell"],
                "subclass": ["T lineage", "T lineage", None],
            }
        )

    def get_anndata(
        _census: object,
        organism: str,
        *,
        obs_coords: np.ndarray,
        obs_column_names: list[str],
    ) -> ad.AnnData:
        observed["obs_coords"] = obs_coords.tolist()
        observed["obs_column_names"] = obs_column_names
        assert organism == "Homo sapiens"
        return ad.AnnData(
            X=np.ones((1, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {"cell_type": ["T cell"], "subclass": ["T lineage"]},
                index=["c0"],
            ),
            var=pd.DataFrame(index=["G1"]),
        )

    fake_census = SimpleNamespace(
        open_soma=open_soma,
        get_obs=get_obs,
        get_anndata=get_anndata,
    )
    monkeypatch.setitem(sys.modules, "cellxgene_census", fake_census)
    result_queue = _FakeQueue()
    output = tmp_path / "worker.h5ad"

    crosschecks_module._acquire_census_worker(
        result_queue,
        destination=str(output),
        organism="Homo sapiens",
        tissue="lung",
        census_version="2025-01-30",
        max_cells=10,
        label_columns=("cell_type", "subclass"),
    )

    assert result_queue.payload == {"status": "complete", "n_cells": 1}
    assert observed == {
        "census_version": "2025-01-30",
        "organism": "Homo sapiens",
        "value_filter": "tissue_general == 'lung'",
        "column_names": ["soma_joinid", "cell_type", "subclass"],
        "obs_coords": [101],
        "obs_column_names": ["cell_type", "subclass"],
    }
    assert output.is_file()
