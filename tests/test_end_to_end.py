from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import SKILL_ROOT

from cell_curator.config import CellCuratorConfig, load_config
from cell_curator.contracts import ContractError, sha256
from cell_curator.guided import approve_assumptions, construct_plan
from cell_curator.review import (
    build_critic_packet,
    build_evidence_cards,
    build_html_report,
    import_critic_findings,
    import_reconciliation,
    resolve_durable_mapping_values,
)
from cell_curator.writeback import validate_writeback, write_back


def _synthetic_input(path: Path) -> None:
    rng = np.random.default_rng(7)
    counts = np.zeros((60, 6), dtype=np.float32)
    counts[:30, :2] = rng.poisson(8, size=(30, 2)) + 1
    counts[30:, 2:4] = rng.poisson(8, size=(30, 2)) + 1
    obs = pd.DataFrame(
        {
            "canonical_cluster": ["0"] * 30 + ["1"] * 30,
            "lineage": ["immune"] * 60,
        },
        index=[f"cell-{index:03d}" for index in range(60)],
    )
    var = pd.DataFrame(
        {"gene_symbol": ["A1", "A2", "B1", "B2", "N1", "N2"]},
        index=["A1", "A2", "B1", "B2", "N1", "N2"],
    )
    adata = ad.AnnData(X=np.log1p(counts), obs=obs, var=var)
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(30, 3)), rng.normal(3, 0.2, size=(30, 3))]
    ).astype(np.float32)
    adata.write_h5ad(path)


def _config(tmp_path: Path, source: Path, markers: Path) -> Path:
    value = yaml.safe_load((SKILL_ROOT / "assets" / "config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "synthetic-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "annotation",
        }
    )
    value["input"]["path"] = str(source)
    value["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["adaptive"]["eligibility"].update(
        {"min_cells": 5, "min_donors": None, "min_captures": None}
    )
    value["adaptive"]["parent_local"]["enabled"] = False
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "synthetic-v1"}
    ]
    value["guidance"] = {
        "interactive": False,
        "preauthorized": True,
        "reviewer": "synthetic-fixture-reviewer",
        "pause_after_each_level": False,
        "assumptions": ["Synthetic marker programs and context are pre-authorized."],
        "hierarchy_levels": ["L1"],
        "context": {
            "organism": "human",
            "tissue": "synthetic tissue",
            "developmental_stage": "adult",
            "condition": "healthy",
            "technology": "dissociated",
            "experimental_context": "synthetic executable integration fixture",
            "composition_covariates": ["lineage"],
            "vocabulary_contract": "local synthetic marker programs",
            "reference_contract": "local_only",
            "visualization_key": "X_pca_harmony",
            "writeback_prefix": "cell_curator",
        },
    }
    value["evidence"].update(
        {
            "min_positive_coverage": 1.0,
            "max_negative_violation": 0.0,
            "min_confidence": 0.6,
            "min_margin": 0.2,
        }
    )
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(value, sort_keys=False))
    return path


def _completed_synthetic_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, CellCuratorConfig, Path]:
    """A corpus run with consensus labels, finalized, reviewed and critic-attested."""

    from consensus_fixture import write_consensus_run
    from test_annotation_pipeline import annotation_fixture

    from cell_curator.annotation import collect_level, finalize_annotation
    from cell_curator.consensus_import import import_consensus_level

    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch)
    collect_level(config_path, "L1")
    write_consensus_run(
        config_path, root, {cluster: values[0] for cluster, values in expected.items()}
    )
    import_consensus_level(config_path, "L1")
    finalize_annotation(config_path)
    config = load_config(config_path)
    build_evidence_cards(config, root)
    build_critic_packet(config, root)
    build_html_report(config, root, strict=True)
    _record_critics(config, root, tmp_path)
    return Path(config.input.path), config, root


def _rebuild_strict_review(config: CellCuratorConfig, root: Path) -> None:
    build_evidence_cards(config, root)
    build_critic_packet(config, root)
    build_html_report(config, root, strict=True, level="L1")
    build_html_report(config, root, strict=True)


def _record_critics(config: CellCuratorConfig, root: Path, tmp_path: Path) -> None:
    """Every configured pass attests to finding nothing, and a reviewer signs it off."""

    passes = config.review.critic_passes
    findings = tmp_path / "critic-findings.tsv"
    pd.DataFrame(
        {
            "finding_id": [f"none-{name}" for name in passes],
            "target_kind": "run",
            "target_id": config.run.run_id,
            "critic_pass": passes,
            "severity": "NONE",
            "claim": "no issue found",
            "evidence_path": "review/critic_packet.json",
            "proposed_resolution": "none",
        }
    ).to_csv(findings, sep="\t", index=False)
    import_critic_findings(root, findings)
    reconciliation = tmp_path / "critic-reconciliation.tsv"
    pd.DataFrame(
        {
            "finding_id": [f"none-{name}" for name in passes],
            "status": "DISMISSED",
            "resolution": "attestation accepted",
            "reviewer": "fixture-reviewer",
        }
    ).to_csv(reconciliation, sep="\t", index=False)
    import_reconciliation(root, reconciliation)
    build_html_report(config, root, strict=True, level="L1")
    build_html_report(config, root, strict=True)


def _typed(prefix_from: dict[str, object]):
    """A person typing the diff prefix shown for ``prefix_from``."""

    return lambda _preview: str(prefix_from["approval_token"])


def _confirm(preview: dict[str, object]) -> str:
    return str(preview["approval_token"])


def test_durable_mapping_resolution_uses_deepest_parent_scoped_review() -> None:
    mapping = pd.DataFrame(
        {
            "barcode": ["a", "b"],
            "canonical_parent": ["0", "1"],
            "leaf_id": ["global::0::candidate::child-a", "global::1"],
        }
    )
    labels = pd.DataFrame(
        [
            {
                "level": "L1",
                "parent_scope": "ROOT",
                "cluster_id": "0",
                "label": "Parent A",
                "cell_state": "normal",
            },
            {
                "level": "L1",
                "parent_scope": "ROOT",
                "cluster_id": "1",
                "label": "Parent B",
                "cell_state": "normal",
            },
            {
                "level": "L2",
                "parent_scope": "Parent A",
                "cluster_id": "child-a",
                "label": "Reviewed child",
                "cell_state": "cycling",
            },
        ]
    )
    l1 = resolve_durable_mapping_values(mapping, labels, max_level="L1")
    full = resolve_durable_mapping_values(mapping, labels)
    assert list(l1["cell_type"]) == ["Parent A", "Parent B"]
    assert list(full["cell_type"]) == ["Reviewed child", "Parent B"]
    assert list(full["cell_state"]) == ["cycling", "normal"]


def test_writeback_confirmation_binds_the_previewed_diff(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, root = _completed_synthetic_run(tmp_path, monkeypatch)
    output = tmp_path / "bound-output.h5ad"
    preview = validate_writeback(config, root, output_path=output)
    # A prefix typed for a different target is a different diff and refuses.
    other = validate_writeback(config, root, output_path=tmp_path / "other.h5ad")
    assert other["approval_token"] != preview["approval_token"]
    with pytest.raises(ContractError, match="does not match"):
        write_back(config, root, output_path=output, confirm=lambda _p: "0" * 12)
    assert not output.exists()
    write_back(config, root, output_path=output, confirm=_typed(preview))
    assert output.is_file()


def test_writeback_refuses_labels_older_than_a_later_approved_assumption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from cell_curator.api import record_annotation_assumption

    _, config, root = _completed_synthetic_run(tmp_path, monkeypatch)
    config_path = tmp_path / "config.yaml"
    output = tmp_path / "after-new-assumption.h5ad"
    validate_writeback(config, root, output_path=output)
    record_annotation_assumption(
        config_path,
        scope="L1::ROOT",
        statement="Cluster 1 is a known contaminant and must stay Unknown.",
    )
    approve_assumptions(config_path, reviewer="fixture-reviewer", preauthorized=True)
    with pytest.raises(ContractError, match="stale assumptions"):
        validate_writeback(config, root, output_path=output)


def test_a_finer_level_abstention_keeps_the_supported_parent_label() -> None:
    mapping = pd.DataFrame({"barcode": ["b1"], "leaf_id": ["c1"], "canonical_parent": ["c1"]})
    rows = [
        ("L1", "ROOT", "c1", "T cell"),
        ("L2", "T cell", "c1", "Unknown"),
    ]
    labels = pd.DataFrame(
        [
            {
                "level": lv,
                "parent_scope": ps,
                "cluster_id": cid,
                "label": lab,
                "cell_state": "normal",
            }
            for lv, ps, cid, lab in rows
        ]
    )
    resolved = resolve_durable_mapping_values(mapping, labels)
    assert list(resolved["cell_type"]) == ["T cell"]


def test_ontology_backoff_is_a_plan_assumption_only_when_an_ontology_is_configured(
    tmp_path: Path,
) -> None:
    from cell_curator.guided import load_assumptions
    from cell_curator.ontology import BACKOFF_ASSUMPTION

    source = tmp_path / "input.h5ad"
    _synthetic_input(source)
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "myeloid cell": {"pos": ["A1"], "cl_id": "CL:0000763"},
                "monocyte": {"pos": ["B1"], "cl_id": "CL:0000576"},
            }
        )
    )
    obo = tmp_path / "cl.obo"
    obo.write_text(
        "format-version: 1.2\n\n[Term]\nid: CL:0000000\nname: cell\n\n"
        "[Term]\nid: CL:0000763\nname: myeloid cell\nis_a: CL:0000000\n\n"
        "[Term]\nid: CL:0000576\nname: monocyte\nis_a: CL:0000763 ! myeloid cell\n"
    )
    for run_id, with_ontology in (("cl-free", False), ("with-cl", True)):
        config_path = _config(tmp_path, source, markers)
        value = yaml.safe_load(config_path.read_text())
        value["run"]["run_id"] = run_id
        if with_ontology:
            value["knowledge"]["providers"].append(
                {"kind": "ontology", "path": str(obo), "version": "fixture"}
            )
        config_path.write_text(yaml.safe_dump(value, sort_keys=False))
        construct_plan(config_path)
        rows = load_assumptions(tmp_path / "results" / run_id)
        statements = {row["statement"] for row in rows if row["scope"] == "global"}
        assert (BACKOFF_ASSUMPTION in statements) is with_ontology
        # A label and its own ancestor as siblings is put to the reviewer.
        assert (
            any(
                row["statement"].startswith(
                    "myeloid cell competes with its Cell Ontology descendants monocyte"
                )
                for row in rows
            )
            is with_ontology
        )


def test_a_cli_failure_is_one_json_object_with_its_own_exit_code(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    from cell_curator import cli

    def refuse(_args: object) -> None:
        raise ContractError("Gate A blocks biological labels")

    monkeypatch.setattr(cli, "dispatch", refuse)
    monkeypatch.setattr(
        "sys.argv", ["cell-curator", "progress", "status", "--config", "c.yaml"]
    )
    with pytest.raises(SystemExit) as stopped:
        cli.main()
    assert stopped.value.code == 2
    error = json.loads(capsys.readouterr().err.strip().splitlines()[-1])
    assert error == {
        "error": "ContractError",
        "message": "Gate A blocks biological labels",
        "exit_code": 2,
    }


def test_a_full_run_and_write_back_open_no_socket(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Offline is the default, and the whole path works with the network gone."""

    import socket

    def no_network(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("cell-curator opened a network connection")

    monkeypatch.setattr(socket.socket, "connect", no_network)
    monkeypatch.setattr(socket.socket, "connect_ex", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    _, config, root = _completed_synthetic_run(tmp_path, monkeypatch)
    assert config.knowledge.network == "offline"
    output = tmp_path / "offline.h5ad"
    write_back(config, root, output_path=output, confirm=_confirm)
    assert set(ad.read_h5ad(output).obs["cell_curator_type"].astype(str)).difference(
        {"Unknown"}
    )


def test_a_write_back_interrupted_before_its_receipt_completes_on_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A crash after the output lands but before the receipt does must not wedge the run."""

    import os

    from cell_curator import writeback

    _, config, root = _completed_synthetic_run(tmp_path, monkeypatch)
    output = tmp_path / "interrupted.h5ad"
    real_replace = os.replace

    def crash_on_receipt(source: object, target: object) -> None:
        if str(target).endswith("receipt.json"):
            raise OSError("power lost")
        real_replace(source, target)

    monkeypatch.setattr(writeback.os, "replace", crash_on_receipt)
    with pytest.raises(OSError, match="power lost"):
        write_back(config, root, output_path=output, confirm=_confirm)
    assert output.is_file() and not (root / "writeback/receipt.json").exists()

    monkeypatch.setattr(writeback.os, "replace", real_replace)
    assert write_back(config, root, output_path=output, confirm=None) == output
    receipt = json.loads((root / "writeback/receipt.json").read_text())
    assert receipt["output_sha256"] == sha256(output)
    assert not (root / "writeback/receipt.json.pending").exists()


def test_repeating_a_completed_write_back_returns_the_same_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, root = _completed_synthetic_run(tmp_path, monkeypatch)
    output = tmp_path / "once.h5ad"
    write_back(config, root, output_path=output, confirm=_confirm)
    first = sha256(output)
    assert write_back(config, root, output_path=output, confirm=None) == output
    assert sha256(output) == first


def test_a_filesystem_that_cannot_fsync_a_directory_still_completes_write_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import os
    import stat

    from cell_curator import writeback

    _, config, root = _completed_synthetic_run(tmp_path, monkeypatch)
    real_fsync = os.fsync

    def no_directory_fsync(fd: int) -> None:
        if stat.S_ISDIR(os.fstat(fd).st_mode):
            raise OSError(22, "Invalid argument")
        real_fsync(fd)

    monkeypatch.setattr(writeback.os, "fsync", no_directory_fsync)
    output = tmp_path / "fuse.h5ad"
    write_back(config, root, output_path=output, confirm=_confirm)
    assert json.loads((root / "writeback/receipt.json").read_text())["output_sha256"] == sha256(
        output
    )
