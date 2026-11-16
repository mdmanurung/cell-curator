"""Corpus hierarchy previews preserve human review and never apply a write."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from consensus_fixture import write_consensus_run
from test_annotation_pipeline import annotation_fixture
from test_end_to_end import _record_critics

from cell_curator.annotation import collect_level, finalize_annotation
from cell_curator.config import load_config
from cell_curator.consensus_import import import_consensus_level
from cell_curator.contracts import ContractError, sha256
from cell_curator.pipeline import validate
from cell_curator.review import build_critic_packet, build_evidence_cards, build_html_report
from cell_curator.writeback import (
    _align_mapping,
    _load_durable_labels,
    _writeback_values,
    validate_writeback,
)


def _completed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
    return config, root, expected


def test_corpus_preview_has_all_identity_and_state_layers_without_writing(tmp_path, monkeypatch):
    config, root, expected = _completed(tmp_path, monkeypatch)
    source_hash = sha256(Path(config.input.path))
    target = tmp_path / "annotated.h5ad"
    validation = validate(tmp_path / "config.yaml")
    assert validation["status"] == "complete" and validation["writeback_permitted"] is False
    with pytest.raises(ContractError, match="critic pass"):
        validate_writeback(config, root, output_path=target)
    _record_critics(config, root, tmp_path)
    preview = validate_writeback(config, root, output_path=target)
    assert set(preview["output_columns"]) == {
        "cell_curator_type", "cell_curator_state", "cell_curator_leaf",
        "cell_curator_L1", "cell_curator_L4",
    }
    assert preview["critics_reconciled"] is True
    assert preview["n_cells"] == 160
    assert preview["annotation_cells_sha256"] == sha256(root / "annotation" / "cell_labels.tsv")
    diff = pd.read_csv(preview["diff_path"], sep="\t", keep_default_na=False)
    cells = pd.read_csv(root / "annotation" / "cell_labels.tsv", sep="\t", keep_default_na=False)
    first = cells[cells["cluster_id"].astype(str).eq("0")].iloc[0]
    observed = diff[diff["barcode"].eq(first["barcode"])].set_index("column")["new_value"]
    assert observed["cell_curator_L1"] == expected["0"][0]
    assert observed["cell_curator_L4"] == expected["0"][0] + " S+"
    assert observed["cell_curator_state"] == "S+"
    assert not target.exists()
    assert not (root / "writeback" / "receipt.json").exists()
    assert sha256(Path(config.input.path)) == source_hash
    report = (root / "review" / "report.html").read_text()
    assert "<svg" in report and "<img" not in report


def test_corpus_preview_rejects_changed_sidecar(tmp_path, monkeypatch):
    config, root, _ = _completed(tmp_path, monkeypatch)
    _record_critics(config, root, tmp_path)
    path = root / "annotation" / "cell_labels.tsv"
    path.write_text(path.read_text().replace("Lineage A", "Wrong label"))
    with pytest.raises(ContractError, match="stale annotation artifact"):
        validate_writeback(config, root, output_path=tmp_path / "annotated.h5ad")


def test_corpus_values_must_match_durable_labels_even_with_rehashed_sidecar(tmp_path, monkeypatch):
    config, root, _ = _completed(tmp_path, monkeypatch)
    path = root / "annotation" / "cell_labels.tsv"
    cells = pd.read_csv(path, sep="\t", keep_default_na=False)
    cells.loc[0, "L1"] = "Wrong label"
    cells.to_csv(path, sep="\t", index=False)
    manifest_path = root / "annotation" / "annotation.manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["artifacts"][str(path)] = sha256(path)
    manifest_path.write_text(json.dumps(manifest))
    mapping = pd.read_csv(root / "mapping" / "parent_child_mapping.tsv", sep="\t", keep_default_na=False)
    labels, _ = _load_durable_labels(root)
    mapping = _align_mapping(mapping, pd.Index(cells["barcode"].astype(str)))
    with pytest.raises(ContractError, match="L1 values differ"):
        _writeback_values(config, root, mapping, labels, "cell_curator")
