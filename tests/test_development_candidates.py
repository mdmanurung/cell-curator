"""annotation.development_candidates: unreviewed records run on development splits, never written back."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml
from test_annotation_pipeline import annotation_fixture
from test_corpus_isolation import SPLIT
from test_corpus_writeback import _completed

from cell_curator.annotation import annotation_contract
from cell_curator.config import load_config
from cell_curator.contracts import ContractError
from cell_curator.corpus import CorpusLabel, CorpusSnapshot, CorpusSource
from cell_curator.writeback import validate_writeback

STATUSES = ("reviewed", "imported", "proposed", "disputed", "superseded")


def _snapshot() -> CorpusSnapshot:
    labels = []
    for status in STATUSES:
        extra = {"reviewer": "r", "rationale": "x"} if status == "reviewed" else {}
        if status == "superseded":
            extra = {"superseded_by": "reviewed"}
        labels.append(CorpusLabel(
            label_id=status, name=status.title(), description="Synthetic.", source_id="s",
            source_label="s", species=["human"], tissues=["blood"], positive=["G"],
            status=status, **extra,
        ))
    return CorpusSnapshot(
        corpus_id="c", revision="r",
        sources=(CorpusSource(source_id="s", name="s", url="https://example.org", version="1",
                              sha256="a" * 64, retrieved_at="2026-11-10"),),
        labels=tuple(labels),
    )


def test_flag_off_selects_only_reviewed() -> None:
    ids = {label.label_id for label in _snapshot().select({"species": "human", "tissue": "blood"})}
    assert ids == {"reviewed"}


def test_flag_on_admits_imported_and_proposed_not_disputed_or_superseded() -> None:
    chosen = _snapshot().select({"species": "human", "tissue": "blood"}, development_candidates=True)
    assert {label.label_id for label in chosen} == {"reviewed", "imported", "proposed"}


def _flagged_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, isolation: bool) -> Path:
    config_path, _, _ = annotation_fixture(tmp_path, monkeypatch)
    value = yaml.safe_load(config_path.read_text())
    value["annotation"]["development_candidates"] = True
    if isolation:
        value["annotation"]["isolation"] = {"donor_split": "donor_split.tsv", "test_split": "B"}
        (config_path.parent / "donor_split.tsv").write_text(SPLIT)
    config_path.write_text(yaml.safe_dump(value))
    return config_path


def test_flag_without_isolation_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config_path = _flagged_config(tmp_path, monkeypatch, isolation=False)
    with pytest.raises(ContractError, match="requires annotation.isolation"):
        annotation_contract(load_config(config_path))


def test_flag_with_isolation_is_recorded_in_the_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config_path = _flagged_config(tmp_path, monkeypatch, isolation=True)
    contract = annotation_contract(load_config(config_path))
    assert contract["development_candidates"] is True
    assert contract["rules"]["development_candidates"] is True
    assert contract["corpora"][0]["isolation"]["status"] == "isolated"


def test_writeback_refuses_a_flagged_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config, root, _ = _completed(tmp_path, monkeypatch)
    target = tmp_path / "annotated.h5ad"
    flagged = config.model_copy(update={
        "annotation": config.annotation.model_copy(update={"development_candidates": True})})
    with pytest.raises(ContractError, match="can never be written back"):
        validate_writeback(flagged, root, output_path=target)
    # A config edited after the run still meets the flag recorded in the plan.
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    plan_path.write_text(json.dumps({**plan, "development_candidates": True}))
    with pytest.raises(ContractError, match="can never be written back"):
        validate_writeback(config, root, output_path=target)
