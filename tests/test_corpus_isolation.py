"""Decision row 21: a corpus learned from test donors is refused."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from test_annotation_pipeline import annotation_fixture

from cell_curator.annotation import annotation_contract
from cell_curator.config import load_config
from cell_curator.contracts import ContractError
from cell_curator.corpus import CorpusLabel, CorpusSnapshot, CorpusSource, check_test_isolation

SPLIT = "donor_id\tsplit\tn_cells\nd1\tA\t5\nd2\tA\t5\nd3\tB\t5\n"


def _snapshot(*derived: list[dict]) -> CorpusSnapshot:
    labels = tuple(
        CorpusLabel(
            label_id=f"l{i}", name=f"L{i}", description="Synthetic.", source_id="s", source_label="s",
            species=["human"], tissues=["blood"], positive=["G"], derived_from=items,
        )
        for i, items in enumerate(derived)
    )
    return CorpusSnapshot(
        corpus_id="c", revision="r",
        sources=(CorpusSource(source_id="s", name="s", url="https://example.org", version="1",
                              sha256="a" * 64, retrieved_at="2026-11-08"),),
        labels=labels,
    )


@pytest.fixture
def split(tmp_path: Path) -> Path:
    path = tmp_path / "donor_split.tsv"
    path.write_text(SPLIT)
    return path


def test_records_learned_from_split_a_pass_and_pin_their_donors(split: Path) -> None:
    snapshot = _snapshot([{"dataset_id": "aida", "donors": ["d1", "d2"]}], [])
    result = check_test_isolation(snapshot, split, test_split="B", dataset_id="aida")
    assert result["status"] == "isolated"
    assert snapshot.learned_donors() == ("d1", "d2")


def test_a_test_donor_is_refused(split: Path) -> None:
    snapshot = _snapshot([{"dataset_id": "aida", "donors": ["d2", "d3"]}])
    with pytest.raises(ContractError, match="learned from test donors"):
        check_test_isolation(snapshot, split, test_split="B", dataset_id="aida")


def test_citing_the_test_dataset_without_donors_is_refused(split: Path) -> None:
    snapshot = _snapshot([{"dataset_id": "aida"}])
    with pytest.raises(ContractError, match="without naming donors"):
        check_test_isolation(snapshot, split, test_split="B", dataset_id="aida")
    # A record from another dataset needs no donors.
    check_test_isolation(_snapshot([{"dataset_id": "other"}]), split, test_split="B", dataset_id="aida")


def test_unknown_donor_of_the_test_dataset_is_refused(split: Path) -> None:
    snapshot = _snapshot([{"dataset_id": "aida", "donors": ["d9"]}])
    with pytest.raises(ContractError, match="absent from the split file"):
        check_test_isolation(snapshot, split, test_split="B", dataset_id="aida")


def test_config_enforces_isolation_before_gate_a(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, split: Path
) -> None:
    config_path, _, _ = annotation_fixture(tmp_path, monkeypatch)
    value = yaml.safe_load(config_path.read_text())
    value["annotation"]["isolation"] = {"donor_split": "donor_split.tsv", "test_split": "B"}
    (config_path.parent / "donor_split.tsv").write_text(SPLIT)
    config_path.write_text(yaml.safe_dump(value))
    contract = annotation_contract(load_config(config_path))
    assert contract["corpora"][0]["isolation"]["status"] == "isolated"
