from __future__ import annotations

from pathlib import Path

import pytest

from cell_curator.annotation import annotation_contract
from cell_curator.config import CellCuratorConfig, ProviderConfig
from cell_curator.contracts import ContractError, sha256
from cell_curator.corpus import CorpusLabel, CorpusSnapshot, CorpusSource
from cell_curator.knowledge import (
    CorpusProvider,
    collect_programs,
    corpus_context,
    merge_corpus_context,
)
from cell_curator.provenance import publish_json


def provider(tmp_path: Path, identifier: str, *, name: str = "Type A", positive: str = "A1", context: dict | None = None) -> CorpusProvider:
    snapshot = CorpusSnapshot(
        corpus_id=identifier, revision="v1",
        sources=[CorpusSource(source_id=identifier, name="Synthetic source", url="https://example.org", version="v1", sha256="a" * 64, retrieved_at="2026-11-07")],
        labels=[CorpusLabel(label_id=f"{identifier}:a", name=name, description="Synthetic definition.", source_id=identifier, source_label=name, species=["Homo sapiens"], tissues=["blood"], positive=[positive], status="reviewed", reviewer="fixture-curator", rationale="Synthetic fixture.")],
    )
    path = tmp_path / (identifier + ".json")
    publish_json(path, snapshot.model_dump(mode="json"))
    config = ProviderConfig(kind="corpus", path=str(path), sha256=sha256(path), version="v1", context=context or {})
    return CorpusProvider(config, context={"species": "human", "tissue": "blood"})


def test_corpus_provider_cannot_override_dataset_biology(tmp_path: Path) -> None:
    with pytest.raises(ContractError, match="contradicts dataset species"):
        provider(tmp_path, "a", context={"species": "mouse", "tissue": "blood"})
    assert len(provider(tmp_path, "b", context={"species": "Homo sapiens", "tissue": "blood"}).programs()) == 1
    for key, before, after in (("tissue", "blood", "lung"), ("stage", "adult", "embryonic"), ("condition", "healthy", "tumor"), ("study_id", "study-a", "study-b")):
        with pytest.raises(ContractError, match="contradicts dataset"):
            merge_corpus_context({key: before}, {key: after})
    for unknown in ("", "unknown"):
        with pytest.raises(ContractError):
            merge_corpus_context({"tissue": unknown}, {"tissue": "blood"})
    assert merge_corpus_context({"assay": ""}, {"assay": "rna"})["assay"] == "rna"


@pytest.mark.parametrize("reverse", [False, True])
def test_cross_corpus_labels_cannot_be_merged_without_reviewed_equivalence(tmp_path: Path, reverse: bool) -> None:
    a = provider(tmp_path, "a", positive="A1")
    b = provider(tmp_path, "b", positive="OTHER_MARKER")
    with pytest.raises(ContractError, match="explicit exact mapping"):
        collect_programs([b, a] if reverse else [a, b])


def test_noncolliding_corpus_sources_remain_distinct(tmp_path: Path) -> None:
    a = provider(tmp_path, "a")
    b = provider(tmp_path, "b", name="Type B", positive="B1")
    programs = collect_programs([a, b])
    assert [item.label for item in programs] == ["Type A", "Type B"]
    assert programs[0].source != programs[1].source


def test_annotation_plan_uses_the_same_context_guard(config: dict, tmp_path: Path) -> None:
    corpus = provider(tmp_path, "a")
    path = tmp_path / "a.json"
    value = CellCuratorConfig.model_validate(config)
    value.guidance.context.organism = "mouse"
    value.knowledge.providers = [ProviderConfig(kind="corpus", path=str(path), sha256=sha256(path), version=corpus.snapshot.revision, context={"species": "human", "tissue": "blood"})]
    with pytest.raises(ContractError, match="contradicts dataset species"):
        annotation_contract(value)


def test_provider_rejects_changed_snapshot_or_revision(tmp_path: Path) -> None:
    provider(tmp_path, "a")
    path = tmp_path / "a.json"
    with pytest.raises(ContractError, match="checksum mismatch"):
        CorpusProvider(ProviderConfig(kind="corpus", path=str(path), sha256="b" * 64, version="v1"))
    with pytest.raises(ContractError, match="revision differs"):
        CorpusProvider(ProviderConfig(kind="corpus", path=str(path), sha256=sha256(path), version="wrong"))


def test_annotation_context_uses_a_stable_assay_scope(config: dict) -> None:
    value = CellCuratorConfig.model_validate(config)
    value.annotation.assay = "transcriptome"
    assert corpus_context(value)["assay"] == "rna"
    value.annotation.corpus_assay = "10x-3prime-v3"
    assert corpus_context(value)["assay"] == "10x-3prime-v3"
