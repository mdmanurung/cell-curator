from __future__ import annotations

import json
from pathlib import Path

import pytest

from cell_curator.contracts import ContractError, sha256
from cell_curator.corpus import CorpusSnapshot, load_corpus
from cell_curator.corpus_import import import_corpus_records
from cell_curator.provenance import publish_json


@pytest.fixture
def empty_corpus(tmp_path: Path) -> Path:
    path = tmp_path / "base.json"
    snapshot = CorpusSnapshot(corpus_id="tests", revision="v0", sources=(), labels=())
    publish_json(path, snapshot.model_dump(mode="json"))
    return path


def import_saved(base: Path, source: Path, out: Path, **kwargs) -> dict:
    return import_corpus_records(
        base, source, out, source_kind="cap", source_id="source-a", name="Source A",
        url="https://example.org/source-a", version="1", retrieved_at="2026-11-07", **kwargs,
    )


def test_cap_import_preserves_reported_evidence_and_never_approves(empty_corpus: Path, tmp_path: Path) -> None:
    source = tmp_path / "cap.json"
    source.write_text(json.dumps({
        "context": {"species": "Homo sapiens", "tissue": "blood"},
        "labels": [{"id": "nk", "cell_label": "Natural killer (NK) cells",
                    "cell_fullname": "Natural killer cell", "rationale": "Source definition.",
                    "canonical_marker_genes": "GNLY,NKG7", "marker_gene_evidence": "NKG7;PRF1",
                    "ontology_term_id": "CL:0000623", "synonyms": ["NK cell"]}],
    }))
    output = tmp_path / "v1.json"
    receipt = import_saved(empty_corpus, source, output)
    snapshot = load_corpus(output)
    record = snapshot.labels[0]
    assert record.label_id == "source-a:nk"
    assert record.source_label == "Natural killer (NK) cells"
    assert record.positive == ("GNLY", "NKG7")
    assert record.negative == ()
    assert record.assertions[0].kind == "reported_marker"
    assert record.assertions[0].features == ("NKG7", "PRF1")
    assert record.ontology_relation == "unresolved"
    assert record.status == "imported"
    assert snapshot.programs({"species": "human", "tissue": "blood"}) == []
    assert receipt["source_sha256"] == sha256(source)
    assert output.with_suffix(".json.source.json").read_bytes() == source.read_bytes()
    assert receipt["review_required"] is True
    assert receipt["network_attempted"] is False
    assert receipt["status"] == "partial"
    assert "negative_markers_unreported" in receipt["omissions"][0]["fields"]


def test_import_source_levels_do_not_imply_canonical_hierarchy(empty_corpus: Path, tmp_path: Path) -> None:
    source = tmp_path / "allen.json"
    source.write_text(json.dumps({"labels": [
        {"id": "t", "name": "T cell", "level": "1", "description": "Source lineage.", "markers": ["TRAC"]},
        {"id": "cd8", "name": "Naive CD8 T cell", "level": "2", "parent_id": "t", "description": "Source child.", "markers": ["CCR7"]},
    ], "context": {"species": "human", "tissue": "blood"}}))
    output = tmp_path / "v1.json"
    receipt = import_corpus_records(empty_corpus, source, output, source_kind="allen", source_id="allen-example", name="Example", url="https://example.org/allen", version="1", retrieved_at="2026-11-07")
    child = next(item for item in load_corpus(output).labels if item.name == "Naive CD8 T cell")
    assert child.source_level == "2"
    assert child.source_parent_id == "t"
    assert child.parent_id is None
    assert child.level == "L1"
    assert any("canonical_parent_unreviewed" in item["fields"] for item in receipt["omissions"])


def test_partial_import_keeps_invalid_and_unknown_fields_in_source_copy(empty_corpus: Path, tmp_path: Path) -> None:
    source = tmp_path / "partial.json"
    source.write_text(json.dumps([
        {"id": "one", "name": "Type A", "description": "Definition.", "positive": ["A"], "unknown_source_field": 4},
        {"id": "bad", "name": ""},
        "invalid row",
    ]))
    output = tmp_path / "v1.json"
    receipt = import_saved(empty_corpus, source, output)
    assert receipt["imported_count"] == 1
    assert receipt["failed_count"] == 2
    assert receipt["input_count"] == 3
    saved = json.loads(output.with_suffix(".json.source.json").read_text())
    assert saved[0]["unknown_source_field"] == 4
    assert len(saved) == 3


def test_sources_have_equal_status_and_same_names_are_not_merged(empty_corpus: Path, tmp_path: Path) -> None:
    source = tmp_path / "one.json"
    source.write_text(json.dumps([{"id": "t", "name": "T cell", "description": "Definition.", "positive": ["CD3D"]}]))
    first = tmp_path / "v1.json"
    import_saved(empty_corpus, source, first)
    second = tmp_path / "v2.json"
    import_corpus_records(first, source, second, source_kind="allen", source_id="source-b", name="Source B", url="https://example.org/b", version="1", retrieved_at="2026-11-07")
    snapshot = load_corpus(second)
    assert len(snapshot.sources) == len(snapshot.labels) == 2
    assert {item.source_id for item in snapshot.labels} == {"source-a", "source-b"}
    assert {item.status for item in snapshot.labels} == {"imported"}
    assert snapshot.mappings == ()


def test_import_limits_and_invalid_json_fail_before_output(empty_corpus: Path, tmp_path: Path) -> None:
    source = tmp_path / "source.json"
    source.write_text('[{"name":"A"},{"name":"B"}]')
    output = tmp_path / "out.json"
    for bounds in ({"max_records": 1}, {"max_bytes": 2}, {"max_records": 0}):
        with pytest.raises(ContractError, match="limit|exceeds"):
            import_saved(empty_corpus, source, output, **bounds)
    assert not output.exists()
    source.write_text('{"data":{"labels":[]}}')
    with pytest.raises(ContractError, match="label list"):
        import_saved(empty_corpus, source, output)
    source.write_bytes(b"\xff")
    with pytest.raises(ContractError, match="UTF-8 JSON"):
        import_saved(empty_corpus, source, output)


def test_import_refuses_input_overwrite_duplicate_ids_and_source_redefinition(empty_corpus: Path, tmp_path: Path) -> None:
    source = tmp_path / "source.json"
    source.write_text('[{"id":"a","name":"A"},{"id":"a","name":"B"}]')
    with pytest.raises(ContractError, match="IDs repeat"):
        import_saved(empty_corpus, source, tmp_path / "v1.json")
    source.write_text('[{"id":"a","name":"A"}]')
    with pytest.raises(ContractError, match="cannot overwrite"):
        import_saved(empty_corpus, source, source)
    first = tmp_path / "v1.json"
    import_saved(empty_corpus, source, first)
    source.write_text('[{"id":"b","name":"B"}]')
    with pytest.raises(ContractError, match="different content"):
        import_saved(first, source, tmp_path / "v2.json")


def test_bundled_corpus_is_unreviewed_multi_source_example() -> None:
    root = Path(__file__).resolve().parents[1]
    path = root / "skills/cell-curator/assets/corpus/example.json"
    corpus = load_corpus(path)
    assert len(corpus.sources) >= 2
    assert {item.status for item in corpus.labels} == {"imported"}
    assert corpus.programs({"species": "human", "tissue": "blood"}) == []
    for source in corpus.sources:
        saved = path.parent / (source.source_id + ".json")
        assert sha256(saved) == source.sha256
