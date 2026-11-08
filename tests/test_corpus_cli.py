from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from cell_curator.cli import main
from cell_curator.corpus import CorpusLabel, CorpusSnapshot, CorpusSource, load_corpus
from cell_curator.provenance import publish_json


@pytest.fixture
def corpus_path(tmp_path: Path) -> Path:
    source = CorpusSource(source_id="fixture", name="Synthetic test source", url="https://example.org", version="1", sha256="a" * 64, retrieved_at="2026-11-07")
    label = CorpusLabel(label_id="a", name="Type A", description="Synthetic identity.", source_id="fixture", source_label="source A", species=["human"], positive=["A1"], negative=["B1"], status="reviewed", reviewer="fixture-human", rationale="Synthetic test fixture.")
    snapshot = CorpusSnapshot(corpus_id="test", revision="v1", sources=[source], labels=[label])
    path = tmp_path / "corpus.json"
    publish_json(path, snapshot.model_dump(mode="json"))
    return path


def run_cli(monkeypatch, capsys, *args: object) -> dict:
    monkeypatch.setattr(sys, "argv", ["cell-curator", *map(str, args)])
    main()
    output = capsys.readouterr()
    assert output.err == ""
    return json.loads(output.out)


def test_corpus_cli_validate_and_export(monkeypatch, capsys, corpus_path: Path, tmp_path: Path) -> None:
    result = run_cli(monkeypatch, capsys, "corpus", "validate", "--input", corpus_path)
    assert result["status"] == "complete"
    assert result["reviewed"] == result["labels"] == 1
    assert result["content_digest"] == load_corpus(corpus_path).digest
    context = tmp_path / "context.json"
    context.write_text('{"species":"Homo sapiens"}')
    output = tmp_path / "markers.json"
    exported = run_cli(monkeypatch, capsys, "corpus", "export", "--input", corpus_path, "--context", context, "--measurement", "gene_expression", "--output", output)
    assert exported["programs"] == 1
    programs = json.loads(output.read_text())
    assert programs[0]["label"] == "Type A"
    assert programs[0]["positive"] == ["A1"]
    assert programs[0]["source_version"] == "v1"


def test_corpus_cli_propose_review_diff_preserve_base(monkeypatch, capsys, corpus_path: Path, tmp_path: Path) -> None:
    original = corpus_path.read_bytes()
    changes = tmp_path / "changes.json"
    label = load_corpus(corpus_path).labels[0].model_dump(mode="json")
    label["description"] = "Corrected synthetic definition."
    changes.write_text(json.dumps({"labels": [label]}))
    proposal = tmp_path / "proposal.json"
    proposed = run_cli(monkeypatch, capsys, "corpus", "propose", "--input", corpus_path, "--changes", changes, "--output", proposal, "--reason", "Definition correction.", "--proposer", "analyst")
    assert proposed["labels"][0]["status"] == "proposed"
    revised = tmp_path / "v2.json"
    review = run_cli(monkeypatch, capsys, "corpus", "review", "--input", corpus_path, "--proposal", proposal, "--output", revised, "--decision", "accept", "--revision", "v2", "--reviewer", "fixture-human", "--rationale", "Checked correction.")
    assert review["revision"] == "v2"
    assert review["label_changes"][0]["label_id"] == "a"
    diff = run_cli(monkeypatch, capsys, "corpus", "diff", "--input", corpus_path, "--other", revised)
    assert diff["changed"] == ["a"]
    assert diff["added"] == diff["removed"] == []
    assert corpus_path.read_bytes() == original
    assert load_corpus(revised).labels[0].description == "Corrected synthetic definition."


def test_corpus_cli_import_is_local_and_unreviewed(monkeypatch, capsys, corpus_path: Path, tmp_path: Path) -> None:
    records = tmp_path / "saved-cap.json"
    records.write_text(json.dumps({"context": {"species": "human"}, "labels": [{"id": "b", "cell_label": "Type B", "rationale": "Synthetic source definition.", "canonical_marker_genes": ["B1"]}]}))
    output = tmp_path / "imported.json"
    result = run_cli(monkeypatch, capsys, "corpus", "import", "--input", corpus_path, "--records", records, "--output", output, "--source-kind", "cap", "--source-id", "second", "--name", "Second synthetic source", "--url", "https://example.org/second", "--source-version", "1", "--retrieved-at", "2026-11-07")
    assert result["imported_count"] == 1
    assert result["review_required"] is True
    assert result["network_attempted"] is False
    snapshot = load_corpus(output)
    assert len(snapshot.labels) == 2
    assert next(item for item in snapshot.labels if item.source_id == "second").status == "imported"
    assert len(snapshot.programs({"species": "human"})) == 1


@pytest.mark.parametrize("kind", ["bad-changes", "bad-context"])
def test_corpus_cli_refuses_invalid_trust_boundary_inputs(monkeypatch, capsys, corpus_path: Path, tmp_path: Path, kind: str) -> None:
    invalid = tmp_path / "invalid.json"
    output = tmp_path / "out.json"
    if kind == "bad-changes":
        invalid.write_text('{"reviewer":"agent","labels":[]}')
        args = ["propose", "--changes", invalid, "--proposer", "a", "--reason", "test"]
    else:
        invalid.write_text('["human"]')
        args = ["export", "--context", invalid]
    monkeypatch.setattr(sys, "argv", ["cell-curator", "corpus", *map(str, args), "--input", str(corpus_path), "--output", str(output)])
    with pytest.raises(SystemExit) as stopped:
        main()
    assert stopped.value.code == 2
    response = capsys.readouterr()
    assert json.loads(response.err)["error"] == "ContractError"
    assert not output.exists()
