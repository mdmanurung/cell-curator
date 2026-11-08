"""scGate_models import: decision rows 23 to 25, 28 and 29."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from cell_curator.contracts import ContractError
from cell_curator.corpus import CorpusSnapshot, CorpusSource, load_corpus
from cell_curator.corpus_import import import_corpus_records
from cell_curator.corpus_scgate import scgate_payload
from cell_curator.provenance import publish_json


def _git(folder: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-C", str(folder), "-c", "user.name=t", "-c", "user.email=t@t", *args],
        check=True, capture_output=True,
    )


@pytest.fixture
def clone(tmp_path: Path) -> Path:
    root = tmp_path / "scGate_models"
    folder = root / "human" / "PBMC"
    folder.mkdir(parents=True)
    (root / "CellOntology_dictionary.tsv").write_text(
        "scGate_multi\tCellOntology_name\tCellOntology_ID\nNK\tnatural killer cell\tCL:0000623\n"
    )
    (folder / "master_table.tsv").write_text(
        "name\tsignature\nImmune\tPTPRC\nNK\tKLRD1;NKG7;CD3D-\nEpithelial\tKRT8\n"
    )
    (folder / "NK_scGate_Model.tsv").write_text(
        "levels\tuse_as\tname\tsignature\n"
        "level1\tpositive\tImmune\t\nlevel1\tnegative\tEpithelial\t\n"
        "level2\tpositive\tNK\t\n"
    )
    _git(root, "init", "-q")
    _git(root, "add", ".")
    _git(root, "commit", "-q", "-m", "models")
    _git(root, "tag", "v0.16")
    return root


def _base(tmp_path: Path) -> Path:
    base = CorpusSnapshot(
        corpus_id="c", revision="r0",
        sources=(CorpusSource(source_id="seed", name="seed", url="https://example.org", version="1",
                              sha256="a" * 64, retrieved_at="2026-11-08"),),
        labels=(),
    )
    path = tmp_path / "base.json"
    publish_json(path, base.model_dump(mode="json"))
    return path


def test_model_becomes_one_label_with_signed_genes_and_its_exclusions(clone: Path) -> None:
    payload = json.loads(scgate_payload(clone, ["PBMC"], species="human", tissue="blood", tag="v0.16"))
    [label] = payload["labels"]
    assert (label["name"], label["positive"], label["negative"]) == ("NK", ["KLRD1", "NKG7"], ["CD3D"])
    assert label["cl_id"] == "CL:0000623"  # decision row 29: a first guess for the ontology check
    assert label["scgate_excludes"] == ["Epithelial"]
    # Decision row 30: the excluded population travels with its own signature genes.
    assert label["exclusions"] == [{"name": "Epithelial", "positive": ["KRT8"], "negative": []}]
    assert label["scgate_gate_path"]["level1"] == {"positive": ["Immune"], "negative": ["Epithelial"]}
    assert payload["context"] == {"species": "human", "tissue": "blood"}


def test_clone_must_sit_at_the_named_tag_and_a_tissue_is_required(clone: Path) -> None:
    _git(clone, "commit", "-q", "--allow-empty", "-m", "newer")
    with pytest.raises(ContractError, match="not v0.16"):
        scgate_payload(clone, ["PBMC"], species="human", tissue="blood", tag="v0.16")
    with pytest.raises(ContractError, match="tissue scope"):
        scgate_payload(clone, ["PBMC"], species="human", tissue="", tag="v0.16")


def test_import_creates_unreviewed_labels_and_keeps_provenance(clone: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.json"
    receipt = import_corpus_records(
        _base(tmp_path), clone, out, source_kind="scgate", source_id="scgate-pbmc",
        name="scGate_models PBMC", url="https://github.com/carmonalab/scGate_models",
        version="v0.16", retrieved_at="2026-11-08",
        scgate={"collections": ["PBMC"], "tag": "v0.16", "species": "human", "tissue": "blood"},
    )
    assert receipt["review_required"] is True and receipt["imported_count"] == 1
    snapshot = load_corpus(out)
    [label] = [item for item in snapshot.labels if item.source_id == "scgate-pbmc"]
    assert label.status == "imported"
    assert snapshot.programs({"species": "human", "tissue": "blood"}) == []  # not callable until reviewed
    kept = json.loads(Path(str(out) + ".source.json").read_text())["labels"][0]
    assert len(kept["scgate_file_sha256"]) == 64
