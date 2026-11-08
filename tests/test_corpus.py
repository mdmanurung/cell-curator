from __future__ import annotations

import json
from pathlib import Path

import pytest

from cell_curator.contracts import ContractError, sha256
from cell_curator.corpus import (
    CorpusAssertion,
    CorpusLabel,
    CorpusMapping,
    CorpusSnapshot,
    CorpusSource,
    load_corpus,
    propose_corpus,
    review_corpus,
)
from cell_curator.ontology import CellOntology, OntologyTerm
from cell_curator.provenance import publish_json


def source(source_id: str = "study-a", **updates) -> CorpusSource:
    return CorpusSource.model_validate({
        "source_id": source_id, "name": "Study A", "url": "https://example.org/study-a",
        "version": "1", "sha256": "a" * 64, "retrieved_at": "2026-11-07",
        **updates,
    })


def label(label_id: str = "T", **updates) -> CorpusLabel:
    return CorpusLabel.model_validate({
        "label_id": label_id, "name": "T cell", "description": "Measured T cell identity.",
        "source_id": "study-a", "source_label": "T cells", "species": ["Homo sapiens"],
        "positive": ["CD3D", "TRAC"], "negative": ["MS4A1"],
        "status": "reviewed", "reviewer": "curator", "rationale": "Reviewed source criteria.",
        **updates,
    })


def snapshot(*labels: CorpusLabel, **updates) -> CorpusSnapshot:
    return CorpusSnapshot.model_validate({
        "corpus_id": "example", "revision": "v1", "sources": [source()],
        "labels": list(labels) or [label()], **updates,
    })


def save(path: Path, value: CorpusSnapshot) -> Path:
    publish_json(path, value.model_dump(mode="json"))
    return path


def test_snapshot_roundtrip_checksum_and_immutable_records(tmp_path: Path) -> None:
    base = snapshot()
    path = save(tmp_path / "v1.json", base)
    assert load_corpus(path, sha256(path)) == base
    with pytest.raises(ContractError, match="checksum mismatch"):
        load_corpus(path, "b" * 64)
    with pytest.raises(ValueError, match="frozen"):
        base.labels[0].name = "Changed"
    with pytest.raises(ContractError, match="does not exist"):
        load_corpus(tmp_path / "missing.json")
    path.write_text('{"surprise": true}')
    with pytest.raises(ContractError, match="invalid corpus"):
        load_corpus(path)


@pytest.mark.parametrize("changes, message", [
    ({"label_id": "  "}, "at least 1 character"),
    ({"name": " "}, "at least 1 character"),
    ({"reviewer": ""}, "reviewer and rationale"),
    ({"species": []}, "explicit species"),
    ({"positive": ["CD3D"], "negative": ["cd3d"]}, "contradictory"),
    ({"positive": [""]}, "blank strings"),
    ({"cl_id": "CL:bad"}, "pattern"),
    ({"ontology_relation": "exact"}, "requires cl_id"),
    ({"study_id": "study-local"}, "accepted parent"),
    ({"measurement": "protein_expression"}, "feature_space"),
    ({"extra_field": True}, "Extra inputs"),
])
def test_label_validation(changes: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        label(**changes)


def test_duplicate_unknown_and_cyclic_references() -> None:
    with pytest.raises(ValueError, match="duplicate corpus"):
        snapshot(label(), label())
    with pytest.raises(ValueError, match="unknown corpus source"):
        snapshot(label(source_id="missing"))
    with pytest.raises(ValueError, match="unknown parent"):
        snapshot(label(parent_id="missing", level="L2"))
    with pytest.raises(ValueError, match="cycle"):
        snapshot(label("A", parent_id="B"), label("B", parent_id="A"))
    with pytest.raises(ValueError, match="child level"):
        snapshot(label("A"), label("B", parent_id="A"))
    with pytest.raises(ValueError, match="parents must be cell identities"):
        snapshot(label("A", axis="cell_state"), label("B", parent_id="A", level="L2"))
    with pytest.raises(ValueError, match="unknown source"):
        snapshot(sources=[source(depends_on=["missing"])])


def test_context_scope_and_ancestor_review_are_required() -> None:
    child = label("CD8", name="CD8 T cell", parent_id="T", level="L2", tissues=["blood"])
    base = snapshot(label(), child)
    assert base.select() == []
    assert [item.label_id for item in base.select({"species": "human"})] == ["T"]
    assert len(base.select({"species": "NCBITaxon:9606", "tissue": "blood"})) == 2
    assert base.select({"species": "mouse", "tissue": "blood"}) == []
    imported_parent = label(status="imported", reviewer="", rationale="")
    assert snapshot(imported_parent, child).select({"species": "human", "tissue": "blood"}) == []
    assert len(base.select({"species": "human"}, callable_only=False)) == 1
    with pytest.raises(ContractError, match="unreviewed"):
        base.programs({"species": "human"}, callable_only=False)


def test_deterministic_export_preserves_parent_and_source_ids() -> None:
    child = label("CD8", name="CD8 T cell", parent_id="T", level="L2", positive=["CD8B"])
    a, b = snapshot(label(), child), snapshot(child, label())
    assert a.digest == b.digest
    assert a.programs({"species": "human"}) == b.programs({"species": "human"})
    programs = {item.label: item for item in a.programs({"species": "human"})}
    assert programs["CD8 T cell"].parent == "T cell"
    assert programs["T cell"].positive == ("CD3D", "TRAC")
    assert programs["T cell"].source == f"corpus:example:{a.digest}:T"


def test_exact_ontology_export_is_checked_and_broader_mapping_stays_rich() -> None:
    ontology = CellOntology({"CL:0000084": OntologyTerm("CL:0000084", "T cell")}, source="fixture")
    exact = snapshot(label(cl_id="CL:0000084", ontology_relation="exact"))
    assert exact.programs({"species": "human"}, ontology=ontology)[0].cl_id == "CL:0000084"
    assert exact.programs({"species": "human"})[0].cl_id is None
    wrong = snapshot(label(name="B cell", cl_id="CL:0000084", ontology_relation="exact"))
    with pytest.raises(ContractError, match="inconsistent"):
        wrong.programs({"species": "human"}, ontology=ontology)
    broader = snapshot(label(cl_id="CL:0000084", ontology_relation="broader"))
    assert broader.labels[0].cl_id == "CL:0000084"
    assert broader.programs({"species": "human"}, ontology=ontology)[0].cl_id is None


def test_ontology_parent_must_be_semantically_nested() -> None:
    parent = label(cl_id="CL:0000084", ontology_relation="exact")
    child = label("B", name="B cell", parent_id="T", level="L2", cl_id="CL:0000236", ontology_relation="exact")
    ontology = CellOntology({
        "CL:0000084": OntologyTerm("CL:0000084", "T cell"),
        "CL:0000236": OntologyTerm("CL:0000236", "B cell"),
    }, source="fixture")
    with pytest.raises(ContractError, match="contradicts its parent"):
        snapshot(parent, child).programs({"species": "human"}, ontology=ontology)


def test_exact_mapping_is_explicit_and_retains_all_provenance() -> None:
    other = label("other:T", source_id="study-b", positive=["CD3G"])
    mapping = CorpusMapping(source_label_id="other:T", target_label_id="T", relation="exact", rationale="Same criteria.", reviewer="expert")
    base = snapshot(label(), other, sources=[source(), source("study-b")], mappings=[mapping])
    programs = base.programs({"species": "human"})
    assert len(programs) == 1
    assert programs[0].positive == ("CD3D", "CD3G", "TRAC")
    assert len(base.labels) == 2
    assert len(base.sources) == 2
    with pytest.raises(ContractError, match="explicit exact mapping"):
        snapshot(label(), other, sources=base.sources).programs({"species": "human"})
    with pytest.raises(ValueError, match="reviewer"):
        snapshot(label(), other, sources=base.sources, mappings=[mapping.model_copy(update={"reviewer": ""})])


def test_conflicting_exact_mapping_cannot_overrule_negative_evidence() -> None:
    other = label("other:T", source_id="study-b", positive=["MS4A1"], negative=[])
    mapping = CorpusMapping(source_label_id="other:T", target_label_id="T", relation="exact", rationale="Test conflict.", reviewer="expert")
    base = snapshot(label(), other, sources=[source(), source("study-b")], mappings=[mapping])
    with pytest.raises(ContractError, match="conflicting reviewed"):
        base.programs({"species": "human"})


def test_study_extension_and_measurement_are_explicit() -> None:
    state = label("local:state", name="PD1+TIGIT+", parent_id="T", level="L4", axis="cell_state", study_id="s1", positive=["PD1", "TIGIT"], negative=[], measurement="protein_expression", feature_space="protein")
    base = snapshot(label(), state)
    assert len(base.select({"species": "human"})) == 1
    context = {"species": "human", "study_id": "s1"}
    programs = base.programs(context, measurements=["gene_expression", "protein_expression"])
    assert programs[0].axis == "cell_type"
    assert programs[1].axis == "cell_state"
    with pytest.raises(ContractError, match="measurement is unavailable"):
        base.programs(context, measurements=["gene_expression"])
    isoform = label(measurement="transcript_isoform", feature_space="transcript")
    with pytest.raises(ContractError, match="no executable marker backend"):
        snapshot(isoform).programs({"species": "human"})


def test_observed_marker_and_threshold_cannot_masquerade_as_known_marker() -> None:
    with pytest.raises(ValueError, match="run_id"):
        CorpusAssertion(assertion_id="a", source_id="study-a", statement="Observed.", kind="observed_marker")
    with pytest.raises(ValueError, match="units"):
        CorpusAssertion(assertion_id="a", source_id="study-a", statement="Threshold.", kind="threshold", threshold=0.5)
    assertion = CorpusAssertion(assertion_id="a", source_id="study-a", statement="Observed.", kind="observed_marker", run_id="run1", features=["CD3D"])
    assert snapshot(label(assertions=[assertion])).labels[0].assertions[0].kind == "observed_marker"


def test_proposal_never_changes_base_and_requires_explicit_review(tmp_path: Path) -> None:
    base = save(tmp_path / "v1.json", snapshot())
    original = base.read_bytes()
    proposal_path = tmp_path / "proposal.json"
    proposal = propose_corpus(base, proposal_path, labels=[label("CD8", name="CD8 T cell", parent_id="T", level="L2")], proposer="analyst", reason="Independent source supplied.")
    assert proposal.labels[0].status == "proposed"
    assert base.read_bytes() == original
    assert len(load_corpus(base).programs({"species": "human"})) == 1
    output = tmp_path / "v2.json"
    receipt = review_corpus(base, proposal_path, output, reviewer="human", rationale="Reviewed source and measurements.", decision="accept", revision="v2")
    revised = load_corpus(output)
    assert len(revised.programs({"species": "human"})) == 2
    assert revised.parent_revision == "v1"
    assert revised.review and revised.review.reviewer == "human"
    assert base.read_bytes() == original
    assert receipt["output_sha256"] == sha256(output)
    assert output.with_suffix(".json.receipt.json").is_file()
    assert review_corpus(base, proposal_path, output, reviewer="human", rationale="Reviewed source and measurements.", decision="accept", revision="v2") == receipt
    with pytest.raises(ContractError, match="immutable"):
        review_corpus(base, proposal_path, output, reviewer="someone-else", rationale="Other review.", decision="accept", revision="v2")


@pytest.mark.parametrize("decision, count, status", [("reject", 1, "reviewed"), ("dispute", 2, "disputed")])
def test_rejected_and_disputed_proposals_do_not_enter_callable_vocabulary(tmp_path: Path, decision: str, count: int, status: str) -> None:
    base = save(tmp_path / "v1.json", snapshot())
    proposal_path = tmp_path / "proposal.json"
    propose_corpus(base, proposal_path, labels=[label("B", name="B cell")], proposer="analyst", reason="Test review.")
    output = tmp_path / "v2.json"
    review_corpus(base, proposal_path, output, reviewer="human", rationale="Evidence does not resolve it.", decision=decision, revision="v2")
    revised = load_corpus(output)
    assert len(revised.labels) == count
    assert len(revised.programs({"species": "human"})) == 1
    if decision == "dispute":
        assert next(item for item in revised.labels if item.label_id == "B").status == status


def test_stale_proposal_and_input_overwrite_fail_closed(tmp_path: Path) -> None:
    base = save(tmp_path / "v1.json", snapshot())
    proposal = tmp_path / "proposal.json"
    propose_corpus(base, proposal, labels=[label(description="Changed definition.")], proposer="analyst", reason="Correction.")
    with pytest.raises(ContractError, match="cannot overwrite"):
        review_corpus(base, proposal, base, reviewer="human", rationale="Review.", decision="accept", revision="v2")
    with pytest.raises(ContractError, match="new nonblank revision"):
        review_corpus(base, proposal, tmp_path / "new.json", reviewer="human", rationale="Review.", decision="accept", revision="v1")
    data = json.loads(base.read_text())
    data["labels"][0]["description"] = "Changed after proposal."
    base.write_text(json.dumps(data))
    with pytest.raises(ContractError, match="different or changed base"):
        review_corpus(base, proposal, tmp_path / "new.json", reviewer="human", rationale="Review.", decision="accept", revision="v2")


def test_alias_and_scope_conflicts_cannot_be_silently_normalized() -> None:
    with pytest.raises(ContractError, match="ambiguous corpus alias"):
        snapshot(label(), label("B", name="B cell", synonyms=["T cell"])).programs({"species": "human"})
    mapping = CorpusMapping(source_label_id="mouse:T", target_label_id="T", relation="exact", rationale="Test incompatible scope.", reviewer="expert")
    with pytest.raises(ValueError, match="incompatible biological scope"):
        snapshot(label(), label("mouse:T", species=["Mus musculus"]), mappings=[mapping])


def test_proposal_cannot_overwrite_base_or_repeat_changes(tmp_path: Path) -> None:
    base = save(tmp_path / "v1.json", snapshot())
    with pytest.raises(ContractError, match="cannot overwrite"):
        propose_corpus(base, base, labels=[label()], reason="Test.", proposer="analyst")
    with pytest.raises(ValueError, match="duplicate changes"):
        propose_corpus(base, tmp_path / "p.json", labels=[label(), label()], reason="Test.", proposer="analyst")


def test_reviewed_identity_levels_cannot_drift_from_exported_parent_depth() -> None:
    with pytest.raises(ValueError, match="hierarchy depth"):
        snapshot(label(level="L2"))
    with pytest.raises(ValueError, match="hierarchy depth"):
        snapshot(label(), label("child", name="Child identity", level="L3", parent_id="T"))
    state = label("state", name="State", axis="cell_state", level="L4", parent_id="T")
    assert len(snapshot(label(), state).labels) == 2
    imported = label("source-level", level="L8", status="imported", reviewer="", rationale="")
    assert snapshot(imported).labels[0].level == "L8"
