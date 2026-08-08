from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

import build_review_packet
import validate_final_mapping
from contracts import ContractError, atomic_write_text, sha256


def _mapping_inputs():
    canonical = pd.DataFrame(
        {"barcode": list("abcdefgh"), "scope": ["s"] * 8, "parent_id": ["0"] * 6 + ["1"] * 2}
    )
    registry = pd.DataFrame(
        {
            "scope": ["s", "s", "s"],
            "candidate_key": ["pass1", "fail", "pass2"],
            "parent_id": ["0", "0", "0"],
            "candidate_child_uid": ["s::0::candidate_1", "s::0::candidate_2", "s::0::candidate_3"],
        }
    )
    memberships = []
    selected = {"pass1": {"a", "b"}, "fail": {"c"}, "pass2": {"d", "e"}}
    for candidate, members in selected.items():
        for barcode in list("abcdef"):
            memberships.append(
                {"scope": "s", "candidate_key": candidate, "parent_id": "0", "barcode": barcode,
                 "candidate_membership": barcode in members}
            )
    membership = pd.DataFrame(memberships)
    decisions = pd.DataFrame(
        {
            "scope": ["s", "s"],
            "parent_id": ["0", "1"],
            "decision_outcome": ["ACCEPT SPLIT", "RETAIN PARENT"],
            "accepted_child_ids": ["s::0::candidate_1;s::0::candidate_3", ""],
        }
    )
    return canonical, decisions, registry, membership


def test_mapping_uses_each_accepted_child_and_folds_failed_candidate_to_remainder() -> None:
    canonical, decisions, registry, membership = _mapping_inputs()
    mapping = validate_final_mapping.construct_mapping(
        canonical=canonical, decisions=decisions, registry=registry, membership=membership
    )
    indexed = mapping.set_index("barcode")
    assert indexed.loc["a", "leaf_id"] == "s::0::candidate_1"
    assert indexed.loc["d", "leaf_id"] == "s::0::candidate_3"
    assert indexed.loc["c", "leaf_id"] == "s::0::remainder"
    assert indexed.loc["f", "leaf_id"] == "s::0::remainder"
    assert indexed.loc["g", "leaf_id"] == "s::1"
    result = validate_final_mapping.validate_mapping(canonical, mapping)
    assert result["n_cells"] == 8
    assert result["n_global_leaf_identities"] == 4
    assert result["invariants"]["parent_counts_conserved"] is True


def test_mapping_rejects_overlapping_accepted_children() -> None:
    canonical, decisions, registry, membership = _mapping_inputs()
    membership.loc[
        membership["candidate_key"].eq("pass2") & membership["barcode"].eq("a"),
        "candidate_membership",
    ] = True
    with pytest.raises(ContractError, match="overlap"):
        validate_final_mapping.construct_mapping(
            canonical=canonical, decisions=decisions, registry=registry, membership=membership
        )


def _minimal_packet_run(root: Path, config: dict, *, label: str = "") -> None:
    text_files = {
        "context.md": "context\n",
        "assumptions.md": "approved assumptions\n",
        "plan.md": "plan\n",
    }
    for relative, text in text_files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    state = {
        "mode": config["run"]["mode"],
        "gate_a_approved": False,
        "knowledgebase_available": False,
        "critic_reconciled": True,
        "human_approval": False,
    }
    (root / "run_state.json").write_text(json.dumps(state))
    json_files = (
        "input_contract.manifest.json",
        "discovery/discovery.manifest.json",
        "gpu_markers/gpu_receipts.validated.json",
        "final_provenance.manifest.json",
    )
    for relative in json_files:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"status": "complete"}))
    generic = pd.DataFrame({"scope": ["s"], "parent_id": ["0"]})
    generic.to_csv(root / "adaptive_parent_audit.tsv", sep="\t", index=False)
    registry = pd.DataFrame(
        {"scope": ["s"], "candidate_key": ["c"], "parent_id": ["0"],
         "biological_label_L1": [label], "cl_id": [""]}
    )
    registry_path = root / "discovery/candidate_registry.tsv"
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    registry.to_csv(registry_path, sep="\t", index=False)
    pd.DataFrame({"scope": ["s"], "parent_id": ["0"], "triage_action": ["FULL_TESTING"]}).to_csv(
        root / "discovery/parent_triage.tsv", sep="\t", index=False
    )
    (root / "contrasts").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"lane_id": ["x"]}).to_csv(
        root / "contrasts/comparator_manifest.tsv", sep="\t", index=False
    )
    for relative in (
        "parent_decision_summary.tsv",
        "child_evidence_summary.tsv",
        "annotation_review_template.tsv",
    ):
        path = root / relative
        pd.DataFrame({"scope": ["s"], "parent_id": ["0"], "biological_label_L1": [label], "cl_id": [""]}).to_csv(
            path, sep="\t", index=False
        )
    pd.DataFrame({"barcode": ["a"], "leaf_id": ["s::0"]}).to_csv(
        root / "parent_child_mapping.tsv", sep="\t", index=False
    )
    pd.DataFrame({"scope": ["s"], "leaf_total": [1]}).to_csv(
        root / "leaf_identity_totals.tsv", sep="\t", index=False
    )
    pd.DataFrame({"finding_id": ["f1"], "severity": ["BLOCKER"]}).to_csv(
        root / "critic_findings.tsv", sep="\t", index=False
    )
    pd.DataFrame({"finding_id": ["f1"], "status": ["RESOLVED"]}).to_csv(
        root / "critic_reconciliation.tsv", sep="\t", index=False
    )


def test_evidence_only_packet_prohibits_labels_and_preserves_sources(
    tmp_path: Path, config: dict
) -> None:
    run = tmp_path / "run"
    run.mkdir()
    _minimal_packet_run(run, config)
    input_hash = sha256(run / "input_contract.manifest.json")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    packet = build_review_packet.build_packet(config_path, run)

    assert (packet / "review_packet.manifest.json").is_file()
    assert sha256(run / "input_contract.manifest.json") == input_hash
    assert json.loads((packet / "review_packet.manifest.json").read_text())["writeback_permitted"] is False
    packet_hash = sha256(packet / "review_packet.manifest.json")
    with pytest.raises(ContractError, match="overwrite"):
        build_review_packet.build_packet(config_path, run)
    assert sha256(packet / "review_packet.manifest.json") == packet_hash


def test_evidence_only_packet_rejects_biological_labels(tmp_path: Path, config: dict) -> None:
    run = tmp_path / "run"
    run.mkdir()
    _minimal_packet_run(run, config, label="forbidden")
    with pytest.raises(ContractError, match="must remain blank"):
        build_review_packet.validate_review_gate(config, run)


def test_atomic_writer_refuses_resume_overwrite(tmp_path: Path) -> None:
    path = tmp_path / "immutable.txt"
    atomic_write_text(path, "first\n")
    with pytest.raises(ContractError, match="overwrite"):
        atomic_write_text(path, "second\n")
    assert path.read_text() == "first\n"


def test_atomic_writer_can_publish_after_unrelated_stale_partial(tmp_path: Path) -> None:
    path = tmp_path / "artifact.txt"
    stale = tmp_path / ".artifact.txt.partial.previous-job"
    stale.write_text("incomplete\n")

    atomic_write_text(path, "complete\n")

    assert path.read_text() == "complete\n"
    assert stale.read_text() == "incomplete\n"
