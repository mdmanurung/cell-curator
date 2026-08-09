from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from cell_curator import discover_candidates, validate_candidate_registry
from cell_curator.contracts import ContractError


def _write(frame: pd.DataFrame, path: Path) -> Path:
    frame.to_csv(path, sep="\t", index=False)
    return path


def _parents() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": [f"c{i}" for i in range(16)],
            "scope": ["immune"] * 16,
            "parent_id": ["0"] * 6 + ["1"] * 6 + ["2"] * 4,
            "donor_id": ["d1", "d2"] * 8,
            "capture_id": ["x1", "x2"] * 8,
            "closest_sibling": ["1"] * 6 + ["0"] * 6 + ["1"] * 4,
        }
    )


def test_persistent_stored_child_and_local_child_share_registry(
    tmp_path: Path, config: dict
) -> None:
    parents = _parents()
    stored_rows = []
    for resolution, child in ((1.0, "a"), (1.2, "b")):
        for barcode in ["c0", "c1"]:
            stored_rows.append(
                {
                    "barcode": barcode,
                    "scope": "immune",
                    "resolution": resolution,
                    "child_id": child,
                }
            )
        for barcode in ["c2", "c3", "c4", "c5"]:
            stored_rows.append(
                {
                    "barcode": barcode,
                    "scope": "immune",
                    "resolution": resolution,
                    "child_id": "rest0",
                }
            )
        for barcode in [f"c{i}" for i in range(6, 16)]:
            stored_rows.append(
                {
                    "barcode": barcode,
                    "scope": "immune",
                    "resolution": resolution,
                    "child_id": f"p{parents.set_index('barcode').loc[barcode, 'parent_id']}",
                }
            )
    local_rows = []
    # Resolution 0.1 is stable but biologically uninformative (one cluster).
    for barcode in [f"c{i}" for i in range(6, 12)]:
        local_rows.append(
            {
                "barcode": barcode,
                "scope": "immune",
                "parent_id": "1",
                "resolution": 0.1,
                "local_child_id": "0",
                "stable": True,
                "biologically_resolving": False,
                "technical_dominated": False,
                "eligible": True,
            }
        )
    for index, barcode in enumerate([f"c{i}" for i in range(6, 12)]):
        local_rows.append(
            {
                "barcode": barcode,
                "scope": "immune",
                "parent_id": "1",
                "resolution": 0.2,
                "local_child_id": "0" if index < 4 else "1",
                "stable": True,
                "biologically_resolving": True,
                "technical_dominated": False,
                "eligible": True,
            }
        )
    audit = pd.DataFrame(
        {
            "scope": ["immune"] * 3,
            "parent_id": ["0", "1", "2"],
            "audit_flagged": [True, True, True],
        }
    )
    paths = {
        "parents": _write(parents, tmp_path / "parents.tsv"),
        "stored": _write(pd.DataFrame(stored_rows), tmp_path / "stored.tsv"),
        "local": _write(pd.DataFrame(local_rows), tmp_path / "local.tsv"),
        "audit": _write(audit, tmp_path / "audit.tsv"),
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    manifest = discover_candidates.run_discovery(
        config_path=config_path,
        parents_path=paths["parents"],
        stored_path=paths["stored"],
        local_path=paths["local"],
        audit_path=paths["audit"],
        output_dir=tmp_path / "discovery",
    )
    registry = pd.read_csv(tmp_path / "discovery/candidate_registry.tsv", sep="\t", dtype=str)
    membership = pd.read_csv(
        tmp_path / "discovery/candidate_membership.tsv", sep="\t", dtype=str
    )
    triage = pd.read_csv(tmp_path / "discovery/parent_triage.tsv", sep="\t", dtype=str)

    assert set(registry["source_kind"]) == {"stored_resolution", "parent_local"}
    local = registry[registry["source_kind"].eq("parent_local")]
    assert set(local["source_resolution"]) == {"0.2"}
    assert manifest["invariants"]["stable_one_cluster_rejected"] is True
    assert (
        triage.set_index("parent_id").loc["2", "triage_action"]
        == "REJECT_NO_DISCRETE_CANDIDATE_EVIDENCE"
    )
    validated = validate_candidate_registry.validate_registry(
        config, audit.astype(str), registry, membership
    )
    assert validated["invariants"]["biological_labels_blank"] is True


def test_registry_rejects_prefilled_biological_label(config: dict) -> None:
    audit = pd.DataFrame({"scope": ["s"], "parent_id": ["0"], "audit_flagged": [True]})
    registry = pd.DataFrame(
        {
            "scope": ["s"],
            "candidate_key": ["x"],
            "parent_id": ["0"],
            "candidate_child_uid": ["s::0::candidate_1"],
            "source_kind": ["stored_resolution"],
            "n_cells_in_parent": [2],
            "n_parent_remainder": [2],
            "n_donors": [2],
            "n_captures": [2],
            "stability_pass": [True],
            "eligibility_counts_pass": [True],
            "closest_sibling": ["1"],
            "membership_sha256": ["wrong"],
            "triage_action": ["FULL_TESTING"],
            "biological_label_L1": ["forbidden"],
        }
    )
    membership = pd.DataFrame(
        {
            "scope": ["s"] * 4,
            "candidate_key": ["x"] * 4,
            "parent_id": ["0"] * 4,
            "barcode": ["a", "b", "c", "d"],
            "candidate_membership": [True, True, False, False],
        }
    )
    with pytest.raises(ContractError, match="must remain blank"):
        validate_candidate_registry.validate_registry(config, audit, registry, membership)
