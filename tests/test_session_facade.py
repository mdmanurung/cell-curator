"""The object-oriented entry point must configure real runs, not a toy config."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from cell_curator import CellCurator
from cell_curator.contracts import ContractError


def _adata(seed: int = 5, n_per_cluster: int = 40) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    total = n_per_cluster * 2
    counts = np.zeros((total, 6), dtype=np.float32)
    counts[:n_per_cluster, :2] = rng.poisson(9, size=(n_per_cluster, 2)) + 1
    counts[n_per_cluster:, 2:4] = rng.poisson(9, size=(n_per_cluster, 2)) + 1
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    value = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {
                "leiden": ["0"] * n_per_cluster + ["1"] * n_per_cluster,
                "lineage": ["immune"] * total,
            },
            index=[f"cell-{index:03d}" for index in range(total)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    value.layers["counts"] = counts
    value.layers["logcounts"] = np.log1p(counts)
    value.obsm["X_pca_harmony"] = np.vstack(
        [
            rng.normal(-3, 0.2, size=(n_per_cluster, 3)),
            rng.normal(3, 0.2, size=(n_per_cluster, 3)),
        ]
    ).astype(np.float32)
    return value


MARKERS = {
    "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
    "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
}


def _curator(tmp_path: Path, **kwargs: object) -> CellCurator:
    defaults: dict[str, object] = {
        "adata": _adata(),
        "organism": "human",
        "tissue": "synthetic tissue",
        "run_id": "facade-v1",
        "cluster_key": "leiden",
        "markers": MARKERS,
        "output_root": str(tmp_path / "results"),
        "layer": "logcounts",
        "counts_layer": "counts",
    }
    defaults.update(kwargs)
    return CellCurator(**defaults)  # type: ignore[arg-type]


def _unattended(tmp_path: Path, **kwargs: object) -> CellCurator:
    """A run with every gate declared up front, as an unattended run requires."""

    defaults: dict[str, object] = {
        "preauthorized": True,
        "reviewer": "facade-fixture-reviewer",
        "assumptions": ["The synthetic marker programs and context are authorized."],
        "developmental_stage": "adult",
        "condition": "healthy",
        "experimental_context": "synthetic facade fixture",
    }
    defaults.update(kwargs)
    return _curator(tmp_path, **defaults)


def test_constructor_writes_a_configuration_that_strictly_validates(tmp_path: Path) -> None:
    curator = _curator(tmp_path)
    config = curator.validate()

    assert config.run.run_id == "facade-v1"
    assert config.input.keys.canonical_parent == "leiden"
    assert config.guidance.context.organism == "human"
    assert config.guidance.context.tissue == "synthetic tissue"
    # The declared partition must match what is actually in obs, not a template guess.
    partition = next(iter(config.canonical_partitions.values()))
    assert partition.expected_parent_count == 2
    assert "example_scope" not in config.canonical_partitions


def test_in_memory_input_is_written_and_frozen(tmp_path: Path) -> None:
    curator = _curator(tmp_path)
    source = Path(curator.validate().input.path)
    assert source.is_file(), "an in-memory AnnData must be materialized before it can be hashed"

    inspected = curator.inspect()
    assert inspected

    frozen = curator.freeze()
    assert frozen["qc"]
    assert (curator.run_directory / "prepared" / "cells.tsv").is_file()


def test_markers_mapping_is_persisted_as_a_local_provider(tmp_path: Path) -> None:
    curator = _curator(tmp_path)
    providers = curator.validate().knowledge.providers
    assert [provider.kind for provider in providers] == ["local_markers"]
    written = json.loads(Path(providers[0].path).read_text())
    assert written == MARKERS


def test_evidence_runs_through_the_same_pipeline_as_the_cli(tmp_path: Path) -> None:
    curator = _unattended(tmp_path)
    curator.prepare()
    evidence, bottom_up = curator.evidence()
    assert not evidence.empty
    assert not bottom_up.empty
    assert set(evidence["cluster_id"].astype(str)) == {"0", "1"}


def test_adaptive_refinement_audits_parents_and_proposes_bounded_candidates(
    tmp_path: Path,
) -> None:
    curator = _unattended(tmp_path)
    curator.prepare()

    audited = curator.audit_parents()
    assert len(audited) == 2
    assert {"parent_id", "n_cells"}.issubset(audited.columns)

    registry, _, _ = curator.propose_subclusters()
    # Canonical parents stay frozen: refinement proposes candidates, never
    # rewrites the parent partition.
    assert curator.validate().input.keys.canonical_parent == "leiden"
    assert (curator.run_directory / "comparators" / "comparator_manifest.tsv").is_file()
    assert registry is not None


def test_exactly_one_of_adata_or_path_is_required(tmp_path: Path) -> None:
    with pytest.raises(ContractError, match="exactly one of adata= or path="):
        CellCurator(
            organism="human", tissue="t", run_id="r", output_root=str(tmp_path / "a")
        )
    with pytest.raises(ContractError, match="exactly one of adata= or path="):
        CellCurator(
            adata=_adata(),
            path=tmp_path / "nope.h5ad",
            organism="human",
            tissue="t",
            run_id="r",
            output_root=str(tmp_path / "b"),
        )


def test_a_missing_cluster_key_fails_before_any_run_directory_work(tmp_path: Path) -> None:
    with pytest.raises(ContractError, match="absent from obs"):
        _curator(tmp_path, cluster_key="not_a_column")


def test_reusing_a_run_id_with_a_new_object_is_refused(tmp_path: Path) -> None:
    _curator(tmp_path)
    with pytest.raises(ContractError, match="single-use"):
        _curator(tmp_path)


def test_overrides_win_over_generated_values(tmp_path: Path) -> None:
    curator = _curator(
        tmp_path,
        overrides={"evidence": {"min_confidence": 0.99}, "run": {"mode": "evidence-only"}},
    )
    assert curator.validate().evidence.min_confidence == pytest.approx(0.99)


def test_generated_config_is_plain_yaml_a_reviewer_can_read(tmp_path: Path) -> None:
    curator = _curator(tmp_path)
    value = yaml.safe_load(curator.config_path.read_text())
    assert value["schema_version"] == 3
    assert value["run"]["immutable"] is True
    # Gates are preserved, not silently disabled by the convenience layer.
    assert value["review"]["require_explicit_human_approval"] is True
    assert value["knowledge"]["require_local_grounding_for_labels"] is True


def test_annotate_reaches_a_frozen_review_packet(tmp_path: Path) -> None:
    """The object must be able to finish a run, not just start one.

    The granular methods stop at decide_subclusters(); before annotate() existed
    the only way to a review packet was to drop down to the function API.
    """

    curator = _unattended(tmp_path, mode="annotation")
    packet = curator.annotate()
    assert packet.is_dir()
    assert curator.validate_run()["status"] == "complete"


def test_annotate_refuses_a_multi_level_run_with_a_useful_message(tmp_path: Path) -> None:
    curator = _unattended(
        tmp_path,
        mode="annotation",
        hierarchy_levels=("L1", "L2"),
    )
    with pytest.raises(ContractError, match="execute_hierarchy_manifest"):
        curator.annotate()


def test_the_added_methods_all_delegate_and_exist(tmp_path: Path) -> None:
    """Guard against a method being documented but never wired."""

    curator = _curator(tmp_path)
    for name in (
        "map_cells",
        "freeze_review_packet",
        "annotate",
        "execute_hierarchy_manifest",
        "run_evidence_lane",
        "profile_refined_clusters",
        "assign_labels",
        "confirm_context",
        "list_assumptions",
        "add_assumption",
        "list_capabilities",
        "validate_knowledge",
        "provision_ontology",
        "assemble_markers",
        "build_evidence_cards",
        "build_report",
        "check_report",
        "build_critic_packet",
        "import_critic_findings",
        "import_critic_reconciliation",
        "validate_critic_reconciliation",
    ):
        assert callable(getattr(curator, name)), name


def test_capability_and_knowledge_inspection_work_without_running_anything(
    tmp_path: Path,
) -> None:
    curator = _curator(tmp_path)
    capabilities = curator.list_capabilities()
    assert capabilities and all(item.assay == "transcriptome" for item in capabilities)

    knowledge = curator.validate_knowledge()
    assert knowledge["status"] == "complete"
    assert sorted(knowledge["labels"]) == ["Type A", "Type B"]
