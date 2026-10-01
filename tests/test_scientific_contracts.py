from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from cell_curator.capabilities import build_capability_registry
from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import (
    ContractError,
    require_safe_path_component,
    safe_path_component,
)
from cell_curator.data import feature_ids_for
from cell_curator.evidence import (
    aggregate_multimodal_evidence,
    compute_bottom_up_markers,
    compute_gpu_bottom_up_markers,
    score_marker_programs,
)
from cell_curator.hierarchy import identity_vocabulary, validate_parent_vocabulary
from cell_curator.models import MarkerProgram
from cell_curator.provenance import atomic_publish
from cell_curator.refinement import _stored_evaluation


def _model(config: dict) -> CellCuratorConfig:
    value = yaml.safe_load(yaml.safe_dump(config))
    value["run"]["mode"] = "annotation"
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": "markers.json", "version": "fixture"}
    ]
    value["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "scipy_mannwhitney"
    value["markers"].update(
        {
            "rapids_singlecell_version": "",
            "rapids_distribution": "",
            "cuda_major": None,
            "allowed_gpu_contracts": [],
        }
    )
    return CellCuratorConfig.model_validate(value)


def test_bottom_up_bh_correction_covers_all_features(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pvalues = iter([0.01, *([0.5] * 99)] * 2)

    class Result:
        def __init__(self, pvalue: float) -> None:
            self.pvalue = pvalue

    monkeypatch.setattr(
        "cell_curator.evidence.stats.mannwhitneyu",
        lambda *_args, **_kwargs: Result(next(pvalues)),
    )
    matrix = np.zeros((8, 100), dtype=float)
    matrix[:4, 0] = 10
    result = compute_bottom_up_markers(
        matrix,
        pd.Series(["target"] * 4 + ["reference"] * 4),
        pd.Index([f"g{index:03d}" for index in range(100)]),
        detection=matrix,
        top_n=1,
    )
    target = result[result["cluster_id"].eq("target")].iloc[0]
    assert target["feature"] == "g000"
    assert target["pvalue"] == pytest.approx(0.01)
    assert target["adjusted_pvalue"] == pytest.approx(0.5)


def test_missing_detection_and_unmeasured_program_features_force_abstention(
    config: dict,
) -> None:
    model = _model(config)
    capabilities = build_capability_registry(model)
    matrix = np.array([[5.0], [4.0], [0.0], [0.0]])
    groups = pd.Series(["0", "0", "1", "1"])
    program = MarkerProgram(
        label="Target",
        positive=("measured", *tuple(f"missing-{index}" for index in range(99))),
    )
    evidence = score_marker_programs(
        matrix,
        groups,
        pd.Index(["measured"]),
        [program],
        detection=None,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    assert evidence["detection_valid"].eq(False).all()
    assert evidence["positive_coverage"].isna().all()
    assert evidence["negative_violation"].isna().all()
    assert evidence["program_coverage"].eq(0.01).all()
    _, decisions = aggregate_multimodal_evidence(
        evidence,
        capabilities,
        model,
        scope="global",
        cluster_ids=["0", "1"],
    )
    assert all(item.abstained and item.cell_type == "Unknown" for item in decisions)


def test_unmeasured_expected_negative_forces_abstention(config: dict) -> None:
    model = _model(config)
    capabilities = build_capability_registry(model)
    matrix = np.array([[5.0], [4.0], [0.0], [0.0]])
    groups = pd.Series(["0", "0", "1", "1"])
    program = MarkerProgram(
        label="Target",
        positive=("measured",),
        negative=("unmeasured-negative",),
    )
    evidence = score_marker_programs(
        matrix,
        groups,
        pd.Index(["measured"]),
        [program],
        detection=matrix,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    assert evidence["negative_violation"].isna().all()
    _, decisions = aggregate_multimodal_evidence(
        evidence,
        capabilities,
        model,
        scope="global",
        cluster_ids=["0", "1"],
    )
    assert all(item.abstained and item.cell_type == "Unknown" for item in decisions)


def test_axisless_evidence_is_refused_rather_than_called_cell_type(config: dict) -> None:
    """Losing the axis column must stop the run, not relabel state as identity.

    The aggregator used to default a missing axis to `cell_type`, so an evidence
    frame that had been round-tripped through a table without the column — a
    cached export, a hand-edited TSV — turned every cycling, stress, and
    low-quality program into a candidate identity without saying so.
    """

    model = _model(config)
    capabilities = build_capability_registry(model)
    matrix = np.array([[5.0], [4.0], [0.0], [0.0]])
    groups = pd.Series(["0", "0", "1", "1"])
    program = MarkerProgram(label="Cycling", positive=("measured",), axis="cell_state")
    evidence = score_marker_programs(
        matrix,
        groups,
        pd.Index(["measured"]),
        [program],
        detection=matrix,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    assert evidence["axis"].eq("cell_state").all()

    with pytest.raises(ContractError, match="axis"):
        aggregate_multimodal_evidence(
            evidence.drop(columns=["axis"]),
            capabilities,
            model,
            scope="global",
            cluster_ids=["0", "1"],
        )


def _vocabulary_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"level": "L1", "parent_scope": "ROOT", "label": "Type A", "axis": "cell_type"},
            {"level": "L1", "parent_scope": "ROOT", "label": "Cycling", "axis": "cell_state"},
        ]
    )


def test_a_state_program_is_refused_as_an_identity_candidate() -> None:
    """An operator-named vocabulary asserts its rows are identity candidates.

    `validate_parent_vocabulary` checked only level, parent, and label, so a
    cell state program sitting in the vocabulary was assignable as a cell's
    type. Nothing downstream re-checked the axis.
    """

    with pytest.raises(ContractError, match="cell_type"):
        validate_parent_vocabulary(_vocabulary_rows())

    identity_only = _vocabulary_rows().iloc[[0]].reset_index(drop=True)
    validate_parent_vocabulary(identity_only)


def test_a_vocabulary_without_an_axis_column_stays_legal() -> None:
    """The documented minimal contract is level, parent, and label.

    A frame without the column carries no state rows to admit, so refusing it
    would break the contract without closing anything.
    """

    bare = _vocabulary_rows().drop(columns=["axis"]).iloc[[0]].reset_index(drop=True)
    validate_parent_vocabulary(bare)
    filtered, excluded = identity_vocabulary(bare)
    assert excluded == []
    assert filtered.equals(bare)


def test_filtering_refuses_a_type_program_descended_from_a_state_program() -> None:
    """`MarkerProgram.parent` is not constrained to one axis, so this is reachable.

    Dropping the state rows would leave the survivor pointing at a parent that
    no longer exists, which silently produces a broken hierarchy instead of a
    stated error.
    """

    crossed = pd.concat(
        [
            _vocabulary_rows(),
            pd.DataFrame(
                [
                    {
                        "level": "L2",
                        "parent_scope": "Cycling",
                        "label": "Type B",
                        "axis": "cell_type",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    with pytest.raises(ContractError, match="descends from a cell state"):
        identity_vocabulary(crossed)

    filtered, excluded = identity_vocabulary(_vocabulary_rows())
    assert excluded == ["Cycling"]
    assert filtered["label"].tolist() == ["Type A"]


def test_config_rejects_missing_comparator_and_duplicate_lane_pair(config: dict) -> None:
    missing_comparator = yaml.safe_load(yaml.safe_dump(config))
    missing_comparator["comparators"]["roles"] = ["parent_remainder"]
    with pytest.raises(ValueError, match="closest_sibling"):
        CellCuratorConfig.model_validate(missing_comparator)

    duplicated = yaml.safe_load(yaml.safe_dump(config))
    duplicated["evidence"]["lanes"].append(
        {
            "assay": "transcriptome",
            "representation": "logcounts",
            "lane_key": "duplicate-key",
        }
    )
    with pytest.raises(ValueError, match="pair is duplicated"):
        CellCuratorConfig.model_validate(duplicated)


def test_unsafe_parent_identifiers_have_distinct_reversible_path_components() -> None:
    left = safe_path_component("parent/a")
    right = safe_path_component("parent?a")
    assert left != right
    assert safe_path_component("/") != safe_path_component("u-Lw")
    assert left != safe_path_component("u-cGFyZW50L2E")
    assert "/" not in left
    assert "?" not in right
    assert safe_path_component("already-safe") == "already-safe"
    with pytest.raises(ContractError, match="safe path component"):
        require_safe_path_component("../../outside", field="hierarchy level")


def test_raw_representation_uses_raw_feature_axis() -> None:
    query = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(
            {"gene_symbol": ["current-a", "current-b"]},
            index=["current-a", "current-b"],
        ),
    )
    query.raw = ad.AnnData(
        X=np.ones((2, 3)),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(
            {"gene_symbol": ["raw-a", "raw-b", "raw-c"]},
            index=["raw-a", "raw-b", "raw-c"],
        ),
    )
    assert list(
        feature_ids_for(query, "gene_symbol", declaration={"source": "raw"})
    ) == ["raw-a", "raw-b", "raw-c"]
    assert list(feature_ids_for(query, "gene_symbol", declaration={"source": "X"})) == [
        "current-a",
        "current-b",
    ]


def test_evidence_feature_ids_apply_declared_gene_mapping(tmp_path) -> None:
    query = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["ENSG000001.1", "ENSG000002.4"]),
    )
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("ensembl_id\tgene_symbol\nENSG000001\tA\nENSG000002\tB\n")
    assert list(
        feature_ids_for(
            query,
            "",
            declaration={"source": "X"},
            gene_mapping_path=mapping,
        )
    ) == ["A", "B"]


def test_gpu_marker_lane_rejects_one_group_before_claiming_methods() -> None:
    with pytest.raises(ContractError, match="no marker methods were executed"):
        compute_gpu_bottom_up_markers(
            np.ones((2, 1)),
            pd.Series(["only", "only"]),
            pd.Index(["gene"]),
            detection=None,
            rsc=object(),
            marker_config={},
        )


def test_stored_resolution_excludes_transient_children(config: dict) -> None:
    value = yaml.safe_load(yaml.safe_dump(config))
    value["input"]["keys"]["stored_resolutions"] = ["stored_low", "stored_high"]
    model = _model(value)
    parent_obs = pd.DataFrame(
        {
            "stored_low": ["A"] * 6 + ["B"] * 6,
            "stored_high": ["A"] * 5 + ["transient"] + ["B"] * 6,
        },
        index=[f"cell-{index}" for index in range(12)],
    )
    embedding = np.vstack(
        [np.column_stack((np.arange(6), np.zeros(6))),
         np.column_stack((np.arange(6), np.ones(6) * 10))]
    )
    evaluation = _stored_evaluation(parent_obs, embedding, model)
    assert evaluation is not None
    assert evaluation.n_clusters == 2
    assert evaluation.eligible_children is not None
    selected_names = {
        str(pd.Categorical(parent_obs["stored_high"]).categories[index])
        for index in evaluation.eligible_children
    }
    assert selected_names == {"A", "B"}


def test_immutable_atomic_publisher_resumes_identical_and_refuses_overwrite(
    tmp_path,
) -> None:
    path = tmp_path / "artifact.txt"
    stale = tmp_path / ".artifact.txt.partial.previous-job"
    stale.write_text("incomplete\n")
    first = atomic_publish(path, b"complete\n")
    assert atomic_publish(path, b"complete\n") == first
    with pytest.raises(ContractError, match="non-identical"):
        atomic_publish(path, b"replacement\n")
    assert path.read_text() == "complete\n"
    assert stale.read_text() == "incomplete\n"


def test_confidence_is_the_cluster_share_voting_for_the_label(config: dict) -> None:
    """A cluster that is 80% T cells is called T at confidence 0.8, whatever else is in the
    vocabulary, and the 20% B cells in it do not count against T's expected negatives."""

    model = _model(config)
    capabilities = build_capability_registry(model)
    features = pd.Index(["CD3E", "MS4A1", "LYZ", *[f"X{i}" for i in range(5)]])
    t, b, m = (np.eye(len(features))[i] * 5.0 for i in range(3))
    matrix = np.vstack([t] * 80 + [b] * 20 + [b] * 100 + [m] * 100)
    groups = pd.Series(["0"] * 100 + ["1"] * 100 + ["2"] * 100)
    core = [
        MarkerProgram(label="T", positive=("CD3E",), negative=("MS4A1", "LYZ")),
        MarkerProgram(label="B", positive=("MS4A1",), negative=("CD3E", "LYZ")),
        MarkerProgram(label="Mono", positive=("LYZ",), negative=("CD3E", "MS4A1")),
    ]
    padding = [MarkerProgram(label=f"Pad{i}", positive=(f"X{i}",)) for i in range(5)]

    def call(programs: list[MarkerProgram]):
        evidence = score_marker_programs(
            matrix,
            groups,
            features,
            programs,
            detection=matrix,
            run_id=model.run.run_id,
            assay="transcriptome",
            representation="logcounts",
            capability=capabilities[0],
        )
        _, decisions = aggregate_multimodal_evidence(
            evidence, capabilities, model, scope="global", cluster_ids=["0", "1", "2"]
        )
        return evidence, {item.cluster_id: item for item in decisions}

    evidence, decisions = call(core)
    t_in_0 = evidence[evidence["cluster_id"].eq("0") & evidence["label"].eq("T")].iloc[0]
    assert t_in_0["vote_fraction"] == pytest.approx(0.8)
    assert t_in_0["n_voters"] == 80
    assert t_in_0["negative_violation"] == 0.0
    assert [decisions[c].cell_type for c in "012"] == ["T", "B", "Mono"]
    assert decisions["0"].confidence == pytest.approx(0.8)
    assert decisions["0"].margin == pytest.approx(0.6)

    _, padded = call(core + padding)
    assert padded["0"].cell_type == "T"
    assert padded["0"].confidence == pytest.approx(0.8)


def test_an_exact_tie_reports_bounded_entropy() -> None:
    from cell_curator.evidence import _normalized_entropy

    # Five-way tie: unclipped, rounding gives 1.0000000000000002, which the model refuses.
    assert _normalized_entropy(np.full(5, 1 / 5)) <= 1.0
