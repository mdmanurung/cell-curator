from __future__ import annotations

from pathlib import Path

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
from cell_curator.evidence import compute_bottom_up_markers, score_marker_programs
from cell_curator.hierarchy import identity_vocabulary, validate_parent_vocabulary
from cell_curator.models import MarkerProgram
from cell_curator.provenance import atomic_publish


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






def test_a_partly_measured_negative_list_is_judged_on_its_measured_genes(config: dict) -> None:
    """Protocol v5 abstained everywhere: one unmeasured negative voided the whole list."""

    model = _model(config)
    capabilities = build_capability_registry(model)
    # Columns: positive, measured negative. Cluster 0 expresses the negative.
    matrix = np.array([[5.0, 3.0], [4.0, 2.0], [0.0, 0.0], [0.0, 0.0]])
    groups = pd.Series(["0", "0", "1", "1"])
    program = MarkerProgram(
        label="Target",
        positive=("measured",),
        negative=("measured-negative", "unmeasured-negative"),
    )
    common = dict(
        detection=matrix,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    features = pd.Index(["measured", "measured-negative"])
    judged = score_marker_programs(
        matrix, groups, features, [program], min_negative_coverage=0.5, **common
    )
    assert judged.set_index("cluster_id")["negative_violation"].to_dict() == {"0": 1.0, "1": 0.0}
    strict = score_marker_programs(
        matrix, groups, features, [program], min_negative_coverage=0.75, **common
    )
    assert strict["negative_violation"].isna().all()




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
    assert list(feature_ids_for(query, "gene_symbol", declaration={"source": "raw"})) == [
        "raw-a",
        "raw-b",
        "raw-c",
    ]
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






def test_a_gene_cannot_be_both_expected_and_forbidden() -> None:
    from pydantic import ValidationError

    from cell_curator.knowledge import collect_programs

    with pytest.raises(ValidationError, match="NKG7"):
        MarkerProgram(label="CD8 T", positive=("CD8A", "NKG7"), negative=("NKG7",))

    class Source:
        def __init__(self, program: MarkerProgram) -> None:
            self.program = program

        def programs(self) -> list[MarkerProgram]:
            return [self.program]

    with pytest.raises(ValidationError, match="NKG7"):
        collect_programs(
            [
                Source(MarkerProgram(label="CD8 T", positive=("NKG7",))),
                Source(MarkerProgram(label="CD8 T", negative=("NKG7",))),
            ]
        )


def test_a_two_group_comparison_ranks_programs_by_effect_not_by_name(config: dict) -> None:
    """Two group means z-score to exactly +/-1, so a 3-unit and a 0.25-unit shift
    used to tie and 'Alpha' won alphabetically."""

    model = _model(config)
    capabilities = build_capability_registry(model)
    rng = np.random.default_rng(0)
    candidate = np.column_stack([rng.normal(3.0, 0.3, 50), rng.normal(1.25, 0.3, 50)])
    reference = np.column_stack([rng.normal(0.0, 0.3, 50), rng.normal(1.0, 0.3, 50)])
    matrix = np.clip(np.vstack([candidate, reference]), 0, None)
    evidence = score_marker_programs(
        matrix,
        pd.Series(["candidate"] * 50 + ["reference"] * 50),
        pd.Index(["ZE1", "AL1"]),
        [
            MarkerProgram(label="Alpha subtype", positive=("AL1",)),
            MarkerProgram(label="Zeta subtype", positive=("ZE1",)),
        ],
        detection=matrix,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    target = evidence[evidence["cluster_id"].eq("candidate")].sort_values(
        ["score", "label"], ascending=[False, True]
    )
    assert target.iloc[0]["label"] == "Zeta subtype"
    assert target.iloc[0]["score"] > target.iloc[1]["score"]
    assert target.iloc[0]["vote_fraction"] > 0.5


CL_SNAPSHOT = Path(__file__).parents[1] / "benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cl.obo"




def test_cells_without_marker_signal_do_not_vote(config: dict) -> None:
    """Half the cluster expresses nothing; those cells must not be handed to a program."""

    model = _model(config)
    capabilities = build_capability_registry(model)
    features = pd.Index(["CD3E", "MS4A1", "LYZ"])
    t, b, m = (np.eye(3)[i] * 5.0 for i in range(3))
    matrix = np.vstack([t] * 50 + [np.zeros(3)] * 50 + [b] * 100 + [m] * 100)
    evidence = score_marker_programs(
        matrix,
        pd.Series(["0"] * 100 + ["1"] * 100 + ["2"] * 100),
        features,
        [
            MarkerProgram(label="B", positive=("MS4A1",)),
            MarkerProgram(label="Mono", positive=("LYZ",)),
            MarkerProgram(label="T", positive=("CD3E",)),
        ],
        detection=matrix,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    cluster = evidence[evidence["cluster_id"].eq("0")].set_index("label")
    assert cluster.loc["T", "vote_fraction"] == pytest.approx(0.5)
    assert cluster.loc["B", "n_voters"] == 0
    assert cluster.loc["Mono", "n_voters"] == 0


def test_backoff_never_leaves_the_cell_ontology() -> None:
    """AIDA v3 cluster 14: HSC, T cell, and ILC programs tied, and their deepest
    shared ancestor in the snapshot was UBERON:0000061 "anatomical structure"."""

    from cell_curator.ontology import CellOntology, backoff_term

    ontology = CellOntology.load(CL_SNAPSHOT)
    tied = [
        ("CL:0000624", 0.30),
        ("CL:0000037", 0.25),
        ("CL:0000625", 0.22),
        ("CL:0001065", 0.20),
    ]
    assert backoff_term(ontology, tied, min_margin=0.2, min_confidence=0.6) is None


def test_a_blank_cell_in_a_marker_table_is_absent_not_the_string_nan(tmp_path: Path) -> None:
    from cell_curator.knowledge import CachedTableProvider

    table = tmp_path / "markers.tsv"
    table.write_text(
        "label\tpositive\tnegative\tparent\tcl_id\n"
        "T cell\tCD3E;CD3D\t\t\t\n"
        "B cell\tMS4A1\tCD3E\t\tCL:0000236\n"
    )
    programs = {item.label: item for item in CachedTableProvider(table, version="x").programs()}
    assert programs["T cell"].parent is None and programs["T cell"].cl_id is None
    assert programs["T cell"].negative == ()
    assert programs["B cell"].cl_id == "CL:0000236"
