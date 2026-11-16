from pathlib import Path

import pandas as pd
import pytest
from consensus_fixture import write_consensus_run
from test_annotation_pipeline import annotation_fixture

from cell_curator.annotation import collect_level, finalize_annotation
from cell_curator.consensus_import import import_consensus_level
from cell_curator.contracts import ContractError
from cell_curator.guided import approve_assumptions


def test_identity_approvals_do_not_approve_a_study_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch, approved=False)
    vocabulary = pd.read_csv(root / "plan" / "vocabulary.tsv", sep="\t")
    # Only the configured level's identity scope: an L1-only run never reaches the L2 identity scopes.
    identity = vocabulary[vocabulary["axis"].eq("cell_type") & vocabulary["level"].eq("L1")]
    for level, parent in identity[["level", "parent_scope"]].drop_duplicates().itertuples(index=False):
        approve_assumptions(config_path, reviewer="synthetic-reviewer", scope=f"{level}::{parent}")
    collect_level(config_path, "L1")
    write_consensus_run(
        config_path, root, {cluster: values[0] for cluster, values in expected.items()}
    )
    import_consensus_level(config_path, "L1")
    with pytest.raises(ContractError, match="Gate A|approved assumption"):
        finalize_annotation(config_path)
    assert not (root / "annotation" / "cell_labels.tsv").exists()
    state = vocabulary[vocabulary["axis"].eq("cell_state")].iloc[0]
    approve_assumptions(config_path, reviewer="synthetic-reviewer", scope=f"{state['level']}::{state['parent_scope']}")
    finalize_annotation(config_path)
    cells = pd.read_csv(root / "annotation" / "cell_labels.tsv", sep="\t")
    assert cells["L4"].str.endswith("S+").sum() == 20
