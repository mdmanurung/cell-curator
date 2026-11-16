"""Approved call rules (decision rows 7, 8, 14) in plans/annotation-corpus-decisions.md."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cell_curator.annotation import _classifier_counts, _marker_family, _signature_support
from cell_curator.evidence import benjamini_hochberg


def _frame(role: str, pvalues: list[float], effects: list[float]) -> pd.DataFrame:
    n = len(pvalues)
    return pd.DataFrame(
        {
            "cluster_id": "c1",
            "comparator_role": role,
            "method": "wilcoxon",
            "feature": [f"{role}_g{i}" for i in range(n)],
            "effect": effects,
            "pvalue": pvalues,
            "adjusted_pvalue": np.nan,
        }
    )


def test_one_bh_family_spans_every_comparator_in_the_scope() -> None:
    frames = [_frame("parent_remainder", [0.01, 0.5], [1.0, 2.0]), _frame("closest_sibling", [0.02, 0.9], [3.0, 4.0])]
    family = _marker_family(frames, adjusted_pvalue_max=0.05)
    expected = benjamini_hochberg(np.array([0.01, 0.5, 0.02, 0.9]))
    assert expected.tolist() == [0.04, 0.6666666666666666, 0.04, 0.9]
    assert set(family["feature"]) == {"parent_remainder_g0", "closest_sibling_g0"}


def test_top_genes_are_ranked_by_effect_after_the_p_filter() -> None:
    frame = _frame("parent_remainder", [1e-9] * 5, [0.1, 0.5, 0.3, 0.9, 0.7])
    family = _marker_family([frame], adjusted_pvalue_max=0.05, top_n=3)
    assert family["feature"].tolist() == ["parent_remainder_g3", "parent_remainder_g4", "parent_remainder_g1"]
    assert family["rank"].tolist() == [1, 2, 3]


def test_coarser_classifier_label_abstains_and_unmapped_label_is_dropped() -> None:
    # Decision rows 12 and 13.
    parents = {"CD8 T": "T", "T": "ROOT", "B": "ROOT"}
    canonical = {"CD8 T", "T", "B"}
    raw = pd.Series(["CD8 T", "T", "B", "raw-x", "raw-y"])
    counts = _classifier_counts(raw, {"raw-y": "B"}, canonical, parents, "CD8 T")
    assert counts == {"agree": 1, "disagree": 2, "coarser": 1, "unmapped": 1}


def test_signature_support_uses_the_state_rule() -> None:
    from types import SimpleNamespace

    config = SimpleNamespace(
        evidence=SimpleNamespace(min_positive_coverage=0.5),
        annotation=SimpleNamespace(min_state_positive_fraction=0.5),
    )
    detection = np.array([[1, 1, 0], [1, 0, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
    groups = pd.Series(["a", "a", "b", "b"])
    lookup = {"g0": 0, "g1": 1, "g2": 2}
    support, measured, coverage = _signature_support(detection, groups, lookup, ("g0", "g1"), ("g2",), config)
    assert measured and coverage == 1.0
    # cluster a: both cells detect >= half the positives and no g2; cluster b: neither passes.
    assert support == {"a": 1.0, "b": 0.0}
    _, measured, _ = _signature_support(detection, groups, lookup, ("nope",), (), config)
    assert not measured
