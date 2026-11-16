"""A measurement that was never taken must survive a trip through disk.

`aggregate_multimodal_evidence` reads `program_evidence.tsv` with
`keep_default_na=False`. An uncallable program has no negative-marker violation,
because the panel gate could not measure one; in memory that is `NaN`, which
pandas writes as an empty field, and the reader gets that field back as `""`.
A bare `float()` then raised `could not convert string to float: ''`.

The existing contract tests assert `isna()` on a frame that never touches disk,
so they pass while the run crashes. The round trip is the only place the two
spellings of absent disagree, which makes it the only place worth testing, and
the first real benchmark run on a 2,000-gene panel is what found this: three of
fourteen PBMC marker programs fell below the panel gate, and the run died on the
first of them.

Two failures are pinned here. The crash, and a quieter one underneath it: the
label summary reduced `positive_coverage`, `negative_violation`, and
`program_coverage` with pandas, on columns that are text at that point, so the
maximum violation of one program could be reported as another program's value.
That one produces a number rather than a traceback.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cell_curator.contracts import measured_float


def _capability() -> object:
    from cell_curator.models import AssayKind, Capability, CapabilityLevel

    return Capability(
        assay="transcriptome",
        representation="logcounts",
        assay_kind=AssayKind("rna"),
        feature_space="gene",
        level=CapabilityLevel.ANNOTATION,
        claims=("cell_identity",),
        backend="scipy_mannwhitney",
    )


def _config(tmp_path):
    import yaml

    from cell_curator.config import load_config

    template = (
        Path(__file__).resolve().parents[1]
        / "skills"
        / "cell-curator"
        / "assets"
        / "config.template.yaml"
    )
    value = yaml.safe_load(template.read_text())
    value["run"].update(
        {
            "run_id": "unmeasured-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "annotation",
        }
    )
    value["input"]["path"] = str(tmp_path / "source.h5ad")
    markers = tmp_path / "markers.json"
    markers.write_text('{"Type A": {"pos": ["A1"]}, "Type B": {"pos": ["B1"]}}')
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "round-trip-v1"}
    ]
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(value, sort_keys=False))
    return load_config(path)


def _row(label: str, *, coverage: object, violation: object, score: float = 1.0) -> dict:
    return {
        "schema_version": 3,
        "cluster_id": "0",
        "assay": "transcriptome",
        "representation": "logcounts",
        "lane_key": "transcriptome__logcounts",
        "label": label,
        "cl_id": "",
        "parent_label": "",
        "axis": "cell_type",
        "positive_coverage": coverage,
        "negative_violation": violation,
        "program_coverage": coverage,
        "score": score,
        "detection_valid": "true",
    }


def _round_tripped(rows: list[dict[str, object]], tmp_path) -> pd.DataFrame:
    """Write the lane table to disk and read it the way the pipeline reads it.

    The round trip is the point of this file. Constructing the frame in memory
    reproduces the existing contract tests, and those pass while a real run
    crashes, because `keep_default_na=False` turns the `NaN` that means "not
    measured" into the empty string that `float()` cannot parse.
    """

    path = tmp_path / "program_evidence.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def test_an_unmeasured_value_is_absent_rather_than_an_error() -> None:
    for absent in ("", "   ", "nan", "NaN", "None", "null", "NA"):
        assert measured_float(absent) is None, absent
    assert measured_float("0.75") == pytest.approx(0.75)
    assert measured_float(0.5) == pytest.approx(0.5)
    assert measured_float(float("nan")) is None
    assert measured_float(float("inf")) is None
    assert measured_float("not a number") is None


