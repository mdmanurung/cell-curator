"""The whole workflow must run as ordinary Python, with no shell in the loop.

Most of this package's capabilities were once reachable only by invoking the CLI,
which makes it unusable from a script or a notebook without `subprocess`. This
test is the standing proof that the situation does not return: it drives a
complete run — freeze, evidence, decide, map, review, validate, write-back —
entirely through the Python surface, with three guards that would catch a
regression to shelling out, to killing the interpreter, or to printing.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from test_session_facade import _adata, _unattended

SUBPROCESS_ENTRY_POINTS = (
    "Popen",
    "run",
    "call",
    "check_call",
    "check_output",
    "getoutput",
    "getstatusoutput",
)


@pytest.fixture
def no_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail loudly if any library call tries to shell out."""

    def forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "cell_curator spawned a subprocess; every capability must be callable "
            f"in-process. Attempted: {args!r}"
        )

    for name in SUBPROCESS_ENTRY_POINTS:
        if hasattr(subprocess, name):
            monkeypatch.setattr(subprocess, name, forbidden)


def test_a_complete_run_needs_no_subprocess_and_prints_nothing(
    tmp_path: Path,
    no_subprocess: None,
    capsys: pytest.CaptureFixture[str],
) -> None:
    curator = _unattended(tmp_path, mode="annotation", run_id="purity-v1")

    # SystemExit is a BaseException: uncaught, it can abort the whole pytest
    # session rather than failing one test, so catch it explicitly.
    try:
        curator.inspect()
        curator.freeze()
        packet = curator.annotate()
        report = curator.check_report()
        validation = curator.validate_run()
    except SystemExit as exc:  # pragma: no cover - guard, not expected to fire
        pytest.fail(f"SystemExit escaped the pure-Python workflow: {exc}")

    assert packet.is_dir()
    assert report["status"]
    assert validation["status"] == "complete"

    captured = capsys.readouterr()
    assert captured.out == "", f"library code wrote to stdout: {captured.out!r}"


def test_the_granular_phases_are_individually_callable(
    tmp_path: Path, no_subprocess: None
) -> None:
    """Each step must stand alone, so a notebook can stop and look between them."""

    curator = _unattended(tmp_path, mode="annotation", run_id="granular-v1")
    curator.prepare()

    assert not curator.audit_parents().empty
    curator.propose_subclusters()
    evidence, bottom_up = curator.evidence()
    assert not evidence.empty and not bottom_up.empty
    decisions, parents = curator.decide_subclusters()
    assert decisions is not None and parents is not None
    assert not curator.map_cells().empty
    assert curator.review_packet().exists()
    assert curator.validate_run()["status"] == "complete"


def test_the_function_api_reaches_the_same_result_without_the_object(
    tmp_path: Path, no_subprocess: None
) -> None:
    """The object is a convenience; the plain functions must remain sufficient."""

    from cell_curator import annotate, load_config, validate_run

    curator = _unattended(tmp_path, mode="annotation", run_id="functions-v1")
    config_path = curator.config_path

    packet = annotate(config_path)
    assert packet.is_dir()
    assert validate_run(config_path)["status"] == "complete"
    assert load_config(config_path).run.run_id == "functions-v1"


def test_read_only_inspection_never_touches_the_caller_object(tmp_path: Path) -> None:
    """A notebook keeps using its AnnData after handing it over."""

    adata = _adata()
    before_obs = list(adata.obs.columns)
    before_shape = adata.shape

    curator = _unattended(tmp_path, adata=adata, run_id="untouched-v1")
    curator.inspect()
    curator.freeze()

    assert list(adata.obs.columns) == before_obs
    assert adata.shape == before_shape
    assert curator.run_directory.is_dir()
