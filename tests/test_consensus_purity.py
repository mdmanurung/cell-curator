"""The consensus flow must run as ordinary Python, quietly, with no shell in the loop.

Ported from the earlier full-workflow purity and logging tests, which drove the removed
pipeline. This drives the flow that replaced it, from collect-level to `validate_run`,
through the public functions only, with guards against spawning a subprocess, printing to
stdout and logging at warning level. The run directory is set up first (scene and plan
may inspect the machine); the guards cover everything after Gate A.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pytest
from consensus_fixture import write_consensus_run
from test_annotation_pipeline import annotation_fixture

from cell_curator import (
    check_review_report,
    collect_annotation_level,
    finalize_corpus_annotation,
    import_consensus_annotation_level,
    load_config,
    validate_run,
)
from cell_curator.pipeline import run_root
from cell_curator.review import build_critic_packet, build_evidence_cards, build_html_report

SUBPROCESS_ENTRY_POINTS = (
    "Popen",
    "run",
    "call",
    "check_call",
    "check_output",
    "getoutput",
    "getstatusoutput",
)


def _forbid_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    def forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError(f"cell_curator spawned a subprocess: {args!r}")

    for name in SUBPROCESS_ENTRY_POINTS:
        if hasattr(subprocess, name):
            monkeypatch.setattr(subprocess, name, forbidden)


def test_the_consensus_flow_needs_no_subprocess_and_stays_quiet(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch)
    capsys.readouterr()
    _forbid_subprocess(monkeypatch)

    # SystemExit is a BaseException: uncaught, it can abort the whole pytest session.
    try:
        with caplog.at_level(logging.WARNING, logger="cell_curator"):
            collect_annotation_level(config_path, "L1")
            write_consensus_run(
                config_path, root, {cluster: values[0] for cluster, values in expected.items()}
            )
            import_consensus_annotation_level(config_path, "L1")
            finalize_corpus_annotation(config_path)
            config = load_config(config_path)
            build_evidence_cards(config, run_root(config))
            build_critic_packet(config, run_root(config))
            build_html_report(config, run_root(config), strict=True)
            report = check_review_report(config_path)
            validation = validate_run(config_path)
    except SystemExit as exc:  # pragma: no cover - guard, not expected to fire
        pytest.fail(f"SystemExit escaped the pure-Python flow: {exc}")

    assert report["status"] == "complete"
    assert validation["status"] == "complete"
    assert capsys.readouterr().out == "", "library code wrote to stdout"
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []
