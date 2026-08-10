"""Library logging must be silent by default and useful when switched on.

A long phase used to produce no output at all, so a notebook cell just sat there
for minutes with no way to tell progress from a hang. Phase-boundary logging
fixes that, but it must not turn into stdout noise for CLI users or for anyone
embedding the package: the library configures no handlers and sets no levels,
which leaves it quiet until a caller opts in.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
from test_session_facade import _unattended


def test_no_records_escape_at_the_default_level(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    curator = _unattended(tmp_path, mode="annotation")
    with caplog.at_level(logging.WARNING, logger="cell_curator"):
        curator.annotate()
    assert [record for record in caplog.records if record.levelno >= logging.WARNING] == []


def test_enabling_info_reveals_the_phase_boundaries(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The opt-in must actually produce something worth reading."""

    curator = _unattended(tmp_path, mode="annotation", run_id="logging-v1")
    with caplog.at_level(logging.INFO, logger="cell_curator"):
        curator.annotate()

    messages = [record.getMessage() for record in caplog.records]
    assert any("run_annotation starting" in message for message in messages)
    for phase in ("auditing parent", "computing evidence", "validating the run"):
        assert any(phase in message for message in messages), phase
    assert any("evidence lane starting" in message for message in messages)


def test_the_library_installs_no_handlers_and_forces_no_level() -> None:
    """Configuring logging is the application's job, not the library's."""

    for name in (
        "cell_curator",
        "cell_curator.pipeline",
        "cell_curator.evidence",
        "cell_curator.markers",
    ):
        logger = logging.getLogger(name)
        assert logger.handlers == [], f"{name} attached a handler"
        assert logger.level == logging.NOTSET, f"{name} forced a level"


def test_nothing_is_written_to_stdout_while_running(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    curator = _unattended(tmp_path, mode="annotation", run_id="stdout-v1")
    curator.annotate()
    captured = capsys.readouterr()
    assert captured.out == "", f"unexpected stdout: {captured.out!r}"
