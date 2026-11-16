"""Library logging must be silent by default and useful when switched on.

A long phase used to produce no output at all, so a notebook cell just sat there
for minutes with no way to tell progress from a hang. Phase-boundary logging
fixes that, but it must not turn into stdout noise for CLI users or for anyone
embedding the package: the library configures no handlers and sets no levels,
which leaves it quiet until a caller opts in.
"""

from __future__ import annotations

import logging


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
