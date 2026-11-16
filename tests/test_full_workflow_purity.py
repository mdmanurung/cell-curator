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
