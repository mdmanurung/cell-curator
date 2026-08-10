"""The example notebook must actually run, and must not carry stored outputs.

A quickstart that has drifted out of date is worse than none: it teaches an API
that no longer exists. Executing it here means a signature change breaks the
build rather than the next reader.

This test starts a Jupyter kernel, which is a subprocess. That is expected and is
a different claim from `test_full_workflow_purity.py`, which asserts the
*library* never shells out. Here the subprocess is the thing under test.
"""

from __future__ import annotations

from pathlib import Path

import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / "examples" / "quickstart_synthetic.ipynb"


def test_the_example_notebook_exists() -> None:
    assert NOTEBOOK.is_file(), f"missing example notebook: {NOTEBOOK}"


def test_the_notebook_is_committed_without_outputs() -> None:
    """Executed outputs would bloat the repository and rot silently."""

    nbformat = pytest.importorskip("nbformat")
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    offenders = [
        index
        for index, cell in enumerate(notebook.cells)
        if cell.cell_type == "code" and (cell.get("outputs") or cell.get("execution_count"))
    ]
    assert not offenders, (
        f"code cells {offenders} carry stored output; clear them before committing"
    )


def test_the_notebook_executes_end_to_end() -> None:
    nbformat = pytest.importorskip("nbformat")
    pytest.importorskip("nbclient")
    pytest.importorskip("ipykernel")
    from nbclient import NotebookClient

    # Executed on an in-memory copy; the committed file is never rewritten.
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    NotebookClient(notebook, timeout=600, kernel_name="python3").execute()
