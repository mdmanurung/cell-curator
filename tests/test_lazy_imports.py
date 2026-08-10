"""Importing the package must not drag in the scientific stack.

`api.py` promises callers can inspect the version and configuration models
"without importing the scientific stack". That promise is easy to break by
accident — one module-level `import pandas` anywhere in the import graph of
`cell_curator/__init__.py` is enough, and nothing else in the suite would
notice.

The check runs in a fresh interpreter on purpose: by the time pytest has
collected this file, numpy and pandas are already in `sys.modules`, so an
in-process assertion would pass no matter what the package does.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

HEAVY = ("numpy", "pandas", "anndata", "scanpy", "mudata", "torch", "scipy")


def test_bare_import_stays_out_of_the_scientific_stack() -> None:
    probe = textwrap.dedent(
        f"""
        import sys

        import cell_curator

        assert cell_curator.__version__
        heavy = {HEAVY!r}
        loaded = sorted(set(heavy) & {{name.split(".")[0] for name in sys.modules}})
        print(",".join(loaded))
        """
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    loaded = [name for name in completed.stdout.strip().split(",") if name]
    assert not loaded, (
        f"`import cell_curator` pulled in {loaded}; keep heavy imports inside "
        "function bodies so the package stays cheap to import"
    )


def test_the_facade_is_still_importable_and_constructible() -> None:
    """Laziness must not be achieved by breaking the public surface."""

    probe = textwrap.dedent(
        """
        from cell_curator import CellCurator

        assert callable(CellCurator)
        print("ok")
        """
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "ok"
