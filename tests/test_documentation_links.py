"""Every relative link in tracked Markdown must resolve from its own directory.

`quadbio-parity.md` shipped 134 links written one directory level too shallow,
plus citations of three test files deleted in the first packaging commit, and
nothing noticed for the life of the repository. Markdown has no compiler; this
is the compiler.
"""

from __future__ import annotations

import re
import subprocess

import pytest
from conftest import REPO_ROOT

LINK = re.compile(r"\]\((?!https?:|mailto:)([^)]+)\)")


def _tracked_markdown() -> list[str]:
    listed = subprocess.run(
        ["git", "ls-files", "*.md"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return sorted(path for path in listed.stdout.split() if path)


@pytest.mark.parametrize("relative_path", _tracked_markdown())
def test_relative_links_resolve(relative_path: str) -> None:
    document = REPO_ROOT / relative_path
    broken = []
    for target in LINK.findall(document.read_text()):
        # A fragment or query is addressing within a target, not a separate file.
        path = target.split("#", 1)[0].split("?", 1)[0].strip()
        if not path:
            continue
        if not (document.parent / path).resolve().exists():
            broken.append(target)
    assert not broken, f"{relative_path} links to missing paths: {broken}"
