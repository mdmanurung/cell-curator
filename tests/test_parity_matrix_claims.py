"""The parity matrix is a claim about the code, so the code gets to check it.

`quadbio-parity.md` names an implementation file, a CLI or Python entry point,
and a test file for each of forty-five capabilities. Every one of those is a
factual assertion that a symbol exists, a command resolves, and a path is real,
and all three were maintained by hand. `test_documentation_links.py` proves the
*paths* resolve; nothing proved the *symbols* are still there, so a renamed
function or a retired subcommand would leave a document that reads as verified
and is not.

The sibling project generates its method-selection document from the router and
fails on drift, which is why this file exists rather than another prose review.
It is a floor and not a substitute for review: a matrix can pass this and still
be wrong about behaviour. What it prevents is the specific failure of a
document accumulating confident references to code that has moved on.

The other tables in the same document are prose about other people's projects and
are deliberately not checked here. Verifying a claim about a third-party
repository would mean either vendoring it or asserting its contents from a
network call, and a test that silently depends on a GitHub round trip is worse
than an honest unchecked sentence.
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path

from conftest import REPO_ROOT

from cell_curator.cli import build_parser

MATRIX = REPO_ROOT / "skills" / "cell-curator" / "references" / "quadbio-parity.md"

# Commands named in the matrix that are not subcommands of `cell-curator`.
# The first pair are the top-level flags, and the rest are Snakemake targets the
# matrix names alongside CLI commands in the same cell.
NON_ACTION_TOKENS = frozenset({"--version", "--help", "all_evidence", "review_packet",
                               "finalization", "writeback"})

BACKTICKED = re.compile(r"`([^`]+)`")
LINK = re.compile(r"\[`([^`]+)`\]\(([^)]+)\)")


def _leaf_actions(
    parser: argparse.ArgumentParser, prefix: tuple[str, ...] = ()
) -> set[str]:
    subparsers = [
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]
    if not subparsers:
        return {" ".join(prefix)}
    actions: set[str] = set()
    for name, child in subparsers[0].choices.items():
        actions |= _leaf_actions(child, (*prefix, name))
    return actions


def _matrix_rows() -> list[list[str]]:
    rows = []
    for line in MATRIX.read_text().splitlines():
        if not line.startswith("| ") or line.count("|") != 6:
            continue
        cells = [cell.strip() for cell in line.split("|")[1:-1]]
        if cells[0] == "QuadBio capability" or set(cells[0]) <= {"-", " "}:
            continue
        rows.append(cells)
    return rows


def _module_symbols() -> set[str]:
    """Every module-level name the package defines, plus its ``__all__`` entries."""

    names: set[str] = set()
    for path in (REPO_ROOT / "src" / "cell_curator").glob("*.py"):
        for node in ast.parse(path.read_text()).body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                names.add(node.name)
            elif isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "__all__"
                for target in node.targets
            ):
                names |= {
                    element.value
                    for element in node.value.elts
                    if isinstance(element, ast.Constant)
                }
    return names


def _member_names() -> set[str]:
    """Qualified names, so a method named in the matrix is checked where it lives."""

    names: set[str] = set()
    for path in (REPO_ROOT / "src" / "cell_curator").glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.ClassDef):
                continue
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    names.add(f"{node.name}.{child.name}")
    return names


def test_the_matrix_still_has_forty_five_mapped_capabilities() -> None:
    """A silently emptied or collapsed matrix would satisfy every other check here.

    Each check below iterates rows, so a table that lost its rows would pass by
    asserting nothing. The floor states the shape the review was performed
    against; raising it is a deliberate act, lowering it needs a reason.
    """

    assert len(_matrix_rows()) == 45


def test_every_implementation_path_in_the_matrix_is_real() -> None:
    rows = _matrix_rows()
    broken = [
        target
        for row in rows
        for _, target in LINK.findall(row[2])
        if not (MATRIX.parent / target).resolve().is_file()
    ]
    assert not broken, f"matrix names implementation paths that do not exist: {broken}"
    assert len(rows) >= 40, "a row without an implementation file would prove nothing"


def test_every_entry_point_in_the_matrix_resolves() -> None:
    """Backticked tokens in the entry-point column are commands or Python symbols."""

    actions = _leaf_actions(build_parser())
    symbols = _module_symbols()
    missing: list[str] = []
    for row in _matrix_rows():
        for token in BACKTICKED.findall(row[3]):
            token = token.strip()
            if token.endswith("()"):
                name = token[:-2]
                if name not in symbols and name not in _member_names():
                    missing.append(f"{name}() (row: {row[0][:40]})")
                continue
            if token and token[0].isupper():
                if token not in symbols:
                    missing.append(f"{token} (row: {row[0][:40]})")
                continue
            # A cell names commands as `group subcommand`, sometimes with a flag
            # that selects a mode. The action set holds the path, so the flag is
            # dropped and the path is compared whole.
            path = " ".join(word for word in token.split() if not word.startswith("-"))
            if path.startswith("cell-curator"):
                continue
            if path not in actions and token not in NON_ACTION_TOKENS:
                missing.append(f"`{token}` (row: {row[0][:40]})")
    assert not missing, f"matrix names entry points that do not resolve: {missing}"


def test_every_test_file_named_in_the_matrix_is_tracked() -> None:
    """A deleted test file is what the first packaging commit did to this document."""

    referenced = {
        Path(target).name
        for row in _matrix_rows()
        for _, target in LINK.findall(row[4])
    }
    existing = {path.name for path in (REPO_ROOT / "tests").glob("test_*.py")}
    assert not referenced - existing, f"matrix cites missing test files: {sorted(referenced - existing)}"
    assert len(referenced) >= 10, "the coverage column would be doing nothing if empty"


def test_the_matrix_states_that_no_comparison_is_published() -> None:
    """The disclaimer is the load-bearing sentence; its removal is a real change.

    Every row in this document is a capability claim, and the reason the
    document is safe to publish is that it explicitly refuses to be read as a
    result. Deleting that would turn a checklist into a superiority argument.
    """

    # The document is hard-wrapped, so the sentences are matched on collapsed
    # whitespace rather than on the source line breaks.
    text = " ".join(MATRIX.read_text().split())
    assert "No accuracy or performance claim is made or implied by any row above" in text
    assert "no such result has been published" in text
