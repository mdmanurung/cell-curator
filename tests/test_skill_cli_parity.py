"""SKILL.md and the CLI must describe the same tool.

`test_cli_api_parity.py` holds the CLI to the Python API, but nothing held the
CLI to the agent skill, and the drift ran both ways: nine optional routes were
documented as bare `--help` probes whose required flags appeared nowhere in the
skill, while sixteen actions existed only in the CLI. An agent reading SKILL.md
can only run what SKILL.md spells out, so an undocumented required flag is a
defect in the same sense a missing function would be.

The sibling test in `test_cli_skill_guided_gates.py` checks that each documented
command *path* resolves. This one checks the flags on it, and the reverse
direction.
"""

from __future__ import annotations

import argparse
import io
import re
import shlex
from contextlib import redirect_stderr, redirect_stdout

import pytest
from conftest import SKILL_ROOT

from cell_curator.cli import build_parser

SKILL_PATH = SKILL_ROOT / "SKILL.md"

# Actions the skill deliberately omits. `pipeline` exposes the individual phases
# that `run execute` already sequences; the rest serve packaging, introspection,
# and benchmarking rather than an annotation run. Documenting an action removes
# it from here — this set only ever shrinks.
INTERNAL_ACTIONS: frozenset[tuple[str, ...]] = frozenset(
    {
        ("benchmark", "evaluate"),
        ("benchmark", "validate-manifest"),
        ("capabilities",),
        ("config", "schema"),
        ("knowledge", "validate"),
        ("ontology", "validate"),
        ("pipeline", "audit"),
        ("pipeline", "candidate-evidence-lane"),
        ("pipeline", "decide"),
        ("pipeline", "evidence"),
        ("pipeline", "evidence-lane"),
        ("pipeline", "map"),
        ("pipeline", "packet"),
        ("pipeline", "refine"),
        ("pipeline", "review"),
        ("run", "init"),
    }
)


def _leaf_actions(
    parser: argparse.ArgumentParser, prefix: tuple[str, ...] = ()
) -> list[tuple[str, ...]]:
    subparsers = [
        action for action in parser._actions if isinstance(action, argparse._SubParsersAction)
    ]
    if not subparsers:
        return [prefix]
    leaves: list[tuple[str, ...]] = []
    for action in subparsers:
        for name, sub in action.choices.items():
            leaves.extend(_leaf_actions(sub, (*prefix, name)))
    return leaves


def _documented_invocations() -> list[str]:
    """Every `cell-curator ...` line in the skill, backslash continuations joined."""

    text = SKILL_PATH.read_text().replace("\\\n", " ")
    return sorted({line.strip() for line in re.findall(r"^cell-curator .*$", text, re.MULTILINE)})


def _documented_actions() -> set[tuple[str, ...]]:
    actions = set()
    for invocation in _documented_invocations():
        words = shlex.split(invocation)[1:]
        actions.add(tuple(word for word in words[:2] if not word.startswith("-")))
    return actions


@pytest.mark.parametrize("invocation", _documented_invocations())
def test_every_documented_invocation_parses(invocation: str) -> None:
    """A flag the skill writes must be a flag the parser accepts.

    `--help` short-circuits argparse before it validates anything, so a skill
    that documents a route only as `<route> --help` would pass while hiding
    every required flag. Reject those outright.
    """

    argv = shlex.split(invocation)[1:]
    assert "--help" not in argv, (
        f"{invocation!r} documents a route as a help probe. "
        "Write the invocation with its required flags instead."
    )
    buffer = io.StringIO()
    try:
        with redirect_stdout(buffer), redirect_stderr(buffer):
            build_parser().parse_args(argv)
    except SystemExit as exit_:  # pragma: no cover - only on failure
        pytest.fail(f"{invocation!r} does not parse (exit {exit_.code}): {buffer.getvalue()}")


def test_every_cli_action_is_documented_or_explicitly_internal() -> None:
    leaves = set(_leaf_actions(build_parser()))
    documented = _documented_actions()

    assert not documented - leaves, (
        f"SKILL.md documents actions the CLI does not have: {sorted(documented - leaves)}"
    )
    undocumented = leaves - documented - INTERNAL_ACTIONS
    assert not undocumented, (
        f"CLI actions missing from SKILL.md: {sorted(undocumented)}. "
        "Document them, or add them to INTERNAL_ACTIONS with a reason."
    )
    stale = INTERNAL_ACTIONS & documented
    assert not stale, f"INTERNAL_ACTIONS lists actions the skill documents: {sorted(stale)}"


def test_the_documented_action_count_cannot_silently_collapse() -> None:
    """A floor, so deleting the skill's command blocks cannot make parity trivial."""

    assert len(_documented_actions()) >= 44
