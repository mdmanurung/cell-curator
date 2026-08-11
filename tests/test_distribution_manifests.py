"""The plugin and marketplace manifests are a distribution surface with no runtime.

Nothing imports `.claude-plugin/`, so a stale version or a path that stopped
resolving cannot fail any other test — it fails at `/plugin install`, on someone
else's machine, after the tag is already public. These assertions are the only
thing standing between a version bump and that.
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

from conftest import REPO_ROOT, SKILL_ROOT

PLUGIN_ROOT = REPO_ROOT / ".claude-plugin"


def _project_version() -> str:
    return tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())["project"]["version"]


def _manifest(name: str) -> dict:
    return json.loads((PLUGIN_ROOT / name).read_text())


def test_every_declared_version_matches_the_project_version() -> None:
    """One release, one version string, three files that each repeat it."""

    expected = _project_version()
    assert _manifest("plugin.json")["version"] == expected
    assert _manifest("marketplace.json")["metadata"]["version"] == expected


def test_plugin_manifest_skill_path_resolves_to_the_shipped_skill() -> None:
    """`skills` is a path Claude Code follows; a rename here is silent otherwise."""

    plugin = _manifest("plugin.json")
    skills_root = (PLUGIN_ROOT.parent / plugin["skills"]).resolve()
    assert skills_root.is_dir()
    assert SKILL_ROOT.resolve().parent == skills_root

    # The four members the skill contract promises. Packaging repeats this list in
    # pyproject's data-files, so a directory dropped from one and not the other is
    # the exact drift worth catching.
    assert (SKILL_ROOT / "SKILL.md").is_file()
    for directory in ("agents", "assets", "references"):
        assert (SKILL_ROOT / directory).is_dir()
        assert list((SKILL_ROOT / directory).iterdir())


def test_marketplace_entry_points_at_a_directory_holding_the_plugin_manifest() -> None:
    """`marketplace add` reads this; `plugin install` then reads what it points to."""

    marketplace = _manifest("marketplace.json")
    entries = marketplace["plugins"]
    assert [entry["name"] for entry in entries] == [_manifest("plugin.json")["name"]]
    for entry in entries:
        source = (PLUGIN_ROOT.parent / entry["source"]).resolve()
        assert (source / ".claude-plugin" / "plugin.json").is_file()


def test_readme_quickstart_paths_exist() -> None:
    """The README tells users to copy these before anything else works."""

    assert (SKILL_ROOT / "assets" / "config.template.yaml").is_file()
    assert (REPO_ROOT / "workflow" / "Snakefile").is_file()
    assert (REPO_ROOT / "workflow" / "profiles" / "slurm" / "config.yaml").is_file()


def test_no_tracked_path_the_manifests_promise_is_git_ignored() -> None:
    """A shipped file that `.gitignore` swallows reaches nobody who clones."""

    ignored = {line.strip() for line in (REPO_ROOT / ".gitignore").read_text().splitlines()}
    promised = [Path(".claude-plugin"), Path("skills"), Path("workflow"), Path("schemas")]
    for path in promised:
        assert f"{path}/" not in ignored and path.name not in ignored
