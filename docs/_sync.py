"""Copy the repository's Markdown into the docs tree, rewriting relative links.

The Markdown that documents this package has to serve three readers at once: an
agent loading `SKILL.md` from an installed skill directory, someone browsing the
repository on GitHub, and this documentation site. Only the first two can use a
path like ``../../../src/cell_curator/data.py``.

So the site is generated, never authored: the sources stay where the other two
readers need them, and this module rewrites their links on the way in. A link to
another generated page becomes a local reference; a link to anything else in the
repository becomes a GitHub URL at ``BLOB_BASE``; a link to neither raises. That
last branch is the point — a silently passed-through link is how the parity
reference shipped 134 dead ones.
"""

from __future__ import annotations

import re
from pathlib import Path

DOCS = Path(__file__).resolve().parent
REPO_ROOT = DOCS.parent
GENERATED = DOCS / "_generated"
BLOB_BASE = "https://github.com/mdmanurung/cell-curator/blob/main"

LINK = re.compile(r"(?<=\]\()(?!https?:|mailto:|#)([^)]+)(?=\))")

#: Repository Markdown, and the page each becomes. Order is display order.
PAGES: dict[str, str] = {
    "README.md": "readme.md",
    "skills/cell-curator/SKILL.md": "skill.md",
    **{
        f"skills/cell-curator/references/{name}.md": f"reference/{name}.md"
        for name in (
            "input-contract",
            "assay-routing",
            "adaptive-subclustering",
            "multimodal-evidence",
            "decision-gates",
            "artifact-schemas",
            "compute-and-snakemake",
            "validation",
            "quadbio-parity",
        )
    },
}


class UnclassifiableLink(RuntimeError):
    """A relative link points at neither a generated page nor a repository file."""


def rewrite_link(target: str, *, source: str) -> str:
    """Rewrite one relative link found in ``source`` (a repo-relative Markdown path)."""

    path, _, fragment = target.partition("#")
    suffix = f"#{fragment}" if fragment else ""
    if not path:
        return target

    resolved = ((REPO_ROOT / source).parent / path).resolve()
    try:
        repo_relative = resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError as error:  # pragma: no cover - a link escaping the repo
        raise UnclassifiableLink(f"{source}: {target} escapes the repository") from error

    if repo_relative in PAGES:
        here = Path(PAGES[source]).parent
        there = Path(PAGES[repo_relative])
        import os.path

        return os.path.relpath(there, here) + suffix

    if resolved.exists():
        # GitHub serves directories under /tree and redirects /blob to it.
        kind = "tree" if resolved.is_dir() else "blob"
        return f"{BLOB_BASE.replace('/blob', f'/{kind}')}/{repo_relative}{suffix}"

    raise UnclassifiableLink(
        f"{source} links to {target}, which is neither a documentation page "
        f"nor a file in the repository"
    )


def render(source: str) -> str:
    """Return the text of one repository Markdown file with its links rewritten."""

    text = (REPO_ROOT / source).read_text()
    return LINK.sub(lambda match: rewrite_link(match.group(1), source=source), text)


def render_api() -> str:
    """Build the API page from ``cell_curator.__all__``.

    Curating this list by hand would mean a new export is undocumented until
    someone remembers the page exists. Driving it from ``__all__`` means the
    public surface and its documentation cannot disagree.
    """

    import cell_curator

    exported = sorted(cell_curator.__all__)
    # CellCurator has its own section below; autosummary would document it a
    # second time and every method with it.
    classes = [name for name in exported if name[:1].isupper() and name != "CellCurator"]
    functions = [name for name in exported if not name[:1].isupper()]

    def table(names: list[str]) -> str:
        # `currentmodule` rather than fully-qualified entries: autosummary cannot
        # import a bare name without it, and the qualified form would print
        # "cell_curator." in front of all 73 rows.
        entries = "\n".join(f"   {name}" for name in names)
        return (
            "```{eval-rst}\n"
            ".. currentmodule:: cell_curator\n\n"
            ".. autosummary::\n"
            "   :toctree: _api\n"
            "   :nosignatures:\n\n"
            f"{entries}\n"
            "```"
        )

    return f"""# Python API

Everything `cell_curator` exports. The guided object is the usual entry point; the
functions beneath it are the same implementation the CLI and the Snakemake DAG call.

## Guided object

```{{eval-rst}}
.. autoclass:: cell_curator.CellCurator
   :members:
   :member-order: bysource
```

## Configuration and record types

{table(classes)}

## Functions

{table(functions)}
"""


def sync() -> list[Path]:
    """Write every page in :data:`PAGES`, plus the API page, into ``docs/_generated``."""

    written = []
    for source, page in PAGES.items():
        destination = GENERATED / page
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(render(source))
        written.append(destination)

    api = GENERATED / "api.md"
    api.write_text(render_api())
    written.append(api)
    return written


if __name__ == "__main__":
    for path in sync():
        print(path.relative_to(REPO_ROOT))
