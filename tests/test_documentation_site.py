"""The documentation site is generated, so the generator is what needs testing.

A copy step that runs without raising proves nothing: the parity reference was
copied faithfully for the life of the repository and every link in it was dead.
What matters is that no relative link survives into the published tree, and that
the public API cannot grow an export the site does not document.

These tests import the generator rather than running Sphinx, so they cost
milliseconds and run in the same suite as everything else.
"""

from __future__ import annotations

import re
import sys

import pytest
from conftest import REPO_ROOT

sys.path.insert(0, str(REPO_ROOT / "docs"))

import _sync  # noqa: E402

LINK = re.compile(r"(?<=\]\()(?!https?:|mailto:|#)([^)]+)(?=\))")


@pytest.mark.parametrize("source", sorted(_sync.PAGES))
def test_no_relative_link_survives_into_the_published_page(source: str) -> None:
    """Anything left relative is either a broken URL or a 404 on the site."""

    rendered = _sync.render(source)
    page_targets = {value for value in _sync.PAGES.values()}
    surviving = []
    for target in LINK.findall(rendered):
        path = target.split("#", 1)[0]
        # A link to another generated page is relative *by design* and must
        # resolve to a real page, not merely look plausible.
        candidate = (REPO_ROOT / "docs" / "_generated" / _sync.PAGES[source]).parent / path
        try:
            resolved = candidate.resolve().relative_to(REPO_ROOT / "docs" / "_generated")
        except ValueError:
            surviving.append(target)
            continue
        if resolved.as_posix() not in page_targets:
            surviving.append(target)
    assert not surviving, f"{source} keeps relative links that go nowhere: {surviving}"


def test_generator_refuses_a_link_it_cannot_classify() -> None:
    """The raise is the whole design; a pass-through would restore the old failure."""

    with pytest.raises(_sync.UnclassifiableLink):
        _sync.rewrite_link("../does/not/exist.md", source="README.md")


def test_every_public_export_appears_in_the_api_page() -> None:
    """Driving the page from __all__ is only useful if nothing filters names out."""

    import cell_curator

    page = _sync.render_api()
    missing = [name for name in cell_curator.__all__ if not re.search(rf"\b{name}\b", page)]
    assert not missing, f"exported but undocumented: {missing}"


def test_every_reference_document_is_published() -> None:
    """A new contract file must not be able to exist without reaching the site."""

    on_disk = {
        path.name for path in (REPO_ROOT / "skills/cell-curator/references").glob("*.md")
    }
    published = {
        source.rsplit("/", 1)[1]
        for source in _sync.PAGES
        if "references/" in source
    }
    assert on_disk == published
