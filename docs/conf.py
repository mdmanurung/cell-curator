"""Sphinx configuration.

The narrative pages are not written here — they are generated from the
repository's Markdown by :mod:`_sync`, which runs before every build. See that
module for why the site is generated rather than authored.
"""

from __future__ import annotations

import sys
from importlib.metadata import version as _version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _sync  # noqa: E402

_sync.sync()

project = "cell-curator"
author = "cell-curator contributors"
copyright = "cell-curator contributors"
release = _version("cell-curator")
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "deflist", "fieldlist"]
myst_heading_anchors = 3

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": False}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
}

html_theme = "furo"
html_title = f"cell-curator {release}"
html_theme_options = {
    "source_repository": "https://github.com/mdmanurung/cell-curator",
    "source_branch": "main",
    "source_view_link": f"{_sync.BLOB_BASE}/docs/{{filename}}",
}


def _view_source_at_its_real_origin(app, pagename, templatename, context, doctree):
    """Point "View this page" at the file a reader can actually edit.

    Furo derives the link from the page's own path, which for a generated page
    is ``docs/_generated/...`` — ignored by git, so the link 404s. Every page
    here has a real origin; :data:`_sync.PAGES` knows where.
    """

    del app, templatename, doctree
    origins = {f"_generated/{page[:-3]}": source for source, page in _sync.PAGES.items()}
    origins["_generated/api"] = "docs/_sync.py"

    if pagename in origins:
        context["theme_source_view_link"] = f"{_sync.BLOB_BASE}/{origins[pagename]}"
    elif pagename.startswith("_generated/_api/"):
        # An autosummary stub has no source file of its own to view.
        context["theme_source_view_link"] = ""


def setup(app):
    app.connect("html-page-context", _view_source_at_its_real_origin)

# The generated tree carries GitHub URLs for files that are not pages here.
# Checking them on every build would make the docs depend on network reachability.
linkcheck_ignore = [r"https://github\.com/mdmanurung/cell-curator/blob/.*"]
