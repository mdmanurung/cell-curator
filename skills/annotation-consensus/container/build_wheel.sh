#!/usr/bin/env bash
# Purpose: build the cell_curator wheel from a clean `git archive` of a commit, and write cell-curator.txt.
# Run:     bash build_wheel.sh [COMMIT]     (default HEAD; the repo is not touched)
# Outputs: wheelhouse/cell_curator-*.whl, cell-curator.txt (name, version, wheel sha256), _build/wheel_source.txt
# Notes:   SOURCE_DATE_EPOCH is the time of the last commit that changed what the wheel ships, so the
#          wheel hash stays the same across later commits that change nothing it ships.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
COMMIT=${1:-HEAD}
COMMIT=$(git -C "$REPO" rev-parse "$COMMIT")
export SOURCE_DATE_EPOCH=$(git -C "$REPO" log -1 --format=%ct "$COMMIT" -- src pyproject.toml README.md LICENSE schemas workflow)
mkdir -p "$HERE/_build" "$HERE/wheelhouse"
SRC=$(mktemp -d "$HERE/_build/wheelsrc.XXXXXX")
git -C "$REPO" archive "$COMMIT" | tar -x -C "$SRC"
rm -f "$HERE"/wheelhouse/cell_curator-*.whl
(cd "$SRC" && uv build --wheel --python 3.11 --out-dir "$HERE/wheelhouse" .)
rm -rf "$SRC" "$HERE/wheelhouse/.gitignore"  # uv build writes a "*" ignore file here
WHL=$(ls "$HERE"/wheelhouse/cell_curator-*.whl)
NAME=$(basename "$WHL" | cut -d- -f2)
SHA=$(sha256sum "$WHL" | cut -d' ' -f1)
printf 'cell-curator==%s \\\n    --hash=sha256:%s\n' "$NAME" "$SHA" > "$HERE/cell-curator.txt"
{ echo "commit  $COMMIT"; echo "epoch   $SOURCE_DATE_EPOCH"; echo "wheel   $(basename "$WHL")"; echo "sha256  $SHA"; } | tee "$HERE/_build/wheel_source.txt"
