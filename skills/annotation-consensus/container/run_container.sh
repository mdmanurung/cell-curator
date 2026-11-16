#!/usr/bin/env bash
# Purpose: run the skill inside the SIF with a clean environment: no host home, the repo read-only,
#          only OUTDIR's parent writable. The SIF sha256 goes to the manifest (CONSENSUS_SIF_SHA256).
# Run:     bash run_container.sh SIF OUTDIR [skill args...]    (default args: --demo; OUTDIR must not exist)
set -euo pipefail
APPTAINER=${APPTAINER:-$(command -v apptainer || echo /share/software/apptainer/1.2.5/bin/apptainer)}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
SIF=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); shift 2
ARGS=("$@"); [ ${#ARGS[@]} -gt 0 ] || ARGS=(--demo)
export APPTAINER_CACHEDIR=$HERE/_build/cache APPTAINER_TMPDIR=$HERE/_build/tmp
mkdir -p "$(dirname "$OUT")" "$APPTAINER_TMPDIR"
"$APPTAINER" exec --cleanenv --containall --env CONSENSUS_SIF_SHA256="$(sha256sum "$SIF" | cut -d' ' -f1)" \
  --bind "$REPO:/repo:ro" --bind "$(dirname "$OUT"):/out" "$SIF" \
  python /repo/skills/annotation-consensus/annotation_consensus.py "${ARGS[@]}" --output "/out/$(basename "$OUT")"
