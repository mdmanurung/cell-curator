#!/usr/bin/env bash
# Purpose: build _build/NAME.sif from consensus.def with --fakeroot, and record sha256, size, base digest.
# Run:     sbatch build_sif.sbatch [NAME]      (or: bash build_sif.sh [NAME]; run build_wheel.sh first)
# Outputs: _build/NAME.sif, _build/NAME.record.txt, _build/build_NAME.log
set -euo pipefail
APPTAINER=${APPTAINER:-$(command -v apptainer || echo /share/software/apptainer/1.2.5/bin/apptainer)}
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
NAME=${1:-consensus}
export APPTAINER_CACHEDIR=$HERE/_build/cache APPTAINER_TMPDIR=$HERE/_build/tmp
mkdir -p _build/cache _build/tmp
rm -f "_build/$NAME.sif"
"$APPTAINER" build --fakeroot "_build/$NAME.sif" consensus.def > "_build/build_$NAME.log" 2>&1
{ sha256sum "_build/$NAME.sif"; stat -c 'size_bytes %s' "_build/$NAME.sif"; cat _build/wheel_source.txt; date -u +'built_utc %FT%TZ'; } > "_build/$NAME.record.txt"
cat "_build/$NAME.record.txt"
