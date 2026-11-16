#!/usr/bin/env bash
# Purpose: remake lock/requirements.lock (hashes) from lock/requirements.in. Not needed to build or run.
# Run:     bash make_lock.sh      (needs network for PyPI metadata)
# Notes:   versions are held to the host .venv (lock/host_freeze.txt, from `uv pip freeze`).
#          connection-pool 0.0.3 (a snakemake dependency) has no wheel on PyPI. Its sdist hash was
#          checked (bf429e7a...), a wheel was built from it with SOURCE_DATE_EPOCH fixed, stored in
#          wheelhouse/, and the lock line carries the wheel hash instead.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE/lock"
uv pip compile requirements.in --constraints host_freeze.txt --generate-hashes \
  --python-version 3.11 --python-platform x86_64-manylinux_2_28 --no-header --no-annotate \
  -o ../_build/requirements.compiled
CP=$(sha256sum ../wheelhouse/connection_pool-0.0.3-py3-none-any.whl | cut -d' ' -f1)
python3 - "$CP" <<'PY'
import re, sys
txt = open("../_build/requirements.compiled").read()
new, n = re.subn(r"connection-pool==0\.0\.3 \\\n    --hash=sha256:[0-9a-f]+",
                 f"connection-pool==0.0.3 \\\n    --hash=sha256:{sys.argv[1]}", txt)
if n != 1: raise SystemExit("connection-pool line not found")
open("requirements.lock", "w").write(new)
PY
