#!/usr/bin/env bash
# Purpose: take the PBMC3K consensus run into cell-curator, up to the write-back preview.
# Run:     bash 05_package-run.sh RUN PKG_CONFIG   (CONTINUE=1 after the two person steps)
# Inputs:  RUN/out_h5/{consensus.tsv,run_manifest.json,work/cells.tsv}  [01 to 04 above]
#          PKG_CONFIG: a cell-curator config with annotation.enabled: true,
#          guidance.hierarchy_levels: [L1] and annotation.consensus_directory: RUN/out_h5
# Outputs: the run directory named by PKG_CONFIG (labels, review report, write-back preview)
# Four steps stay with a person and the script stops at each: Gate A (assumptions approve),
# corpus review, critic findings and write-back apply. It never runs them.
set -euo pipefail
RUN=$1
CFG=$2
CC="cell-curator"

# The package plans first, because the consensus run needs plan/vocabulary.tsv. -----------
$CC run scene --config "$CFG"
$CC run plan --config "$CFG"
echo "STOP 1: a person runs: $CC assumptions approve --config $CFG --reviewer NAME"
echo "STOP 2: a person reviews the corpus proposal: $CC corpus review ... --reviewer NAME"
echo "Then rerun this script with CONTINUE=1."
[ "${CONTINUE:-0}" = 1 ] || exit 0

# import-consensus refuses labels outside the plan's vocabulary, so no separate check here.

# Evidence, import, finalize, review. -----------------------------------------------------
$CC annotation collect-level --config "$CFG" --level L1
$CC annotation import-consensus --config "$CFG" --level L1
$CC annotation finalize --config "$CFG"
$CC run status --config "$CFG"
echo "STOP 3: critics record findings or attestations, then: $CC run validate --config $CFG"
echo "STOP 4: $CC write-back preview --config $CFG --output OUT.h5ad; a person runs write-back apply"
