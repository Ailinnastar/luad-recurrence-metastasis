#!/usr/bin/env bash
# Step 0a — download metastasis GEO, build the gene support table, freeze top-K lists.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
STEPS="$PIPELINE_DIR/steps"

export COMBO_METASTASIS_INPUTS="${COMBO_METASTASIS_INPUTS:-$PIPELINE_DIR/geo/metastasis}"
export COMBO_METASTASIS_RESULTS="${COMBO_METASTASIS_RESULTS:-$PIPELINE_DIR/metastasis_results}"
export COMBO_MATA_FROZEN_DIR="${COMBO_MATA_FROZEN_DIR:-$PIPELINE_DIR/frozen}"
mkdir -p "$COMBO_METASTASIS_INPUTS" "$COMBO_MATA_FROZEN_DIR" "$COMBO_METASTASIS_RESULTS"

PYTHON="${PYTHON:-python3}"
Rscript --vanilla "$STEPS/step0a_download.R"
"$PYTHON" -m pip install -q -r "$STEPS/requirements-modeling.txt" 2>/dev/null || true
"$PYTHON" "$STEPS/step0a_metastasis.py"
