#!/usr/bin/env bash
# cd "recur inter" && bash pipeline/run.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export COMBO_PROJECT_ROOT="$ROOT"
export COMBO_MATA_FROZEN_DIR="$ROOT/pipeline/frozen"
REC_K="${1:-5,10,15}"

TRAIN_RDS="$ROOT/pipeline/data/training/e_merged_raw.rds"
FROZEN_LIST="$ROOT/pipeline/frozen/frozen_A_sig_top10_genes.tsv"

# Step 0b downloads recurrence + external GEO and writes integrated RDS under pipeline/data/.
# Step 0a (metastasis GEO) runs only when frozen lists are missing.
if [[ ! -f "$TRAIN_RDS" ]]; then
  if [[ -f "$FROZEN_LIST" ]]; then
    echo "Training RDS missing; running step 0b (download + integrate) using pipeline/frozen/ ..."
    Rscript --vanilla "$ROOT/pipeline/steps/step0b_integrate.R"
  else
    echo "Frozen metastasis lists missing; running step 0a then 0b ..."
    bash "$ROOT/pipeline/steps/step0a_metastasis.sh"
    Rscript --vanilla "$ROOT/pipeline/steps/step0b_integrate.R"
  fi
fi

mkdir -p "$ROOT/pipeline/model result"
exec Rscript --vanilla "$ROOT/pipeline/run.R" "$REC_K" 2>&1 | tee "$ROOT/pipeline/model result/pipeline_run.log"
