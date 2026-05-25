#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

POSTERIOR_JOBNAME="hippocampus_DIFF-RF_RETRO_striatum-global-priors"
POSTERIOR_CONFIG="configs/paper/hippocampus_diff_rf_striatum_global_priors.toml"
POSTERIOR_CHAINS="${POSTERIOR_CHAINS:-4}"
SLURM_PARTITION="${SLURM_PARTITION:-long}"
SLURM_TIME="${SLURM_TIME:-5-00:00:00}"
export SLURM_PARTITION SLURM_TIME

POSTERIOR_PRIOR_SOURCE="$PROJECT_DIR/runs/striatum_DIFF-RF_RETRO/posterior.h5"
MERGE_DEPENDENCY=""
if [[ ! -f "$POSTERIOR_PRIOR_SOURCE" ]]; then
  echo "Merged striatum DIFF-RF posterior not found; submitting merge job first." >&2
  MERGE_JOB_ID="$("$PROJECT_DIR/cluster/submit_merge_chains.sh" striatum_DIFF-RF_RETRO striatum_DIFF-RF_RETRO 4)"
  MERGE_DEPENDENCY="afterok:$MERGE_JOB_ID"
  echo "Posterior-prior hippocampus chains will wait for merge job $MERGE_JOB_ID." >&2
fi

echo "Submitting $POSTERIOR_CHAINS posterior-prior hippocampus chains on partition '$SLURM_PARTITION' for $SLURM_TIME."
for CHAIN in $(seq 1 "$POSTERIOR_CHAINS"); do
  RUN_ID="${POSTERIOR_JOBNAME}_C${CHAIN}"
  echo "Submitting $RUN_ID from $POSTERIOR_CONFIG"
  if [[ -n "$MERGE_DEPENDENCY" ]]; then
    SLURM_DEPENDENCY="$MERGE_DEPENDENCY" "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$POSTERIOR_CONFIG" "$RUN_ID"
  else
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$POSTERIOR_CONFIG" "$RUN_ID"
  fi
done

echo "All posterior-prior hippocampus DIFF-RF retrograde jobs submitted."
