#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

JOBNAME="hippocampus_DIFF-RF_RETRO_striatum-global-priors_ignore-seed"
CONFIG="paper-rf/configs/hippocampus_diff_rf_striatum_global_priors_ignore_seed.toml"
CHAINS="${CHAINS:-4}"
CHAIN_START="${CHAIN_START:-1}"
SLURM_PARTITION="${SLURM_PARTITION:-long}"
SLURM_TIME="${SLURM_TIME:-5-00:00:00}"
export SLURM_PARTITION SLURM_TIME

POSTERIOR_PRIOR_SOURCE="$PROJECT_DIR/runs/striatum_DIFF-RF_RETRO/posterior.h5"
MERGE_DEPENDENCY=""
if [[ ! -f "$POSTERIOR_PRIOR_SOURCE" ]]; then
  echo "Merged striatum DIFF-RF posterior not found; submitting merge job first." >&2
  MERGE_JOB_ID="$("$PROJECT_DIR/cluster/submit_merge_chains.sh" striatum_DIFF-RF_RETRO striatum_DIFF-RF_RETRO 4)"
  MERGE_DEPENDENCY="afterok:$MERGE_JOB_ID"
  echo "Ignore-seed hippocampus chains will wait for merge job $MERGE_JOB_ID." >&2
fi

CHAIN_END=$((CHAIN_START + CHAINS - 1))
echo "Submitting ignore-seed posterior-prior hippocampus chains C${CHAIN_START}-C${CHAIN_END} on partition '$SLURM_PARTITION' for $SLURM_TIME."
for CHAIN in $(seq "$CHAIN_START" "$CHAIN_END"); do
  RUN_ID="${JOBNAME}_C${CHAIN}"
  echo "Submitting $RUN_ID from $CONFIG"
  if [[ -n "$MERGE_DEPENDENCY" ]]; then
    SLURM_DEPENDENCY="$MERGE_DEPENDENCY" "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG" "$RUN_ID"
  else
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG" "$RUN_ID"
  fi
done

echo "All ignore-seed posterior-prior hippocampus DIFF-RF retrograde jobs submitted."
