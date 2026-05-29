#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

JOBNAME="striatum_DIFF-RF_RETRO_ignore-seed"
CONFIG="paper-rf/configs/striatum_diff_rf_ignore_seed.toml"
CHAINS="${CHAINS:-4}"
CHAIN_START="${CHAIN_START:-1}"
SLURM_PARTITION="${SLURM_PARTITION:-long}"
SLURM_TIME="${SLURM_TIME:-5-00:00:00}"
export SLURM_PARTITION SLURM_TIME

CHAIN_END=$((CHAIN_START + CHAINS - 1))
echo "Submitting ignore-seed striatum DIFF-RF chains C${CHAIN_START}-C${CHAIN_END} on partition '$SLURM_PARTITION' for $SLURM_TIME."
for CHAIN in $(seq "$CHAIN_START" "$CHAIN_END"); do
  RUN_ID="${JOBNAME}_C${CHAIN}"
  echo "Submitting $RUN_ID from $CONFIG"
  "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG" "$RUN_ID"
done

echo "All ignore-seed striatum DIFF-RF retrograde jobs submitted."
