#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

declare -A CONFIGS
CONFIGS["striatum_DIFF_RETRO"]="configs/paper/striatum_core.toml"
CONFIGS["striatum_DIFF-R_RETRO"]="configs/paper/striatum_diff_r_core.toml"
CONFIGS["striatum_DIFF-RF_RETRO"]="configs/paper/striatum_diff_rf_core.toml"
CONFIGS["hippocampus_DIFF_RETRO"]="configs/paper/hippocampus_core.toml"
CONFIGS["hippocampus_DIFF-R_RETRO"]="configs/paper/hippocampus_diff_r_core.toml"
CONFIGS["hippocampus_DIFF-RF_RETRO"]="configs/paper/hippocampus_diff_rf_core.toml"

for JOBNAME in "${!CONFIGS[@]}"; do
  CONFIG_PATH="${CONFIGS[$JOBNAME]}"
  for CHAIN in {1..4}; do
    RUN_ID="${JOBNAME}_C${CHAIN}"
    echo "Submitting $RUN_ID from $CONFIG_PATH"
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG_PATH" "$RUN_ID"
  done
done

echo "All paper retrograde jobs submitted."
