#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

declare -A CONFIGS
CONFIGS["hippocampus_DIFF-RF_RETRO"]="configs/paper/hippocampus_diff_rf_core.toml"
CONFIGS["hippocampus_DIFF-RF_RETRO_striatum-global-priors"]="configs/paper/hippocampus_diff_rf_striatum_global_priors.toml"

POSTERIOR_PRIOR_SOURCE="$PROJECT_DIR/runs/striatum_DIFF-RF_RETRO/posterior.h5"
if [[ ! -f "$POSTERIOR_PRIOR_SOURCE" ]]; then
  cat >&2 <<EOF
Missing posterior-prior source:
  $POSTERIOR_PRIOR_SOURCE

Create or sync the merged striatum DIFF-RF run before submitting the posterior-prior hippocampus job.
For example, merge striatum_DIFF-RF_RETRO_C1..C4 into runs/striatum_DIFF-RF_RETRO.
EOF
  exit 1
fi

for JOBNAME in "${!CONFIGS[@]}"; do
  CONFIG_PATH="${CONFIGS[$JOBNAME]}"
  for CHAIN in {1..4}; do
    RUN_ID="${JOBNAME}_C${CHAIN}"
    echo "Submitting $RUN_ID from $CONFIG_PATH"
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG_PATH" "$RUN_ID"
  done
done

echo "All hippocampus DIFF-RF retrograde jobs submitted."
