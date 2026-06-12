#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

JOBNAME="hippocampus_LOCAL-RF_striatum-global-priors_ignore-seed"
CONFIG="paper-rf/configs/hippocampus_local_rf_striatum_global_priors_ignore_seed.toml"
CHAINS="${CHAINS:-4}"
CHAIN_START="${CHAIN_START:-1}"
SLURM_PARTITION="${SLURM_PARTITION:-long}"
SLURM_TIME="${SLURM_TIME:-5-00:00:00}"
export SLURM_PARTITION SLURM_TIME

POSTERIOR_PRIOR_SOURCE="$PROJECT_DIR/runs/striatum_DIFF-RF_RETRO_ignore-seed/posterior.h5"
if [[ ! -f "$POSTERIOR_PRIOR_SOURCE" ]]; then
  cat >&2 <<EOF
Missing required posterior-prior source:
  $POSTERIOR_PRIOR_SOURCE

Merge the existing striatum no-seed DIFF-RF chains first, then rerun this script.
This script does not submit merge jobs automatically.
EOF
  exit 1
fi

CHAIN_END=$((CHAIN_START + CHAINS - 1))
echo "Submitting ignore-seed posterior-prior hippocampus LOCAL-RF chains C${CHAIN_START}-C${CHAIN_END} on partition '$SLURM_PARTITION' for $SLURM_TIME."
for CHAIN in $(seq "$CHAIN_START" "$CHAIN_END"); do
  RUN_ID="${JOBNAME}_C${CHAIN}"
  echo "Submitting $RUN_ID from $CONFIG"
  "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG" "$RUN_ID"
done

echo "All ignore-seed posterior-prior hippocampus LOCAL-RF jobs submitted."
