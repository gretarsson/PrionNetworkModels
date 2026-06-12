#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

CHAINS="${CHAINS:-4}"
CHAIN_START="${CHAIN_START:-1}"
SLURM_PARTITION="${SLURM_PARTITION:-long}"
SLURM_TIME="${SLURM_TIME:-5-00:00:00}"
export SLURM_PARTITION SLURM_TIME

SEED_PRIOR_SOURCE="$PROJECT_DIR/runs/striatum_DIFF-RF_RETRO/posterior.h5"
NOSEED_PRIOR_SOURCE="$PROJECT_DIR/runs/striatum_DIFF-RF_RETRO_ignore-seed/posterior.h5"

missing=0
if [[ ! -f "$SEED_PRIOR_SOURCE" ]]; then
  echo "Missing required seed-included striatum prior source:" >&2
  echo "  $SEED_PRIOR_SOURCE" >&2
  missing=1
fi
if [[ ! -f "$NOSEED_PRIOR_SOURCE" ]]; then
  echo "Missing required no-seed striatum prior source:" >&2
  echo "  $NOSEED_PRIOR_SOURCE" >&2
  missing=1
fi
if [[ "$missing" == "1" ]]; then
  cat >&2 <<EOF

Merge the existing striatum chains first, then rerun this script.
This script does not submit merge jobs automatically.
EOF
  exit 1
fi

submit_grid_item() {
  local jobname="$1"
  local config="$2"
  local chain_end=$((CHAIN_START + CHAINS - 1))

  echo "Submitting $jobname chains C${CHAIN_START}-C${chain_end} from $config on partition '$SLURM_PARTITION' for $SLURM_TIME."
  for chain in $(seq "$CHAIN_START" "$chain_end"); do
    local run_id="${jobname}_C${chain}"
    echo "Submitting $run_id"
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$config" "$run_id"
  done
}

# 1. DIFF-RF no-seed, no global prior.
submit_grid_item \
  "hippocampus_DIFF-RF_RETRO_ignore-seed_no-global-prior" \
  "paper-rf/configs/hippocampus_diff_rf_ignore_seed.toml"

# 2. DIFF-RF no-seed, global prior from striatum no-seed DIFF-RF.
submit_grid_item \
  "hippocampus_DIFF-RF_RETRO_striatum-noseed-global-priors_ignore-seed" \
  "paper-rf/configs/hippocampus_diff_rf_striatum_noseed_global_priors_ignore_seed.toml"

# 3. LOCAL-RF no-seed, no global prior.
submit_grid_item \
  "hippocampus_LOCAL-RF_ignore-seed_no-global-prior" \
  "paper-rf/configs/hippocampus_local_rf_ignore_seed.toml"

# 4. LOCAL-RF no-seed, global prior from striatum no-seed DIFF-RF.
submit_grid_item \
  "hippocampus_LOCAL-RF_striatum-noseed-global-priors_ignore-seed" \
  "paper-rf/configs/hippocampus_local_rf_striatum_global_priors_ignore_seed.toml"

# 5. LOCAL-RF seed-included, no global prior.
submit_grid_item \
  "hippocampus_LOCAL-RF_seed-included_no-global-prior" \
  "paper-rf/configs/hippocampus_local_rf_seed_included.toml"

# 6. LOCAL-RF seed-included, global prior from striatum seed-included DIFF-RF.
submit_grid_item \
  "hippocampus_LOCAL-RF_striatum-global-priors_seed-included" \
  "paper-rf/configs/hippocampus_local_rf_striatum_global_priors_seed_included.toml"

echo "All hippocampus sensitivity-grid jobs submitted."
