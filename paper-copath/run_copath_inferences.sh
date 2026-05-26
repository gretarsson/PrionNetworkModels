#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

CHAINS="${LOCAL_RF_CHAINS:-4}"
SLURM_PARTITION="${SLURM_PARTITION:-all}"
SLURM_TIME="${SLURM_TIME:-2-00:00:00}"
COPATH_U0_MODE="${COPATH_U0_MODE:-inferu0}"
export SLURM_PARTITION SLURM_TIME

case "$COPATH_U0_MODE" in
  inferu0|detu0) ;;
  *)
    echo "COPATH_U0_MODE must be either 'inferu0' or 'detu0'." >&2
    exit 1
    ;;
esac

DATASETS=(
  "syn_app:paper-copath/configs/copath_syn_app_${COPATH_U0_MODE}.toml"
  "syn_mapt:paper-copath/configs/copath_syn_mapt_${COPATH_U0_MODE}.toml"
  "tau_app:paper-copath/configs/copath_tau_app_${COPATH_U0_MODE}.toml"
  "tau_mapt:paper-copath/configs/copath_tau_mapt_${COPATH_U0_MODE}.toml"
)

echo "Submitting copath LOCAL-RF ${COPATH_U0_MODE} jobs on partition '$SLURM_PARTITION' for $SLURM_TIME."
for item in "${DATASETS[@]}"; do
  dataset="${item%%:*}"
  config="${item#*:}"
  for chain in $(seq 1 "$CHAINS"); do
    run_id="copath_${dataset}_${COPATH_U0_MODE}_C${chain}"
    echo "Submitting $run_id from $config"
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$config" "$run_id"
  done
done

echo "All copath LOCAL-RF ${COPATH_U0_MODE} jobs submitted."
