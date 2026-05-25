#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

CHAINS="${LOCAL_RF_CHAINS:-4}"
SLURM_PARTITION="${SLURM_PARTITION:-all}"
SLURM_TIME="${SLURM_TIME:-2-00:00:00}"
export SLURM_PARTITION SLURM_TIME

DATASETS=(
  "syn_app:paper-copath/configs/syn_app_local_rf.toml"
  "syn_mapt:paper-copath/configs/syn_mapt_local_rf.toml"
  "tau_app:paper-copath/configs/tau_app_local_rf.toml"
  "tau_mapt:paper-copath/configs/tau_mapt_local_rf.toml"
)

echo "Submitting syn-tau-abeta LOCAL-RF jobs on partition '$SLURM_PARTITION' for $SLURM_TIME."
for item in "${DATASETS[@]}"; do
  dataset="${item%%:*}"
  config="${item#*:}"
  for chain in $(seq 1 "$CHAINS"); do
    run_id="syn_tau_abeta_LOCAL-RF_${dataset}_C${chain}"
    echo "Submitting $run_id from $config"
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$config" "$run_id"
  done
done

echo "All syn-tau-abeta LOCAL-RF jobs submitted."
