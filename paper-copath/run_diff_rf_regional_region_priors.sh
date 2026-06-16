#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

if command -v module >/dev/null 2>&1; then
  module purge
  module load julia
fi

if ! command -v julia >/dev/null 2>&1; then
  echo "Could not find julia on PATH after loading the cluster Julia module." >&2
  exit 1
fi

CHAINS="${DIFF_RF_REGIONAL_CHAINS:-4}"
SLURM_PARTITION="${SLURM_PARTITION:-all}"
SLURM_TIME="${SLURM_TIME:-2-00:00:00}"
export SLURM_PARTITION SLURM_TIME
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PROJECT_DIR/.julia_depot}"

DATASETS=(
  "syn_app:paper-copath/configs/copath_syn_app_diff_rf_regional_region_priors.toml"
  "syn_mapt:paper-copath/configs/copath_syn_mapt_diff_rf_regional_region_priors.toml"
  "tau_app:paper-copath/configs/copath_tau_app_diff_rf_regional_region_priors.toml"
  "tau_mapt:paper-copath/configs/copath_tau_mapt_diff_rf_regional_region_priors.toml"
)

for item in "${DATASETS[@]}"; do
  dataset="${item%%:*}"
  julia --project="$PROJECT_DIR" "$PROJECT_DIR/paper-copath/analyses/make_region_rf_prior_source.jl" \
    --dataset "$dataset" \
    --adjusted-region-rf-dir "$PROJECT_DIR/paper-copath/results/region_rf_iterative_drop_low_likelihood_chains" \
    --out "$PROJECT_DIR/paper-copath/results/diff_rf_regional_priors/${dataset}_region_rf_prior_source.h5"
done

echo "Submitting copath DIFF-RF-REGIONAL jobs with REGION-RF posterior priors on partition '$SLURM_PARTITION' for $SLURM_TIME."
for item in "${DATASETS[@]}"; do
  dataset="${item%%:*}"
  config="${item#*:}"
  for chain in $(seq 1 "$CHAINS"); do
    run_id="copath_${dataset}_DIFF-RF-REGIONAL_region-priors_C${chain}"
    echo "Submitting $run_id from $config"
    "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$config" "$run_id"
  done
done

echo "All copath DIFF-RF-REGIONAL region-prior jobs submitted."
