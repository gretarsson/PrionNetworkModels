#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v julia >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load julia
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PROJECT_DIR/.julia_depot}"

DATASETS=(
  "syn_app:paper-copath/data/syn_pathology_app.csv"
  "syn_mapt:paper-copath/data/syn_pathology_mapt.csv"
  "tau_app:paper-copath/data/tau_pathology_app.csv"
  "tau_mapt:paper-copath/data/tau_pathology_mapt.csv"
)

for spec in "${DATASETS[@]}"; do
  dataset="${spec%%:*}"
  observations="${spec#*:}"
  root="$PROJECT_DIR/runs/region_rf/copath_${dataset}"
  echo "Plotting REGION-RF bundle for $dataset"
  julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/plot_region_rf_bundle.jl" \
    --root "$root" \
    --observations "$PROJECT_DIR/$observations" \
    --network "$PROJECT_DIR/paper-copath/data/network.csv"
done
