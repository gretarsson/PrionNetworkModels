#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v julia >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load julia
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PROJECT_DIR/.julia_depot}"

DATASETS=(
  "striatum:paper-rf/data/striatum/observations.csv:paper-rf/data/striatum/network.csv"
  "hippocampus:paper-rf/data/hippocampus/observations.csv:paper-rf/data/hippocampus/network.csv"
)

for spec in "${DATASETS[@]}"; do
  IFS=: read -r dataset observations network <<< "$spec"
  root="$PROJECT_DIR/runs/region_rf/${dataset}"
  echo "Plotting REGION-RF bundle for $dataset"
  julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/plot_region_rf_bundle.jl" \
    --root "$root" \
    --observations "$PROJECT_DIR/$observations" \
    --network "$PROJECT_DIR/$network"
done
