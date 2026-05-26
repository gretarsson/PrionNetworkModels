#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATASETS=(striatum hippocampus)

if ! command -v julia >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load julia
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PROJECT_DIR/.julia_depot}"

for dataset in "${DATASETS[@]}"; do
  root="$PROJECT_DIR/paper-rf/results/region_rf/${dataset}"
  if [[ ! -d "$root" ]]; then
    echo "Skipping $dataset; output directory does not exist: $root"
    continue
  fi
  echo "Collecting $dataset"
  julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/collect_region_rf.jl" \
    --root "$root"
done
