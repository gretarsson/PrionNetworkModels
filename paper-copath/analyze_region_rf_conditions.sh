#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v julia >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load julia
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PROJECT_DIR/.julia_depot}"

julia --project="$PROJECT_DIR" "$PROJECT_DIR/paper-copath/analyses/compare_region_rf_conditions.jl" \
  --project-root "$PROJECT_DIR"
