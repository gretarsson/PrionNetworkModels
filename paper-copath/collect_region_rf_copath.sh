#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATASET="${COPATH_DATASET:-syn_mapt}"

case "$DATASET" in
  syn_app|syn_mapt|tau_app|tau_mapt)
    ;;
  *)
    echo "COPATH_DATASET must be one of: syn_app, syn_mapt, tau_app, tau_mapt." >&2
    exit 1
    ;;
esac

ROOT="$PROJECT_DIR/runs/region_rf/copath_${DATASET}"

if ! command -v julia >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load julia
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PROJECT_DIR/.julia_depot}"

julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/collect_region_rf.jl" \
  --root "$ROOT"
