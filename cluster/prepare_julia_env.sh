#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEPOT_DIR="$PROJECT_DIR/.julia_depot"
MARKER_FILE="$DEPOT_DIR/.prepared_prionnetworkmodels"

export JULIA_DEPOT_PATH="$DEPOT_DIR"
mkdir -p "$DEPOT_DIR"

if ! command -v julia >/dev/null 2>&1; then
  if command -v module >/dev/null 2>&1; then
    module purge
    module load julia
  fi
fi

if ! command -v julia >/dev/null 2>&1; then
  echo "Could not find julia on PATH while preparing the project environment."
  exit 1
fi

needs_prepare=false
if [[ ! -f "$MARKER_FILE" ]]; then
  needs_prepare=true
elif [[ "$MARKER_FILE" -ot "$PROJECT_DIR/Project.toml" ]]; then
  needs_prepare=true
elif [[ -f "$PROJECT_DIR/Manifest.toml" && "$MARKER_FILE" -ot "$PROJECT_DIR/Manifest.toml" ]]; then
  needs_prepare=true
fi

if [[ "$needs_prepare" == false ]]; then
  echo "Julia environment already prepared in $DEPOT_DIR"
  exit 0
fi

echo "Preparing Julia environment in $DEPOT_DIR"
julia --project="$PROJECT_DIR" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using PrionNetworkModels'
touch "$MARKER_FILE"
echo "Julia environment ready."
