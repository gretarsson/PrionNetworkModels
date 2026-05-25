#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CHAINS="${PNM_SMOKE_CHAINS:-4}"
SAMPLES="${PNM_SMOKE_SAMPLES:-15}"
WARMUP="${PNM_SMOKE_WARMUP:-15}"
PREFIX="${PNM_SMOKE_PREFIX:-smoke-fit-local-rf-multichain}"
CONFIG="$PROJECT_DIR/configs/examples/local_rf.toml"

for CHAIN in $(seq 1 "$CHAINS"); do
  RUN_ID="${PREFIX}_C${CHAIN}"
  echo "LOCAL-RF smoke chain $CHAIN/$CHAINS -> $RUN_ID"
  julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/fit_model.jl" \
    --config "$CONFIG" \
    --run-id "$RUN_ID" \
    --samples "$SAMPLES" \
    --warmup "$WARMUP"
done

julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/merge_chains.jl" \
  --prefix "$PREFIX" \
  --out-run-id "$PREFIX" \
  --chain-count "$CHAINS"

julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/plot_run.jl" \
  --run "$PROJECT_DIR/runs/$PREFIX"
