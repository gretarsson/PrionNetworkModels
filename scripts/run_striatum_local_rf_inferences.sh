#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

JOBNAME="striatum_LOCAL-RF"
CONFIG="configs/paper/striatum_local_rf_core.toml"
LOCAL_RF_CHAINS="${LOCAL_RF_CHAINS:-4}"

echo "Submitting $LOCAL_RF_CHAINS striatum LOCAL-RF chains."
for CHAIN in $(seq 1 "$LOCAL_RF_CHAINS"); do
  RUN_ID="${JOBNAME}_C${CHAIN}"
  echo "Submitting $RUN_ID from $CONFIG"
  "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$CONFIG" "$RUN_ID"
done

echo "All striatum LOCAL-RF jobs submitted."
