#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_RUNS_DIR="${LOCAL_RUNS_DIR:-$PROJECT_DIR/runs}"
LOCAL_LOGS_DIR="${LOCAL_LOGS_DIR:-$PROJECT_DIR/logs}"

REMOTE_HOST="${REMOTE_HOST:-alexanderc@cubic-login.uphs.upenn.edu}"
REMOTE_PROJECT_DIR="${REMOTE_PROJECT_DIR:-~/PrionNetworkModels}"
REMOTE_RUNS_DIR="${REMOTE_RUNS_DIR:-$REMOTE_PROJECT_DIR/runs}"
REMOTE_LOGS_DIR="${REMOTE_LOGS_DIR:-$REMOTE_PROJECT_DIR/logs}"
CONTROL_SOCKET="${CONTROL_SOCKET:-/tmp/prionnetworkmodels-rsync-%r@%h:%p}"
SSH_CMD=(ssh -o ControlMaster=auto -o "ControlPath=$CONTROL_SOCKET" -o ControlPersist=2m)
RSYNC_CMD=(rsync -avP -e "${SSH_CMD[*]}")
FULL_SYNC="${FULL_SYNC:-0}"

mkdir -p "$LOCAL_RUNS_DIR"
mkdir -p "$LOCAL_LOGS_DIR"

echo "Syncing runs from $REMOTE_HOST:$REMOTE_RUNS_DIR/"
echo "Local destination: $LOCAL_RUNS_DIR/"

if [[ "$FULL_SYNC" == "1" ]]; then
  "${RSYNC_CMD[@]}" \
    "$REMOTE_HOST:$REMOTE_RUNS_DIR/" \
    "$LOCAL_RUNS_DIR/"
else
  synced=0
  skipped=0
  while IFS= read -r remote_dir; do
    rel_dir="${remote_dir#./}"
    local_dir="$LOCAL_RUNS_DIR/$rel_dir"
    if [[ -d "$local_dir" ]]; then
      skipped=$((skipped + 1))
      continue
    fi

    mkdir -p "$(dirname "$local_dir")"
    echo "Syncing new run folder: $rel_dir/"
    "${RSYNC_CMD[@]}" \
      "$REMOTE_HOST:$REMOTE_RUNS_DIR/$rel_dir/" \
      "$local_dir/"
    synced=$((synced + 1))
  done < <(
    "${SSH_CMD[@]}" "$REMOTE_HOST" \
      "cd $REMOTE_RUNS_DIR && find . -mindepth 1 -type d | sort"
  )

  echo "Run folder sync complete: $synced new folder(s), $skipped already present."
  echo "Set FULL_SYNC=1 to force the old full rsync behavior."
fi

echo "Syncing logs from $REMOTE_HOST:$REMOTE_LOGS_DIR/"
echo "Local destination: $LOCAL_LOGS_DIR/"

"${RSYNC_CMD[@]}" \
  --ignore-existing \
  "$REMOTE_HOST:$REMOTE_LOGS_DIR/" \
  "$LOCAL_LOGS_DIR/"

echo "Run and log sync complete."
