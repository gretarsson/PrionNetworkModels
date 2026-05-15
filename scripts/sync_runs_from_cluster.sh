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

mkdir -p "$LOCAL_RUNS_DIR"
mkdir -p "$LOCAL_LOGS_DIR"

echo "Syncing runs from $REMOTE_HOST:$REMOTE_RUNS_DIR/"
echo "Local destination: $LOCAL_RUNS_DIR/"

rsync -avP \
  -e "${SSH_CMD[*]}" \
  "$REMOTE_HOST:$REMOTE_RUNS_DIR/" \
  "$LOCAL_RUNS_DIR/"

echo "Syncing logs from $REMOTE_HOST:$REMOTE_LOGS_DIR/"
echo "Local destination: $LOCAL_LOGS_DIR/"

rsync -avP \
  -e "${SSH_CMD[*]}" \
  "$REMOTE_HOST:$REMOTE_LOGS_DIR/" \
  "$LOCAL_LOGS_DIR/"

echo "Run and log sync complete."
