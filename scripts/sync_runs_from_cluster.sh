#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_RUNS_DIR="${LOCAL_RUNS_DIR:-$PROJECT_DIR/runs}"

REMOTE_HOST="${REMOTE_HOST:-alexanderc@cubic-login.uphs.upenn.edu}"
REMOTE_PROJECT_DIR="${REMOTE_PROJECT_DIR:-~/PrionNetworkModels}"
REMOTE_RUNS_DIR="${REMOTE_RUNS_DIR:-$REMOTE_PROJECT_DIR/runs}"

mkdir -p "$LOCAL_RUNS_DIR"

echo "Syncing runs from $REMOTE_HOST:$REMOTE_RUNS_DIR/"
echo "Local destination: $LOCAL_RUNS_DIR/"

rsync -avP \
  "$REMOTE_HOST:$REMOTE_RUNS_DIR/" \
  "$LOCAL_RUNS_DIR/"

echo "Run sync complete."
