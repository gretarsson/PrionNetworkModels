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
SYNC_RUN_SUBFOLDERS="${SYNC_RUN_SUBFOLDERS:-false}"
SYNC_LOGS="${SYNC_LOGS:-true}"

RUN_IDS=()

usage() {
  cat <<EOF
Usage: sync_runs_from_cluster.sh [options] [run_id ...]

By default, sync run bundle files directly under top-level runs/<run_id>/
directories without descending into nested subfolders. Pass one or more run IDs
to sync only those top-level run bundles.

Options:
  --recursive-runs    Sync the full remote runs/ tree, including subfolders.
  --no-logs           Skip syncing logs/.
  --logs-only         Sync logs/ only.
  -h, --help          Show this help.

Environment:
  SYNC_RUN_SUBFOLDERS=true   Same as --recursive-runs.
  SYNC_LOGS=false            Same as --no-logs.
EOF
}

SYNC_RUNS=true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --recursive-runs)
      SYNC_RUN_SUBFOLDERS=true
      shift
      ;;
    --no-logs)
      SYNC_LOGS=false
      shift
      ;;
    --logs-only)
      SYNC_RUNS=false
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        RUN_IDS+=("$1")
        shift
      done
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      RUN_IDS+=("$1")
      shift
      ;;
  esac
done

mkdir -p "$LOCAL_RUNS_DIR"
mkdir -p "$LOCAL_LOGS_DIR"

if [[ "$SYNC_RUNS" == "true" ]]; then
  echo "Syncing runs from $REMOTE_HOST:$REMOTE_RUNS_DIR/"
  echo "Local destination: $LOCAL_RUNS_DIR/"

  if [[ "$SYNC_RUN_SUBFOLDERS" == "true" ]]; then
    if [[ ${#RUN_IDS[@]} -gt 0 ]]; then
      echo "Ignoring explicit run IDs because --recursive-runs syncs the full runs/ tree." >&2
    fi
    rsync -avP \
      -e "${SSH_CMD[*]}" \
      "$REMOTE_HOST:$REMOTE_RUNS_DIR/" \
      "$LOCAL_RUNS_DIR/"
  else
    FILTERS=()
    if [[ ${#RUN_IDS[@]} -gt 0 ]]; then
      for RUN_ID in "${RUN_IDS[@]}"; do
        if [[ "$RUN_ID" == */* ]]; then
          echo "Run IDs must be top-level names under runs/: $RUN_ID" >&2
          exit 1
        fi
        FILTERS+=(--include="/$RUN_ID/")
        FILTERS+=(--exclude="/$RUN_ID/*/")
        FILTERS+=(--include="/$RUN_ID/*")
      done
    else
      FILTERS+=(--include='/*/')
      FILTERS+=(--exclude='/*/*/')
      FILTERS+=(--include='/*/*')
    fi
    FILTERS+=(--exclude='*')

    echo "Run sync mode: shallow top-level run bundles only."
    rsync -avP \
      -e "${SSH_CMD[*]}" \
      "${FILTERS[@]}" \
      "$REMOTE_HOST:$REMOTE_RUNS_DIR/" \
      "$LOCAL_RUNS_DIR/"
  fi
else
  echo "Skipping run sync."
fi

if [[ "$SYNC_LOGS" == "true" ]]; then
  echo "Syncing logs from $REMOTE_HOST:$REMOTE_LOGS_DIR/"
  echo "Local destination: $LOCAL_LOGS_DIR/"

  rsync -avP \
    -e "${SSH_CMD[*]}" \
    "$REMOTE_HOST:$REMOTE_LOGS_DIR/" \
    "$LOCAL_LOGS_DIR/"
else
  echo "Skipping log sync."
fi

echo "Sync complete."
