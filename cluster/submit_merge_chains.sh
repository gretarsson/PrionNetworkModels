#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: submit_merge_chains.sh <run_prefix> [out_run_id] [chain_count] [runs_root]" >&2
  exit 1
fi

RUN_PREFIX="$1"
OUT_RUN_ID="${2:-$RUN_PREFIX}"
CHAIN_COUNT="${3:-4}"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNS_ROOT="${4:-$PROJECT_DIR/runs}"
LOG_DIR="$PROJECT_DIR/logs"
MARKER_FILE="$PROJECT_DIR/.julia_depot/.prepared_prionnetworkmodels"
mkdir -p "$LOG_DIR"

if [[ ! -f "$MARKER_FILE" || "$MARKER_FILE" -ot "$PROJECT_DIR/Project.toml" || ( -f "$PROJECT_DIR/Manifest.toml" && "$MARKER_FILE" -ot "$PROJECT_DIR/Manifest.toml" ) ]]; then
  echo "Preparing Julia environment before merge submission..." >&2
  "$PROJECT_DIR/cluster/prepare_julia_env.sh" >&2
fi

JOB_NAME="merge_${OUT_RUN_ID}"
RAW_JOB_ID="$(
  sbatch --parsable <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --time=04:00:00
#SBATCH --chdir=$PROJECT_DIR
#SBATCH --output=$LOG_DIR/${JOB_NAME}-%j.out
#SBATCH --error=$LOG_DIR/${JOB_NAME}-%j.err
#SBATCH --hint=nomultithread

set -euo pipefail
PROJECT_DIR="$PROJECT_DIR"
module purge
module load julia
ulimit -t unlimited

export JULIA_DEPOT_PATH="$PROJECT_DIR/.julia_depot"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "$PROJECT_DIR/.julia_depot"

exec julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/merge_chains.jl" \
  --prefix "$RUN_PREFIX" \
  --out-run-id "$OUT_RUN_ID" \
  --chain-count "$CHAIN_COUNT" \
  --runs-root "$RUNS_ROOT"
EOF
)"
JOB_ID="${RAW_JOB_ID%%;*}"

echo "Submitted merge job $JOB_ID for $RUN_PREFIX -> $OUT_RUN_ID" >&2
echo "$JOB_ID"
