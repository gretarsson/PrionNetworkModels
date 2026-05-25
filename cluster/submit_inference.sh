#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: submit_inference.sh <config.toml> <run_id> [extra fit_model args...]"
  exit 1
fi

CONFIG_PATH="$1"
RUN_ID="$2"
shift 2

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
MARKER_FILE="$PROJECT_DIR/.julia_depot/.prepared_prionnetworkmodels"
mkdir -p "$LOG_DIR"

if [[ ! -f "$MARKER_FILE" || "$MARKER_FILE" -ot "$PROJECT_DIR/Project.toml" || ( -f "$PROJECT_DIR/Manifest.toml" && "$MARKER_FILE" -ot "$PROJECT_DIR/Manifest.toml" ) ]]; then
  echo "Preparing Julia environment before submission..."
  "$PROJECT_DIR/cluster/prepare_julia_env.sh"
fi

JOB_NAME="${RUN_ID}"
EXTRA_ARGS=("$@")
EXTRA_ARGS_STR=""
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  printf -v EXTRA_ARGS_STR ' %q' "${EXTRA_ARGS[@]}"
fi
SBATCH_PARTITION="${SLURM_PARTITION:-all}"
SBATCH_TIME="${SLURM_TIME:-2-00:00:00}"

submit_job() {
  sbatch "$@" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=$SBATCH_PARTITION
#SBATCH --time=$SBATCH_TIME
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

exec julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/fit_model.jl" \
  --config "$CONFIG_PATH" \
  --run-id "$RUN_ID"$EXTRA_ARGS_STR
EOF
}

if [[ -n "${SLURM_DEPENDENCY:-}" ]]; then
  submit_job --dependency="$SLURM_DEPENDENCY"
else
  submit_job
fi
