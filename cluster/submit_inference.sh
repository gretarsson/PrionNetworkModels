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
mkdir -p "$LOG_DIR"

JOB_NAME="${RUN_ID}"
EXTRA_ARGS=("$@")
EXTRA_ARGS_STR=""
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  printf -v EXTRA_ARGS_STR ' %q' "${EXTRA_ARGS[@]}"
fi

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --time=2-00:00:00
#SBATCH --chdir=$PROJECT_DIR
#SBATCH --output=$LOG_DIR/${JOB_NAME}-%j.out
#SBATCH --error=$LOG_DIR/${JOB_NAME}-%j.err
#SBATCH --hint=nomultithread

set -euo pipefail
PROJECT_DIR="$PROJECT_DIR"
module purge
module load julia

export JULIA_DEPOT_PATH="$PROJECT_DIR/.julia_depot"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "$PROJECT_DIR/.julia_depot"

julia --project="$PROJECT_DIR" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

exec julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/fit_model.jl" \
  --config "$CONFIG_PATH" \
  --run-id "$RUN_ID"$EXTRA_ARGS_STR
EOF
