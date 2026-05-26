#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

DATASET="${COPATH_DATASET:-syn_mapt}"
CHAINS="${REGION_RF_CHAINS:-4}"
SAMPLES="${REGION_RF_SAMPLES:-1000}"
WARMUP="${REGION_RF_WARMUP:-1000}"
MAXITERS="${REGION_RF_MAXITERS:-10000}"
U0_PRIOR_SD="${REGION_RF_U0_PRIOR_SD:-0.01}"
ALPHA_PRIOR_SD="${REGION_RF_ALPHA_PRIOR_SD:-0.5}"
MAX_CONCURRENT="${REGION_RF_MAX_CONCURRENT:-40}"
SLURM_PARTITION="${SLURM_PARTITION:-all}"
SLURM_TIME="${SLURM_TIME:-2-00:00:00}"
SLURM_MEM="${SLURM_MEM:-8G}"
MEAN_DATA="${REGION_RF_MEAN_DATA:-0}"
SKIP_TRACES="${REGION_RF_SKIP_TRACES:-0}"

case "$DATASET" in
  syn_app)
    OBSERVATIONS="paper-copath/data/syn_pathology_app.csv"
    ;;
  syn_mapt)
    OBSERVATIONS="paper-copath/data/syn_pathology_mapt.csv"
    ;;
  tau_app)
    OBSERVATIONS="paper-copath/data/tau_pathology_app.csv"
    ;;
  tau_mapt)
    OBSERVATIONS="paper-copath/data/tau_pathology_mapt.csv"
    ;;
  *)
    echo "COPATH_DATASET must be one of: syn_app, syn_mapt, tau_app, tau_mapt." >&2
    exit 1
    ;;
esac

NETWORK="paper-copath/data/network.csv"
RUN_PREFIX="copath_${DATASET}_REGION-RF"
OUT_ROOT="runs/region_rf/copath_${DATASET}"
JOB_NAME="copath_${DATASET}_regionrf"

EXTRA_FLAGS=()
if [[ "$MEAN_DATA" == "1" ]]; then
  EXTRA_FLAGS+=("--mean-data")
fi
if [[ "$SKIP_TRACES" == "1" ]]; then
  EXTRA_FLAGS+=("--skip-traces")
fi
EXTRA_FLAGS_STR=""
if [[ ${#EXTRA_FLAGS[@]} -gt 0 ]]; then
  printf -v EXTRA_FLAGS_STR ' %q' "${EXTRA_FLAGS[@]}"
fi

echo "Submitting REGION-RF array for $DATASET."
echo "Output root: $PROJECT_DIR/$OUT_ROOT"
echo "Array: 1-412%$MAX_CONCURRENT, chains=$CHAINS, samples=$SAMPLES, warmup=$WARMUP"

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=$SLURM_MEM
#SBATCH --partition=$SLURM_PARTITION
#SBATCH --time=$SLURM_TIME
#SBATCH --array=1-412%$MAX_CONCURRENT
#SBATCH --chdir=$PROJECT_DIR
#SBATCH --output=$LOG_DIR/${JOB_NAME}_%A_%a.out
#SBATCH --error=$LOG_DIR/${JOB_NAME}_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail
module purge
module load julia
ulimit -t unlimited

export JULIA_DEPOT_PATH="$PROJECT_DIR/.julia_depot"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

exec julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/fit_region_rf.jl" \\
  --observations "$PROJECT_DIR/$OBSERVATIONS" \\
  --network "$PROJECT_DIR/$NETWORK" \\
  --region-index "\$SLURM_ARRAY_TASK_ID" \\
  --samples "$SAMPLES" \\
  --warmup "$WARMUP" \\
  --chains "$CHAINS" \\
  --maxiters "$MAXITERS" \\
  --u0-prior-sd "$U0_PRIOR_SD" \\
  --alpha-prior-sd "$ALPHA_PRIOR_SD" \\
  --run-prefix "$RUN_PREFIX" \\
  --out-root "$PROJECT_DIR/$OUT_ROOT" \\
  --no-root-summary$EXTRA_FLAGS_STR
EOF
