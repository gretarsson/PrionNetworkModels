#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="$PROJECT_DIR/paper-rf/logs"
mkdir -p "$LOG_DIR"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

CHAINS="${REGION_RF_CHAINS:-4}"
SAMPLES="${REGION_RF_SAMPLES:-1000}"
WARMUP="${REGION_RF_WARMUP:-1000}"
MAXITERS="${REGION_RF_MAXITERS:-50000}"
U0_PRIOR_SD="${REGION_RF_U0_PRIOR_SD:-0.01}"
ALPHA_PRIOR_SD="${REGION_RF_ALPHA_PRIOR_SD:-1.0}"
MAX_CONCURRENT="${REGION_RF_MAX_CONCURRENT:-40}"
SLURM_PARTITION="${SLURM_PARTITION:-all}"
SLURM_TIME="${SLURM_TIME:-2-00:00:00}"
SLURM_MEM="${SLURM_MEM:-8G}"
MEAN_DATA="${REGION_RF_MEAN_DATA:-0}"
SKIP_TRACES="${REGION_RF_SKIP_TRACES:-0}"

DATASETS=(
  "striatum:paper-rf/data/striatum/observations.csv:paper-rf/data/striatum/network.csv"
  "hippocampus:paper-rf/data/hippocampus/observations.csv:paper-rf/data/hippocampus/network.csv"
)

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

echo "Submitting paper-rf REGION-RF arrays."
echo "Output root: $PROJECT_DIR/runs/region_rf"
echo "Array per dataset: 1-412%$MAX_CONCURRENT, chains=$CHAINS, samples=$SAMPLES, warmup=$WARMUP"

for spec in "${DATASETS[@]}"; do
  IFS=: read -r dataset observations network <<< "$spec"
  run_prefix="${dataset}_REGION-RF"
  out_root="runs/region_rf/${dataset}/regional_runs"
  job_name="${dataset}_regionrf"

  echo "Submitting $dataset -> $out_root"
  sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=$SLURM_MEM
#SBATCH --partition=$SLURM_PARTITION
#SBATCH --time=$SLURM_TIME
#SBATCH --array=1-412%$MAX_CONCURRENT
#SBATCH --chdir=$PROJECT_DIR
#SBATCH --output=$LOG_DIR/${job_name}_%A_%a.out
#SBATCH --error=$LOG_DIR/${job_name}_%A_%a.err
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
  --observations "$PROJECT_DIR/$observations" \\
  --network "$PROJECT_DIR/$network" \\
  --region-index "\$SLURM_ARRAY_TASK_ID" \\
  --samples "$SAMPLES" \\
  --warmup "$WARMUP" \\
  --chains "$CHAINS" \\
  --maxiters "$MAXITERS" \\
  --u0-prior-sd "$U0_PRIOR_SD" \\
  --alpha-prior-sd "$ALPHA_PRIOR_SD" \\
  --run-prefix "$run_prefix" \\
  --out-root "$PROJECT_DIR/$out_root" \\
  --no-root-summary$EXTRA_FLAGS_STR
EOF
done

echo "Submitted all paper-rf REGION-RF arrays."
