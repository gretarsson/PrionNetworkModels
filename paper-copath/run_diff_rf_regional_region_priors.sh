#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

"$PROJECT_DIR/cluster/prepare_julia_env.sh"

CHAINS="${DIFF_RF_REGIONAL_CHAINS:-4}"
SLURM_PARTITION="${SLURM_PARTITION:-all}"
SLURM_TIME="${SLURM_TIME:-2-00:00:00}"
PRIOR_SLURM_TIME="${DIFF_RF_REGIONAL_PRIOR_SLURM_TIME:-00:30:00}"
export SLURM_PARTITION SLURM_TIME

DATASETS=(
  "syn_app:paper-copath/configs/copath_syn_app_diff_rf_regional_region_priors.toml"
  "syn_mapt:paper-copath/configs/copath_syn_mapt_diff_rf_regional_region_priors.toml"
  "tau_app:paper-copath/configs/copath_tau_app_diff_rf_regional_region_priors.toml"
  "tau_mapt:paper-copath/configs/copath_tau_mapt_diff_rf_regional_region_priors.toml"
)

PRIOR_JOB_ID="$(
  sbatch --parsable <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=copath_diff_rf_regional_priors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=$SLURM_PARTITION
#SBATCH --time=$PRIOR_SLURM_TIME
#SBATCH --chdir=$PROJECT_DIR
#SBATCH --output=$LOG_DIR/copath_diff_rf_regional_priors-%j.out
#SBATCH --error=$LOG_DIR/copath_diff_rf_regional_priors-%j.err
#SBATCH --hint=nomultithread

set -euo pipefail
PROJECT_DIR="$PROJECT_DIR"
module purge
module load julia

export JULIA_DEPOT_PATH="\$PROJECT_DIR/.julia_depot"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "\$PROJECT_DIR/.julia_depot"

for dataset in syn_app syn_mapt tau_app tau_mapt; do
  julia --project="\$PROJECT_DIR" "\$PROJECT_DIR/paper-copath/analyses/make_region_rf_prior_source.jl" \\
    --dataset "\$dataset" \\
    --adjusted-region-rf-dir "\$PROJECT_DIR/paper-copath/results/region_rf_iterative_drop_low_likelihood_chains" \\
    --out "\$PROJECT_DIR/paper-copath/results/diff_rf_regional_priors/\${dataset}_region_rf_prior_source.h5"
done
EOF
)"
PRIOR_JOB_ID="${PRIOR_JOB_ID%%;*}"
PRIOR_DEPENDENCY="afterok:$PRIOR_JOB_ID"
echo "Submitted REGION-RF prior-source job $PRIOR_JOB_ID."

echo "Submitting copath DIFF-RF-REGIONAL jobs with REGION-RF posterior priors on partition '$SLURM_PARTITION' for $SLURM_TIME."
for item in "${DATASETS[@]}"; do
  dataset="${item%%:*}"
  config="${item#*:}"
  for chain in $(seq 1 "$CHAINS"); do
    run_id="copath_${dataset}_DIFF-RF-REGIONAL_region-priors_C${chain}"
    echo "Submitting $run_id from $config"
    SLURM_DEPENDENCY="$PRIOR_DEPENDENCY" "$PROJECT_DIR/cluster/submit_inference.sh" "$PROJECT_DIR/$config" "$run_id"
  done
done

echo "All copath DIFF-RF-REGIONAL region-prior jobs submitted after dependency $PRIOR_DEPENDENCY."
