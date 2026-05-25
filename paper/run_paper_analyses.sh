#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-$PROJECT_DIR/paper/python/.venv/bin/python}"

if [[ ! -x "$PYTHON" ]]; then
  cat >&2 <<EOF
Missing paper Python environment:
  $PYTHON

Create it with:
  python3 -m venv paper/python/.venv
  paper/python/.venv/bin/python -m pip install -r paper/python/requirements.txt
EOF
  exit 1
fi

cd "$PROJECT_DIR"

STRIATUM_RUN="${STRIATUM_RUN:-runs/striatum_DIFF-RF_RETRO_C1_C3_C4}"
HIPPO_RUN="${HIPPO_RUN:-runs/hippocampus_DIFF-RF_RETRO_C1_C2_C3}"

STRIATUM_PARAMS="paper/results/parameters/striatum_diff_rf"
HIPPO_PARAMS="paper/results/parameters/hippocampus_diff_rf"

julia --project=. paper/analyses/model_parameters/export_parameter_tables.jl \
  --run "$STRIATUM_RUN" \
  --out-dir "$STRIATUM_PARAMS"

julia --project=. paper/analyses/model_parameters/export_parameter_tables.jl \
  --run "$HIPPO_RUN" \
  --out-dir "$HIPPO_PARAMS"

run_filter_level() {
  local level="$1"
  shift
  local pca_args=("$@")
  local root="paper/results/filtering/$level"
  local figures="paper/figures/biological/filtering_levels/$level"

  rm -rf "$root/transcriptomics" "$root/cell_types" "$root/enrichment" "$figures"
  mkdir -p "$root/parameters"
  rm -rf "$root/parameters/striatum_diff_rf" "$root/parameters/hippocampus_diff_rf"
  cp -R "$STRIATUM_PARAMS" "$root/parameters/striatum_diff_rf"
  cp -R "$HIPPO_PARAMS" "$root/parameters/hippocampus_diff_rf"

  "$PYTHON" paper/analyses/transcriptomics/gene_parameter_pca.py \
    --expression paper/data/transcriptomics/avg_Pangea_exp.csv \
    --beta "$STRIATUM_PARAMS/beta.csv" \
    --gamma "$STRIATUM_PARAMS/gamma.csv" \
    ${pca_args[@]+"${pca_args[@]}"} \
    --out-dir "$root/transcriptomics/striatum"

  "$PYTHON" paper/analyses/transcriptomics/gene_parameter_pca.py \
    --expression paper/data/transcriptomics/avg_Pangea_exp.csv \
    --beta "$HIPPO_PARAMS/beta.csv" \
    --gamma "$HIPPO_PARAMS/gamma.csv" \
    ${pca_args[@]+"${pca_args[@]}"} \
    --out-dir "$root/transcriptomics/hippocampus"

  "$PYTHON" paper/analyses/transcriptomics/compare_axes.py \
    --striatum-dir "$root/transcriptomics/striatum" \
    --hippocampus-dir "$root/transcriptomics/hippocampus" \
    --out-dir "$root/transcriptomics/striatum_vs_hippocampus"

  if [[ "${RUN_GSEA:-0}" == "1" ]]; then
    "$PYTHON" paper/analyses/transcriptomics/run_gsea.py \
      --input "$root/transcriptomics/striatum/gene_eta_correlations.csv" \
      --out-dir "$root/enrichment/striatum"

    "$PYTHON" paper/analyses/transcriptomics/run_gsea.py \
      --input "$root/transcriptomics/hippocampus/gene_eta_correlations.csv" \
      --out-dir "$root/enrichment/hippocampus"
  fi

  "$PYTHON" paper/analyses/cell_types/cell_type_axis_associations.py \
    --axis "$root/transcriptomics/striatum/region_axis.csv" \
    --cell-types paper/data/cell_types/connectome_celltype.csv \
    --out-dir "$root/cell_types/striatum"

  "$PYTHON" paper/analyses/cell_types/cell_type_axis_associations.py \
    --axis "$root/transcriptomics/hippocampus/region_axis.csv" \
    --cell-types paper/data/cell_types/connectome_celltype.csv \
    --out-dir "$root/cell_types/hippocampus"

  "$PYTHON" paper/analyses/plotting/plot_biological_figures.py \
    --results-root "$root" \
    --out-dir "$figures"
}

for level in ${FILTER_LEVELS:-all beta_positive updated}; do
  case "$level" in
    all)
      run_filter_level "all" --beta-min=-Inf --no-update-filter
      ;;
    beta_positive)
      run_filter_level "beta_positive" --no-update-filter
      ;;
    updated)
      run_filter_level "updated"
      ;;
    *)
      echo "Unknown filter level: $level" >&2
      exit 1
      ;;
  esac
done

echo "Paper analyses complete:"
echo "  shared parameters: paper/results/parameters"
echo "  filter-level results: paper/results/filtering/{all,beta_positive,updated}"
echo "  filter-level figures: paper/figures/biological/filtering_levels/{all,beta_positive,updated}"
