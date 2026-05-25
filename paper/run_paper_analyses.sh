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

julia --project=. paper/analyses/model_parameters/export_parameter_tables.jl \
  --run "$STRIATUM_RUN" \
  --out-dir paper/results/parameters/striatum_diff_rf

julia --project=. paper/analyses/model_parameters/export_parameter_tables.jl \
  --run "$HIPPO_RUN" \
  --out-dir paper/results/parameters/hippocampus_diff_rf

"$PYTHON" paper/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper/results/parameters/striatum_diff_rf/beta.csv \
  --gamma paper/results/parameters/striatum_diff_rf/gamma.csv \
  --out-dir paper/results/transcriptomics/striatum

"$PYTHON" paper/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper/results/parameters/hippocampus_diff_rf/beta.csv \
  --gamma paper/results/parameters/hippocampus_diff_rf/gamma.csv \
  --out-dir paper/results/transcriptomics/hippocampus

"$PYTHON" paper/analyses/transcriptomics/compare_axes.py \
  --striatum-dir paper/results/transcriptomics/striatum \
  --hippocampus-dir paper/results/transcriptomics/hippocampus \
  --out-dir paper/results/transcriptomics/striatum_vs_hippocampus

if [[ "${RUN_GSEA:-0}" == "1" ]]; then
  "$PYTHON" paper/analyses/transcriptomics/run_gsea.py \
    --input paper/results/transcriptomics/striatum/gene_eta_correlations.csv \
    --out-dir paper/results/enrichment/striatum

  "$PYTHON" paper/analyses/transcriptomics/run_gsea.py \
    --input paper/results/transcriptomics/hippocampus/gene_eta_correlations.csv \
    --out-dir paper/results/enrichment/hippocampus
fi

"$PYTHON" paper/analyses/cell_types/cell_type_axis_associations.py \
  --axis paper/results/transcriptomics/striatum/region_axis.csv \
  --cell-types paper/data/cell_types/connectome_celltype.csv \
  --out-dir paper/results/cell_types/striatum

"$PYTHON" paper/analyses/cell_types/cell_type_axis_associations.py \
  --axis paper/results/transcriptomics/hippocampus/region_axis.csv \
  --cell-types paper/data/cell_types/connectome_celltype.csv \
  --out-dir paper/results/cell_types/hippocampus

"$PYTHON" paper/analyses/plotting/plot_biological_figures.py \
  --results-root paper/results \
  --out-dir paper/figures/biological

echo "Paper analyses complete:"
echo "  results: paper/results"
echo "  figures: paper/figures/biological"
