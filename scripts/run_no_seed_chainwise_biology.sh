#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-$PROJECT_DIR/paper-rf/python/.venv/bin/python}"

if [[ ! -x "$PYTHON" ]]; then
  cat >&2 <<EOF
Missing paper Python environment:
  $PYTHON

Create it with:
  python3 -m venv paper-rf/python/.venv
  paper-rf/python/.venv/bin/python -m pip install -r paper-rf/python/requirements.txt
EOF
  exit 1
fi

cd "$PROJECT_DIR"

CHAINS="${CHAINS:-C1 C2 C3 C4}"
FILTER_LEVELS="${FILTER_LEVELS:-all beta_positive updated}"
RUN_GSEA="${RUN_GSEA:-1}"

RESULTS_ROOT="paper-rf/results/sidequests/no_seed_diff_rf/chainwise"
FIGURES_ROOT="paper-rf/figures/no_seed_diff_rf/chainwise"
EXPR="paper-rf/data/transcriptomics/avg_Pangea_exp.csv"
CELL_TYPES="paper-rf/data/cell_types/connectome_celltype.csv"

run_pca() {
  local beta="$1"
  local gamma="$2"
  local out_dir="$3"
  shift 3
  "$PYTHON" paper-rf/analyses/transcriptomics/gene_parameter_pca.py \
    --expression "$EXPR" \
    --beta "$beta" \
    --gamma "$gamma" \
    "$@" \
    --out-dir "$out_dir"
}

run_filter_level() {
  local chain="$1"
  local level="$2"
  shift 2
  local pca_args=("$@")
  local root="$RESULTS_ROOT/$chain/$level"
  local figures="$FIGURES_ROOT/$chain/$level"
  local striatum_params="$RESULTS_ROOT/$chain/parameters/striatum_diff_rf"
  local hippo_params="$RESULTS_ROOT/$chain/parameters/hippocampus_diff_rf"

  rm -rf "$root" "$figures"

  run_pca "$striatum_params/beta.csv" "$striatum_params/gamma.csv" \
    "$root/transcriptomics/striatum" ${pca_args[@]+"${pca_args[@]}"}

  run_pca "$hippo_params/beta.csv" "$hippo_params/gamma.csv" \
    "$root/transcriptomics/hippocampus" ${pca_args[@]+"${pca_args[@]}"}

  "$PYTHON" paper-rf/analyses/transcriptomics/compare_axes.py \
    --striatum-dir "$root/transcriptomics/striatum" \
    --hippocampus-dir "$root/transcriptomics/hippocampus" \
    --out-dir "$root/transcriptomics/striatum_vs_hippocampus"

  if [[ "$RUN_GSEA" == "1" ]]; then
    "$PYTHON" paper-rf/analyses/transcriptomics/run_gsea.py \
      --input "$root/transcriptomics/striatum/gene_eta_correlations.csv" \
      --out-dir "$root/enrichment/striatum"

    "$PYTHON" paper-rf/analyses/transcriptomics/run_gsea.py \
      --input "$root/transcriptomics/hippocampus/gene_eta_correlations.csv" \
      --out-dir "$root/enrichment/hippocampus"
  fi

  "$PYTHON" paper-rf/analyses/cell_types/cell_type_axis_associations.py \
    --axis "$root/transcriptomics/striatum/region_axis.csv" \
    --cell-types "$CELL_TYPES" \
    --out-dir "$root/cell_types/striatum"

  "$PYTHON" paper-rf/analyses/cell_types/cell_type_axis_associations.py \
    --axis "$root/transcriptomics/hippocampus/region_axis.csv" \
    --cell-types "$CELL_TYPES" \
    --out-dir "$root/cell_types/hippocampus"

  "$PYTHON" paper-rf/analyses/plotting/plot_biological_figures.py \
    --results-root "$root" \
    --out-dir "$figures"

  if [[ "$RUN_GSEA" == "1" ]]; then
    for dataset in striatum hippocampus; do
      mkdir -p "$figures/enrichment/$dataset/gseapy/prerank"
      mv "$figures/$dataset/gsea_dotplot_top_absNES."* "$figures/enrichment/$dataset/"
      mv "$figures/$dataset/category_fisher_log2odds."* "$figures/enrichment/$dataset/"
      mv "$figures/$dataset/category_nes_violin."* "$figures/enrichment/$dataset/"
      mv "$figures/$dataset/category_enrichment_fisher.csv" "$root/enrichment/$dataset/" 2>/dev/null || true
      find "$root/enrichment/$dataset/gseapy/prerank" -maxdepth 1 -type f -name '*.pdf' \
        -exec mv {} "$figures/enrichment/$dataset/gseapy/prerank/" \; 2>/dev/null || true
      rm -rf "$root/enrichment/$dataset/gseapy/prerank"
    done
  fi

  mkdir -p "$figures/transcriptomics" "$figures/cell_types"
  mv "$figures/comparison" "$figures/transcriptomics/" 2>/dev/null || true
  for dataset in striatum hippocampus; do
    mkdir -p "$figures/transcriptomics/$dataset" "$figures/cell_types/$dataset"
    mv "$figures/$dataset/pca_gene_coefficients."* "$figures/transcriptomics/$dataset/" 2>/dev/null || true
    mv "$figures/$dataset/beta_gamma_colored_by_eta."* "$figures/transcriptomics/$dataset/" 2>/dev/null || true
    mv "$figures/$dataset/celltype_eta_correlations."* "$figures/cell_types/$dataset/" 2>/dev/null || true
    mv "$figures/$dataset/monoaminergic_score_vs_eta."* "$figures/cell_types/$dataset/" 2>/dev/null || true
    rmdir "$figures/$dataset" 2>/dev/null || true
  done
}

for chain in $CHAINS; do
  striatum_run="runs/striatum_DIFF-RF_RETRO_ignore-seed_$chain"
  hippo_run="runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_ignore-seed_$chain"
  [[ -d "$striatum_run" ]] || { echo "Missing $striatum_run" >&2; exit 1; }
  [[ -d "$hippo_run" ]] || { echo "Missing $hippo_run" >&2; exit 1; }

  mkdir -p "$RESULTS_ROOT/$chain/parameters"
  rm -rf "$RESULTS_ROOT/$chain/parameters/striatum_diff_rf" "$RESULTS_ROOT/$chain/parameters/hippocampus_diff_rf"

  julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
    --run "$striatum_run" \
    --out-dir "$RESULTS_ROOT/$chain/parameters/striatum_diff_rf"

  julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
    --run "$hippo_run" \
    --out-dir "$RESULTS_ROOT/$chain/parameters/hippocampus_diff_rf"

  for level in $FILTER_LEVELS; do
    case "$level" in
      all)
        run_filter_level "$chain" "$level" --beta-min=-Inf --no-update-filter
        ;;
      beta_positive)
        run_filter_level "$chain" "$level" --no-update-filter
        ;;
      updated)
        run_filter_level "$chain" "$level"
        ;;
      *)
        echo "Unknown filter level: $level" >&2
        exit 1
        ;;
    esac
  done
done

"$PYTHON" - <<'PY'
from pathlib import Path
import pandas as pd

root = Path("paper-rf/results/sidequests/no_seed_diff_rf/chainwise")
rows = []
for chain_dir in sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("C")):
    for level_dir in sorted(p for p in chain_dir.iterdir() if p.is_dir() and p.name in {"all", "beta_positive", "updated"}):
        comp = level_dir / "transcriptomics/striatum_vs_hippocampus/axis_comparison_summary.csv"
        if not comp.exists():
            continue
        stats = pd.read_csv(comp).iloc[0].to_dict()
        for dataset in ["striatum", "hippocampus"]:
            pca = pd.read_csv(level_dir / f"transcriptomics/{dataset}/pca_summary.csv")
            pc1 = pca[pca["component"] == "PC1"].iloc[0]
            rows.append({
                "chain": chain_dir.name,
                "filter_level": level_dir.name,
                "dataset": dataset,
                "pc1_loading_beta": pc1["loading_beta"],
                "pc1_loading_gamma": pc1["loading_gamma"],
                "pc1_explained_variance": pc1["explained_variance_ratio"],
                "n_regions": int(pc1["n_regions"]),
                "n_genes": int(pc1["n_genes"]),
                "striatum_hippocampus_pc_cosine": stats.get("pc_cosine_similarity"),
                "striatum_hippocampus_pc_angle_degrees": stats.get("pc_absolute_angle_degrees"),
                "gene_eta_pearson_r": stats.get("gene_eta_pearson_r"),
                "gene_eta_sign_agreement": stats.get("gene_eta_sign_agreement"),
            })
out = pd.DataFrame(rows)
out_path = root / "chainwise_axis_summary.csv"
out.to_csv(out_path, index=False)
print(out_path)
PY

echo "Chainwise no-seed biology complete:"
echo "  results: $RESULTS_ROOT"
echo "  figures: $FIGURES_ROOT"
