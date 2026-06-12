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

STRIATUM_CHAINS="${STRIATUM_CHAINS:-C1 C2 C3 C4}"
HIPPOCAMPUS_CHAINS="${HIPPOCAMPUS_CHAINS:-C1 C2 C3 C4 C5 C6 C7 C8}"
FILTER_LEVELS="${FILTER_LEVELS:-all beta_positive updated}"

RESULTS_ROOT="paper-rf/results/sidequests/seed_fit_diff_rf/chainwise"
FIGURES_ROOT="paper-rf/figures/seed_fit_diff_rf/chainwise"
EXPR="paper-rf/data/transcriptomics/avg_Pangea_exp.csv"

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

run_chain_dataset() {
  local dataset="$1"
  local chain="$2"
  local run_dir="$3"
  local params_dir="$RESULTS_ROOT/$chain/parameters/${dataset}_diff_rf"

  [[ -d "$run_dir" ]] || { echo "Missing $run_dir" >&2; return 0; }

  rm -rf "$params_dir"
  julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
    --run "$run_dir" \
    --out-dir "$params_dir"

  for level in $FILTER_LEVELS; do
    local out_dir="$RESULTS_ROOT/$chain/$level/transcriptomics/$dataset"
    rm -rf "$out_dir"
    case "$level" in
      all)
        run_pca "$params_dir/beta.csv" "$params_dir/gamma.csv" "$out_dir" --beta-min=-Inf --no-update-filter
        ;;
      beta_positive)
        run_pca "$params_dir/beta.csv" "$params_dir/gamma.csv" "$out_dir" --no-update-filter
        ;;
      updated)
        run_pca "$params_dir/beta.csv" "$params_dir/gamma.csv" "$out_dir"
        ;;
      *)
        echo "Unknown filter level: $level" >&2
        exit 1
        ;;
    esac
  done
}

for chain in $STRIATUM_CHAINS; do
  run_chain_dataset "striatum" "$chain" "runs/striatum_DIFF-RF_RETRO_$chain"
done

for chain in $HIPPOCAMPUS_CHAINS; do
  run_chain_dataset "hippocampus" "$chain" "runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_$chain"
done

"$PYTHON" - <<'PY'
from pathlib import Path
import pandas as pd

root = Path("paper-rf/results/sidequests/seed_fit_diff_rf/chainwise")
rows = []
for chain_dir in sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("C")):
    for level_dir in sorted(p for p in chain_dir.iterdir() if p.is_dir() and p.name in {"all", "beta_positive", "updated"}):
        for dataset in ["striatum", "hippocampus"]:
            pca_path = level_dir / f"transcriptomics/{dataset}/pca_summary.csv"
            if not pca_path.exists():
                continue
            pca = pd.read_csv(pca_path)
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
            })
out = pd.DataFrame(rows)
out_path = root / "chainwise_axis_summary.csv"
out.to_csv(out_path, index=False)
print(out_path)
PY

"$PYTHON" paper-rf/analyses/plotting/plot_no_seed_chainwise_pc1_directions.py \
  --summary "$RESULTS_ROOT/chainwise_axis_summary.csv" \
  --out "$FIGURES_ROOT/pc1_direction_by_filter"

echo "Seed-fit chainwise PC1 directions complete:"
echo "  results: $RESULTS_ROOT"
echo "  figures: $FIGURES_ROOT"
