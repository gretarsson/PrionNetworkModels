#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-$PROJECT_DIR/paper-rf/python/.venv/bin/python}"

cd "$PROJECT_DIR"

"$PYTHON" paper-rf/analyses/figures/rebuild_model_figures.py
if [[ -f "$PROJECT_DIR/paper-rf/figures/pooled_z/all/pc1_direction/pc1_direction_comparison.csv" ]]; then
  "$PYTHON" paper-rf/analyses/rebuild_figures_6_7.py
  "$PYTHON" paper-rf/analyses/plotting/update_figure6_ai_panels.py
else
  echo "Skipping Figure 6/7 AI-panel refresh: pooled-z PC1 direction inputs are not present."
fi
