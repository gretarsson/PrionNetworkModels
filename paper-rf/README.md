# Paper Reproduction Layer

This directory contains the paper-specific analyses that sit on top of the reusable
`PrionNetworkModels` package.

The core package should stay general: model definitions, inference, run bundles,
merging, and model diagnostics live in `src/`, `scripts/`, `configs/`, and `cluster/`.
Analyses that exist to reproduce this manuscript live here.

## Main Analysis Targets

The current manuscript requires these paper-specific analyses:

- posterior parameter export for selected DIFF-RF run bundles
- gene-expression associations with rise (`beta`) and fall (`gamma`) parameters
- PCA of gene-level `beta`/`gamma` regression coefficients
- pre-ranked KEGG/GSEA enrichment along the PCA vulnerability axis
- cell-type and monoaminergic-score associations with the same axis
- striatum/hippocampus comparison of gene-parameter associations

## Inputs

Curated paper input tables are stored in:

- `paper-rf/data/transcriptomics/avg_Pangea_exp.csv`
- `paper-rf/data/cell_types/connectome_celltype.csv`

The modeling inputs remain in `paper-rf/data/`. Large generated run bundles remain in
`runs/` and are not part of the source package API.

The manuscript figures are regenerated from fitted posterior bundles. These are
large artifacts and are intentionally not tracked in git. After downloading the
archived inference outputs, place them under `runs/` with these names:

```text
runs/striatum_DIFF-RF_RETRO_paper/
runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4/
```

## Recommended Run Order

From the repository root:

```bash
julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/striatum_DIFF-RF_RETRO_paper \
  --out-dir paper-rf/results/parameters/striatum_diff_rf

julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4 \
  --out-dir paper-rf/results/parameters/hippocampus_diff_rf
```

Then create gene-parameter PCA outputs:

```bash
python paper-rf/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper-rf/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper-rf/results/parameters/striatum_diff_rf/beta.csv \
  --gamma paper-rf/results/parameters/striatum_diff_rf/gamma.csv \
  --out-dir paper-rf/results/transcriptomics/striatum

python paper-rf/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper-rf/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper-rf/results/parameters/hippocampus_diff_rf/beta.csv \
  --gamma paper-rf/results/parameters/hippocampus_diff_rf/gamma.csv \
  --out-dir paper-rf/results/transcriptomics/hippocampus
```

By default, the figure workflow uses all transcriptomics-matched regions. The
PCA scripts also expose filtering options such as `--beta-min` for sensitivity
checks.

Optional GSEA requires `gseapy` and internet/cache access to the Enrichr library:

```bash
python paper-rf/analyses/transcriptomics/run_gsea.py \
  --input paper-rf/results/transcriptomics/striatum/gene_eta_correlations.csv \
  --out-dir paper-rf/results/enrichment/striatum
```

Cell-type associations:

```bash
python paper-rf/analyses/cell_types/cell_type_axis_associations.py \
  --axis paper-rf/results/transcriptomics/striatum/region_axis.csv \
  --cell-types paper-rf/data/cell_types/connectome_celltype.csv \
  --out-dir paper-rf/results/cell_types/striatum
```

Compare striatal and hippocampal gene-axis structure:

```bash
python paper-rf/analyses/transcriptomics/compare_axes.py \
  --striatum-dir paper-rf/results/transcriptomics/striatum \
  --hippocampus-dir paper-rf/results/transcriptomics/hippocampus \
  --out-dir paper-rf/results/transcriptomics/striatum_vs_hippocampus
```

Create manuscript-style biological figure panels:

```bash
python paper-rf/analyses/plotting/plot_biological_figures.py \
  --results-root paper-rf/results \
  --out-dir paper-rf/figures/biological
```

To regenerate the final manuscript panels for Figures 6 and 7:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/rebuild_figures_6_7.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_figure6_ai_panels.py
```

This writes the independent PDF panels used for manual composition under:

```text
paper-rf/figures/Figure6/
paper-rf/figures/Figure7/
```

To rebuild the current independent panel folders for the main scientific figures,
run:

```bash
bash paper-rf/run_main_figure_panels.sh
```

This writes:

- `paper-rf/figures/Figure2/`: transport WAIC panel.
- `paper-rf/figures/Figure3/`: model WAIC, run-generated predicted-versus-observed panels, run-generated top-4 retrodiction panels, and timepoint agreement panel.
- `paper-rf/figures/Figure4/`: null-model WAIC panel.
- `paper-rf/figures/Figure5/`: held-out predicted-versus-observed and top-4 retrodiction panels.
- `paper-rf/figures/Figure6/` and `paper-rf/figures/Figure7/`: vulnerability-axis panels.

The Figure 2, Figure 4, and Figure 5 source inferences can be imported from the
original `synuclein_spread` artifacts with:

```bash
julia --project=/Users/gretarsson/Desktop/synuclein_spread \
  paper-rf/analyses/figures/import_synuclein_spread_runs.jl
julia --project=. paper-rf/analyses/figures/export_run_bundle_waic.jl
```

After those commands, `paper-rf/analyses/figures/rebuild_model_figures.py`
regenerates the Figure 2-5 scientific panels from `runs/` bundles and run-derived
WAIC tables. The `synuclein_spread` environment is only needed for the one-time
translation/import step; plotting then reads from this repository. Each figure
folder contains a `missing_requirements.md` file; it
should say `None.` when all inputs for that figure are available.

The end-to-end convenience wrapper runs the parameter export, transcriptomics,
cell-type analyses, axis comparison, and plotting separately for three region
filters:

- `all`: all transcriptomics-matched regions, no posterior-update filter
- `beta_positive`: `beta > 0`, no posterior-update filter
- `updated`: `beta > 0` and both `beta`/`gamma` posterior-updated with
  `ks_pvalue < 0.001`

```bash
bash paper-rf/run_paper_analyses.sh
```

To also regenerate full KEGG/GSEA outputs before plotting, run:

```bash
RUN_GSEA=1 bash paper-rf/run_paper_analyses.sh
```

By default, the hippocampus analysis uses the DIFF-RF run initialized with
striatal posterior-derived priors for global parameters (`rho`, `alpha`, and
`sigma`). To use the normal hippocampus merge instead, override `HIPPO_RUN`.

Generated tables are written under `paper-rf/results/`, and generated figures are
written under `paper-rf/figures/`. Both directories are ignored by git.

Appendix vulnerability-axis panels are regenerated with:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/build_appendix_vulnerability_inputs.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_appendix_vulnerability_figures.py
```

Only paper-critical analysis code should live here. General model functionality
should remain in the main package.
