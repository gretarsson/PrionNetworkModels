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

## Recommended Run Order

From the repository root:

```bash
julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/striatum_DIFF-RF_RETRO_C1_C3_C4 \
  --out-dir paper-rf/results/parameters/striatum_diff_rf

julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_partial_C1_C4 \
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

By default, this applies the manuscript-style filters `beta > 0` and
posterior-updated `beta`/`gamma` parameters (`ks_pvalue < 0.001`). For robustness
checks like Fig. S2, pass `--no-update-filter` or adjust `--beta-min`.

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

Filter-specific analyses are emitted under:

- `paper-rf/results/filtering/all/`
- `paper-rf/results/filtering/beta_positive/`
- `paper-rf/results/filtering/updated/`
- `paper-rf/figures/biological/filtering_levels/all/`
- `paper-rf/figures/biological/filtering_levels/beta_positive/`
- `paper-rf/figures/biological/filtering_levels/updated/`

Each filter-level folder contains striatal and hippocampal PCA outputs,
striatum/hippocampus comparison outputs, cell-type outputs, and GSEA outputs when
`RUN_GSEA=1` is used.

Filter-dependent supplementary comparison panels are emitted under each
filter-level figure folder:

- `comparison/` for gene-axis and gene-coefficient comparisons between striatal
  and hippocampal seeding

The raw regional `beta`/`gamma` comparison does not depend on the transcriptomic
filter and is emitted once under:

- `paper-rf/figures/biological/shared/`

## Provenance

This layer ports the relevant analysis logic from:

- `synuclein_spread`: original gene/parameter and PCA exploration
- `gene_enrichment`: pre-ranked GSEA and enrichment plotting
- `cell-type-atlas`: cell-type composition and monoaminergic association analyses

Only the paper-critical pieces should be kept here. General model functionality should
remain in the main package.
