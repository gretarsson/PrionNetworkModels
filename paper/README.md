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

- `paper/data/transcriptomics/avg_Pangea_exp.csv`
- `paper/data/cell_types/connectome_celltype.csv`

The modeling inputs remain in `data/paper/`. Large generated run bundles remain in
`runs/` and are not part of the source package API.

## Recommended Run Order

From the repository root:

```bash
julia --project=. paper/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/striatum_DIFF-RF_RETRO_C1_C3_C4 \
  --out-dir paper/results/parameters/striatum_diff_rf

julia --project=. paper/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_partial_C1_C4 \
  --out-dir paper/results/parameters/hippocampus_diff_rf_striatum_global_priors
```

Then create gene-parameter PCA outputs:

```bash
python paper/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper/results/parameters/striatum_diff_rf/beta.csv \
  --gamma paper/results/parameters/striatum_diff_rf/gamma.csv \
  --out-dir paper/results/transcriptomics/striatum

python paper/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper/results/parameters/hippocampus_diff_rf_striatum_global_priors/beta.csv \
  --gamma paper/results/parameters/hippocampus_diff_rf_striatum_global_priors/gamma.csv \
  --out-dir paper/results/transcriptomics/hippocampus
```

By default, this applies the manuscript-style filters `beta > 0` and
posterior-updated `beta`/`gamma` parameters (`ks_pvalue < 0.001`). For robustness
checks like Fig. S2, pass `--no-update-filter` or adjust `--beta-min`.

Optional GSEA requires `gseapy` and internet/cache access to the Enrichr library:

```bash
python paper/analyses/transcriptomics/run_gsea.py \
  --input paper/results/transcriptomics/striatum/gene_eta_correlations.csv \
  --out-dir paper/results/enrichment/striatum
```

Cell-type associations:

```bash
python paper/analyses/cell_types/cell_type_axis_associations.py \
  --axis paper/results/transcriptomics/striatum/region_axis.csv \
  --cell-types paper/data/cell_types/connectome_celltype.csv \
  --out-dir paper/results/cell_types/striatum
```

Compare striatal and hippocampal gene-axis structure:

```bash
python paper/analyses/transcriptomics/compare_axes.py \
  --striatum-dir paper/results/transcriptomics/striatum \
  --hippocampus-dir paper/results/transcriptomics/hippocampus \
  --out-dir paper/results/transcriptomics/striatum_vs_hippocampus
```

## Provenance

This layer ports the relevant analysis logic from:

- `synuclein_spread`: original gene/parameter and PCA exploration
- `gene_enrichment`: pre-ranked GSEA and enrichment plotting
- `cell-type-atlas`: cell-type composition and monoaminergic association analyses

Only the paper-critical pieces should be kept here. General model functionality should
remain in the main package.
