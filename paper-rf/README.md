# Rise-And-Fall Alpha-Synuclein Paper

This folder contains the analyses used to reproduce the paper figures.

Large fitted inference objects are not stored in git. To reproduce the paper
figures from the finalized inferences, download the archived `runs/` artifact
and place the run folders under the repository-level `runs/` directory.

## Required Inputs

Small curated inputs are tracked here:

```text
paper-rf/data/
  striatum/
  hippocampus/
  transcriptomics/
  cell_types/
```

Large posterior bundles should be placed here:

```text
runs/
```

The main paper analyses use these fitted runs:

```text
runs/striatum_DIFF-RF_RETRO/
runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4/
```

Figures 2-5 also use translated model-comparison and null-model bundles under
`runs/`. The null-model folders store WAIC summaries rather than full posterior
chains.

## Reproduce Main Figure Panels

From the repository root:

```bash
bash paper-rf/run_main_figure_panels.sh
```

This writes independent PDF/PNG panels to:

```text
paper-rf/figures/Figure2/
paper-rf/figures/Figure3/
paper-rf/figures/Figure4/
paper-rf/figures/Figure5/
paper-rf/figures/Figure6/
paper-rf/figures/Figure7/
```

Each folder also contains `missing_requirements.md`. It should say `None.` when
all required run bundles are available.

Generated figures and tables are ignored by git:

```text
paper-rf/figures/
paper-rf/results/
```

## Reproduce Biological Analyses

The main biological workflow exports posterior parameters, runs gene-expression
PCA, runs cell-type associations, compares striatum and hippocampus, and creates
the manuscript panels:

```bash
bash paper-rf/run_paper_analyses.sh
```

To include KEGG/GSEA enrichment:

```bash
RUN_GSEA=1 bash paper-rf/run_paper_analyses.sh
```

The workflow uses all transcriptomics-matched regions by default for the final
figure panels.

## Figures 6 And 7

To rebuild only the vulnerability-axis panels:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/rebuild_figures_6_7.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_figure6_ai_panels.py
```

Outputs:

```text
paper-rf/figures/Figure6/
paper-rf/figures/Figure7/
```

## Appendix Vulnerability Panels

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/build_appendix_vulnerability_inputs.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_appendix_vulnerability_figures.py
```

## Optional: Import Old Inference Files

The final reproduction workflow does not require the old `synuclein_spread`
repository if the archived `runs/` artifact has already been downloaded.

For provenance, this repository includes a one-time importer that converts old
`synuclein_spread` `.jls` inference files into the current run-bundle format:

```bash
julia --project=/path/to/synuclein_spread \
  paper-rf/analyses/figures/import_synuclein_spread_runs.jl
julia --project=. paper-rf/analyses/figures/export_run_bundle_waic.jl
```

After import, the figure scripts read from `runs/` and no longer need
`synuclein_spread`.

## Useful Individual Commands

Export parameter tables:

```bash
julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/striatum_DIFF-RF_RETRO \
  --out-dir paper-rf/results/parameters/striatum_diff_rf

julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4 \
  --out-dir paper-rf/results/parameters/hippocampus_diff_rf
```

Run gene-parameter PCA:

```bash
python paper-rf/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper-rf/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper-rf/results/parameters/striatum_diff_rf/beta.csv \
  --gamma paper-rf/results/parameters/striatum_diff_rf/gamma.csv \
  --out-dir paper-rf/results/transcriptomics/striatum
```

Run cell-type associations:

```bash
python paper-rf/analyses/cell_types/cell_type_axis_associations.py \
  --axis paper-rf/results/transcriptomics/striatum/region_axis.csv \
  --cell-types paper-rf/data/cell_types/connectome_celltype.csv \
  --out-dir paper-rf/results/cell_types/striatum
```
