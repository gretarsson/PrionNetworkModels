# Rise-And-Fall Alpha-Synuclein Paper Analyses

This folder contains the code and small input tables used to reproduce the
scientific panels in the manuscript.

Large posterior inference files are not stored in git. To reproduce the paper
figures exactly, download the archived inference bundles and place the run
folders in the repository-level `runs/` directory.

## Setup

Run commands from the repository root:

```bash
cd /path/to/PrionNetworkModels
```

Install the Julia package environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Create the Python environment used by the paper scripts:

```bash
python3 -m venv paper-rf/python/.venv
paper-rf/python/.venv/bin/python -m pip install -r paper-rf/python/requirements.txt
```

## Inputs

Small inputs tracked in git are under:

```text
paper-rf/data/
  striatum/
  hippocampus/
  transcriptomics/
  cell_types/
```

The final hippocampal analyses use:

```text
paper-rf/data/hippocampus/observations_right_ipsi.csv
```

This table maps the right hemisphere to ipsilateral labels, matching the
hippocampal injection side.

## Inference Bundles

The main manuscript analyses use these fitted run bundles:

```text
runs/striatum_DIFF-RF_RETRO/
runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4/
```

The hippocampal folder name is historical. Its `source_chains.csv` records the
actual retained chains; for the final right-ipsilateral hippocampal inference,
chains 3 and 4 are retained.

Figures 2-5 also use additional model-comparison, null-model, and held-out
evaluation bundles under `runs/`. The null-model folders contain WAIC summaries
rather than full posterior chains.

## Fast Route: Use Archived Inferences

For reproducing the submitted figures, use the archived inference bundles rather
than rerunning MCMC:

1. Download the Zenodo archive of posterior inference chains.
2. Unpack it so that the run folders sit directly under `runs/`.
3. Check that the main folders exist:

```bash
test -f runs/striatum_DIFF-RF_RETRO/posterior.h5
test -f runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4/posterior.h5
```

Then continue with the analysis and plotting commands below.

## Full Route: Refit The Main Inferences

Refitting is computationally expensive. Each command below writes one chain to a
separate folder in `runs/`.

Fit four striatal DIFF-RF retrograde chains:

```bash
for chain in 1 2 3 4; do
  julia --project=. scripts/fit_model.jl \
    --config paper-rf/configs/striatum_diff_rf_core.toml \
    --run-id striatum_DIFF-RF_RETRO_C${chain} \
    --progress
done
```

Merge the striatal chains:

```bash
julia --project=. scripts/merge_chains.jl \
  --prefix striatum_DIFF-RF_RETRO \
  --chains 1,2,3,4 \
  --out-run-id striatum_DIFF-RF_RETRO
```

The hippocampal DIFF-RF run uses posterior priors from the merged striatal
DIFF-RF run, so `runs/striatum_DIFF-RF_RETRO/posterior.h5` must exist first.

Fit four hippocampal DIFF-RF retrograde chains:

```bash
for chain in 1 2 3 4; do
  julia --project=. scripts/fit_model.jl \
    --config paper-rf/configs/hippocampus_diff_rf_striatum_global_priors.toml \
    --run-id hippocampus_DIFF-RF_RETRO_striatum-global-priors_C${chain} \
    --progress
done
```

Merge all hippocampal chains for diagnostics:

```bash
julia --project=. scripts/merge_chains.jl \
  --prefix hippocampus_DIFF-RF_RETRO_striatum-global-priors \
  --chains 1,2,3,4 \
  --out-run-id hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C2_C3_C4
```

Merge the retained high-likelihood hippocampal chains for the main biological
analyses:

```bash
julia --project=. scripts/merge_chains.jl \
  --prefix hippocampus_DIFF-RF_RETRO_striatum-global-priors \
  --chains 3,4 \
  --out-run-id hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4
```

On a SLURM cluster, use the helper scripts in `cluster/` and `scripts/` instead
of running the loop commands interactively. For example:

```bash
bash scripts/run_inferences.sh
bash scripts/run_hippocampus_inferences.sh
```

## Run Analyses

The main biological workflow exports posterior parameters, computes the
vulnerability axes, runs cell-type associations, compares striatal and
hippocampal datasets, and writes intermediate tables:

```bash
bash paper-rf/run_paper_analyses.sh
```

To also run KEGG/GSEA enrichment:

```bash
RUN_GSEA=1 bash paper-rf/run_paper_analyses.sh
```

The final manuscript panels use all transcriptomics-matched regions. Older
filtering levels are available for sensitivity checks but are not needed for the
main figures.

## Make Figure Panels

Create independent PDF/PNG panels for Figures 2-7:

```bash
bash paper-rf/run_main_figure_panels.sh
```

Outputs are written to:

```text
paper-rf/figures/Figure2/
paper-rf/figures/Figure3/
paper-rf/figures/Figure4/
paper-rf/figures/Figure5/
paper-rf/figures/Figure6/
paper-rf/figures/Figure7/
```

Each folder may contain `missing_requirements.md`. It should say `None.` when
all required run bundles are present.

To rebuild only the Figure 6 and Figure 7 vulnerability-axis panels:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/rebuild_figures_6_7.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_figure6_ai_panels.py
```

## Appendix Figures

Rebuild the vulnerability-axis appendix panels:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/build_appendix_vulnerability_inputs.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_appendix_vulnerability_figures.py
```

Rebuild the hippocampal diagnostics appendix panels:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/make_hippocampus_appendix_diagnostics.py
```

Rebuild the curated TCA/OxPhos energy-metabolism panels:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/transcriptomics/tca_oxphos_machinery_eta.py
```

## Useful Individual Commands

Export posterior parameter tables:

```bash
julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/striatum_DIFF-RF_RETRO \
  --out-dir paper-rf/results/parameters/striatum_diff_rf

julia --project=. paper-rf/analyses/model_parameters/export_parameter_tables.jl \
  --run runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4 \
  --out-dir paper-rf/results/parameters/hippocampus_diff_rf
```

Run gene-parameter PCA for one dataset:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/transcriptomics/gene_parameter_pca.py \
  --expression paper-rf/data/transcriptomics/avg_Pangea_exp.csv \
  --beta paper-rf/results/parameters/striatum_diff_rf/beta.csv \
  --gamma paper-rf/results/parameters/striatum_diff_rf/gamma.csv \
  --out-dir paper-rf/results/transcriptomics/striatum
```

Run cell-type associations for one dataset:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/cell_types/cell_type_axis_associations.py \
  --axis paper-rf/results/transcriptomics/striatum/region_axis.csv \
  --cell-types paper-rf/data/cell_types/connectome_celltype.csv \
  --out-dir paper-rf/results/cell_types/striatum
```

## Generated Outputs

Generated figures and result tables are ignored by git:

```text
paper-rf/figures/
paper-rf/results/
runs/
```

These folders can be regenerated from the archived inputs and the scripts above.
