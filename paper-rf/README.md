# Rise-and-fall dynamics reveal a molecular and cellular vulnerability axis in prion-like α-synuclein propagation

This folder contains the scripts used to reproduce the manuscript analyses and
figure panels. The preprint is available on bioRxiv:

```text
https://doi.org/10.64898/2026.03.27.714785
```

Inference is run in Julia. Downstream analyses and plotting are run in Python.

The pathology, network, transcriptomic, and cell-type data are included in
`paper-rf/data/`. The posterior inference chains used in the paper are too large
for git and are archived on Zenodo:

```text
https://doi.org/10.5281/zenodo.21045203
```

After downloading the inference archive, place the run folders in the
repository-level `runs/` folder.

## 1. Run The Inferences

This stage fits the dynamical model parameters to the pathology data with MCMC
and merges the resulting posterior chains.

The paper can be reproduced directly from the archived inference chains. To
rerun the main inferences instead, use:

```bash
bash scripts/run_inferences.sh
bash scripts/run_hippocampus_inferences.sh
```

These scripts submit the model fits defined in `paper-rf/configs/`. Individual
fits can also be run with:

```bash
julia --project=. scripts/fit_model.jl \
  --config paper-rf/configs/striatum_diff_rf_core.toml \
  --run-id striatum_DIFF-RF_RETRO_C1 \
  --progress
```

After fitting separate chains, merge them with:

```bash
julia --project=. scripts/merge_chains.jl \
  --prefix striatum_DIFF-RF_RETRO \
  --chains 1,2,3,4 \
  --out-run-id striatum_DIFF-RF_RETRO

julia --project=. scripts/merge_chains.jl \
  --prefix hippocampus_DIFF-RF_RETRO_striatum-global-priors \
  --chains 1,2,3,4 \
  --out-run-id hippocampus_DIFF-RF_RETRO_striatum-global-priors
```

For any new run, inspect the chain diagnostics before using the merged
posterior. In the manuscript analysis, lower-likelihood hippocampal chains were
excluded after this diagnostic check; the exact posterior bundles used in the
paper are provided in the Zenodo inference archive.

## 2. Run The Analyses

This stage exports posterior parameter estimates and relates them to regional
gene expression and cell-type composition.

Create the Python environment:

```bash
python3 -m venv paper-rf/python/.venv
paper-rf/python/.venv/bin/python -m pip install -r paper-rf/python/requirements.txt
```

Run the manuscript analyses:

```bash
RUN_GSEA=1 bash paper-rf/run_paper_analyses.sh
```

If using different merged run names, set `STRIATUM_RUN` and `HIPPO_RUN` when
running this command.

The analysis outputs are written to `paper-rf/results/`.

## 3. Make The Figure Panels

This stage converts the model-comparison, model-fit, transcriptomic, and
cell-type results into the figure panels used in the manuscript.

Create the independent figure panels for Figures 2-7:

```bash
bash paper-rf/run_main_figure_panels.sh
```

The panels are written to:

```text
paper-rf/figures/Figure2/
paper-rf/figures/Figure3/
paper-rf/figures/Figure4/
paper-rf/figures/Figure5/
paper-rf/figures/Figure6/
paper-rf/figures/Figure7/
```

Appendix panels can be regenerated with:

```bash
paper-rf/python/.venv/bin/python paper-rf/analyses/build_appendix_vulnerability_inputs.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/update_appendix_vulnerability_figures.py
paper-rf/python/.venv/bin/python paper-rf/analyses/plotting/make_hippocampus_appendix_diagnostics.py
paper-rf/python/.venv/bin/python paper-rf/analyses/transcriptomics/tca_oxphos_machinery_eta.py
```

Generated outputs are written to `paper-rf/results/` and `paper-rf/figures/`.
