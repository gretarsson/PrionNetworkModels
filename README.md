# PrionNetworkModels

Julia code for fitting network models of prion-like pathology spread.

The package currently supports four model families:

- `DIFF`
- `DIFF-R`
- `DIFF-RF`
- `LOCAL-RF`

The reusable package code is in `src/`. Command-line entry points are in
`scripts/`. Paper-specific analyses are kept in:

- `paper-rf/`: rise-and-fall alpha-synuclein manuscript
- `paper-copath/`: co-pathology analyses

## Installation

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Most paper plotting scripts use Python. The paper workflow expects the virtual
environment at `paper-rf/python/.venv/`; see [paper-rf/README.md](paper-rf/README.md)
for the manuscript commands.

## Run Bundles

Each model fit is saved as one directory under `runs/`:

```text
runs/example-fit/
  spec.toml
  metadata.json
  posterior.h5
  posterior_summary.csv
  diagnostics.json
  predictions_train.csv
  plots/
```

Important files:

- `spec.toml`: model, data paths, seed choices, and inference settings.
- `posterior.h5`: posterior samples.
- `posterior_summary.csv`: parameter summaries.
- `diagnostics.json`: sampler and convergence metadata.
- `predictions_train.csv`: model predictions at observed timepoints.
- `plots/`: generated diagnostic and fit plots.

Large `runs/` folders are not tracked by git. For paper reproduction, download
the archived run bundles separately and place them in `runs/`.

## Input Data Format

Network files are square CSV matrices with region labels:

```csv
region,r1,r2,r3
r1,0,0.2,0.4
r2,0.2,0,0.1
r3,0.4,0.1,0
```

Observation files contain longitudinal pathology measurements:

```csv
sample_id,timepoint,r1,r2,r3
mouse_1,0.1,0.4,0.01,0.00
mouse_1,0.3,0.3,0.02,0.01
mouse_1,1.0,0.2,0.05,0.03
```

## Fitting A Model

Model settings are stored in TOML config files. A minimal example:

```toml
[model]
name = "DIFF-RF"
transport = "retrograde"
parameter_sharing = "independent"

[data]
network = "data/examples/network.csv"
observations = "data/examples/observations.csv"

[seeding]
seed_indices = [1]
infer_seed = true

[inference]
n_chains = 1
sampler = "NUTS"
n_samples = 150
n_warmup = 150
target_acceptance = 0.8

[holdout]
strategy = "none"
```

Run inference:

```bash
julia --project=. scripts/fit_model.jl \
  --config configs/examples/diff_rf.toml \
  --run-id example-fit
```

Plot a completed run:

```bash
julia --project=. scripts/plot_run.jl --run runs/example-fit
```

## Synthetic Example

Generate a small synthetic dataset:

```bash
julia --project=. scripts/reproduce_core_examples.jl
```

Run a short smoke-test fit:

```bash
julia --project=. scripts/smoke_fit_diff_rf.jl
```

For the uncoupled local model:

```bash
julia --project=. scripts/smoke_fit_local_rf.jl
```

## Multi-Chain Runs

Large fits are usually run as one chain per job and merged afterward.

Merge four chains named `runs/striatum_DIFF-RF_RETRO_C1` through
`runs/striatum_DIFF-RF_RETRO_C4`:

```bash
julia --project=. scripts/merge_chains.jl \
  --prefix striatum_DIFF-RF_RETRO \
  --out-run-id striatum_DIFF-RF_RETRO
```

Merge selected chains:

```bash
julia --project=. scripts/merge_chains.jl \
  --prefix hippocampus_DIFF-RF_RETRO \
  --chains 1,3,4 \
  --out-run-id hippocampus_DIFF-RF_RETRO_C1_C3_C4
```

Each merged bundle writes `source_chains.csv`.

## Cluster Helpers

Cluster scripts live in `cluster/` and `scripts/`.

Common commands:

```bash
scripts/run_inferences.sh
scripts/run_hippocampus_inferences.sh
scripts/sync_runs_from_cluster.sh
```

By default, `sync_runs_from_cluster.sh` performs a shallow sync of each top-level
run bundle. Use `--recursive-runs` to include nested subfolders.

## Paper Reproduction

The main manuscript workflow is documented in [paper-rf/README.md](paper-rf/README.md).

Short version:

1. Clone this repository.
2. Download the archived `runs/` artifact.
3. Place the downloaded run folders under `runs/`.
4. Run:

```bash
bash paper-rf/run_main_figure_panels.sh
```

Generated paper outputs are written to `paper-rf/results/` and
`paper-rf/figures/`; both are ignored by git.

## Repository Layout

```text
src/          Julia package code
scripts/      command-line scripts
cluster/      SLURM helpers
configs/      example configs
data/         example data
paper-rf/     alpha-synuclein paper analyses
paper-copath/ co-pathology analyses
runs/         local/archived model outputs, not tracked
docs/         design notes
```

Additional technical details:

- [docs/architecture.md](docs/architecture.md)
- [docs/run_format.md](docs/run_format.md)
