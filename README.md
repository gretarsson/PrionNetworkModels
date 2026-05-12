# PrionNetworkModels

`PrionNetworkModels` is a Julia toolkit for fitting and analyzing network models of prion-like neurodegenerative spread.

The current public model names are:

- `DIFF`
- `DIFF-R`
- `DIFF-RF`

This repository is intentionally separate from the legacy `synuclein_spread` archive so the new workflow can be cleaner, slimmer, and easier for others to use.

## What The Workflow Looks Like

The intended workflow is:

1. prepare a network CSV and an observations CSV,
2. write a config file that points to those files and chooses a model,
3. run inference,
4. save a structured run bundle,
5. make plots directly from that run bundle.

The important change compared with the old repo is that we do not save a Julia-specific opaque blob and then try to remember later what it meant. Instead, each run has its own folder with a stable layout.

## Run Bundles

Every model fit writes to a folder under `runs/`, for example:

```text
runs/smoke-fit-diff-rf/
  spec.toml
  metadata.json
  posterior.h5
  posterior_summary.csv
  diagnostics.json
  predictions_train.csv
  plots/
```

What each file means:

- `spec.toml`
  - the run recipe you asked for
  - model name, transport mode, data paths, seed choices, and inference settings
- `metadata.json`
  - bundle bookkeeping
  - run id, creation metadata, and a manifest of the files in the bundle
- `posterior.h5`
  - the full posterior output from MCMC
  - this is the main raw inference result
- `posterior_summary.csv`
  - quick human-readable summaries of the posterior parameters
  - useful for inspection without opening the full posterior file
- `diagnostics.json`
  - basic run information such as sampler and iteration counts
  - later this should also hold richer convergence diagnostics
- `predictions_train.csv`
  - model predictions at the observation timepoints
  - useful for predicted-vs-observed comparisons
- `plots/`
  - generated figures for that run
  - for example predicted-vs-observed and one retrodiction plot per region

The main idea is that one fit lives in one place. You should not have to remember which script produced which figure or which posterior belongs to which dataset.

The observed dataset itself is not copied into every run bundle. The run spec points back to the canonical dataset under `data/`, which avoids duplicating the same observation matrix across many runs.

To keep runs portable between cluster and local environments, the run-bundle `spec.toml` stores repo-relative data paths like `data/examples/observations.csv` whenever possible instead of machine-specific absolute paths.

## Input Files

### 1. Network CSV

The network file is a square matrix with region labels.

Expected format:

- first column: region labels
- remaining columns: edge weights
- column headers: region labels in the same order

Example:

```csv
region,r1,r2,r3
r1,0,0.2,0.4
r2,0.2,0,0.1
r3,0.4,0.1,0
```

### 2. Observation CSV

The observation file is longitudinal pathology data.

Expected format:

- first column: sample or replicate ID
- second column: timepoint
- remaining columns: pathology values per region

Example:

```csv
sample_id,timepoint,r1,r2,r3
mouse_1,0.1,0.4,0.01,0.00
mouse_1,0.3,0.3,0.02,0.01
mouse_1,1.0,0.2,0.05,0.03
```

## How To Choose The Model And Run Settings

You do that in a TOML config file under `configs/`.

Example: [diff_r.toml](/Users/gretarsson/Desktop/PrionNetworkModels/configs/examples/diff_r.toml)

```toml
[model]
name = "DIFF-R"
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
target_acceptance = 0.8
sampler = "NUTS"
n_samples = 150
n_warmup = 150

[holdout]
strategy = "none"
```

### What each section means

- `[model]`
  - `name`: which dynamical model to fit
  - `transport`: `retrograde`, `anterograde`, `bidirectional`, or `euclidean`
  - `parameter_sharing`: currently `independent`, with bilateral sharing planned as a first-class mode

- `[data]`
  - `network`: path to the connectivity matrix
  - `observations`: path to the pathology data

- `[seeding]`
  - `seed_indices`: which region indices are the initial seed sites
  - `infer_seed`: whether the initial seed magnitude should be inferred

- `[inference]`
  - `sampler`: currently `NUTS` or `MH`
  - `n_samples`: number of posterior samples to keep
  - `n_warmup`: warmup/adaptation iterations for `NUTS`
  - `target_acceptance`: `NUTS` target acceptance rate

- `[holdout]`
  - reserved for out-of-sample workflows

## How Priors Work Right Now

Right now, priors are defined in code in [inference.jl](/Users/gretarsson/Desktop/PrionNetworkModels/src/inference.jl), in the `default_priors` function.

That means:

- model choice, transport, seeds, and inference settings are user-facing in the config,
- priors are still code-facing for now.

This is an honest temporary state while the core workflow is being stabilized. The next step will be to make priors configurable by profile, and then configurable more directly from TOML.

## Quick Start: Synthetic Example

### Step 1. Generate a tiny synthetic dataset

This makes a random `N=10` network and simulated DIFF-RF observations:

```bash
julia --project=/Users/gretarsson/Desktop/PrionNetworkModels \
  /Users/gretarsson/Desktop/PrionNetworkModels/scripts/reproduce_core_examples.jl
```

This writes:

- [network.csv](/Users/gretarsson/Desktop/PrionNetworkModels/data/examples/network.csv)
- [observations.csv](/Users/gretarsson/Desktop/PrionNetworkModels/data/examples/observations.csv)
- [observations_summary.csv](/Users/gretarsson/Desktop/PrionNetworkModels/data/examples/observations_summary.csv)
- [generating_parameters_diff_rf.csv](/Users/gretarsson/Desktop/PrionNetworkModels/data/examples/generating_parameters_diff_rf.csv)

### Step 2. Run Bayesian inference on that synthetic example

Option A: use the dedicated smoke-test script

```bash
julia --project=/Users/gretarsson/Desktop/PrionNetworkModels \
  /Users/gretarsson/Desktop/PrionNetworkModels/scripts/smoke_fit_diff_rf.jl
```

Option B: use the generic runner

```bash
julia --project=/Users/gretarsson/Desktop/PrionNetworkModels \
  /Users/gretarsson/Desktop/PrionNetworkModels/scripts/fit_model.jl \
  --config /Users/gretarsson/Desktop/PrionNetworkModels/configs/examples/diff_r.toml \
  --run-id my-first-fit \
  --samples 150 \
  --warmup 150
```

### Step 3. Make plots from the run bundle

```bash
julia --project=/Users/gretarsson/Desktop/PrionNetworkModels \
  /Users/gretarsson/Desktop/PrionNetworkModels/scripts/plot_run.jl \
  --run /Users/gretarsson/Desktop/PrionNetworkModels/runs/my-first-fit
```

This writes plots under:

```text
runs/my-first-fit/plots/
  predicted_vs_observed.pdf
  retrodiction/
```

## Cluster Workflow

For larger real datasets, the intended workflow is:

1. submit one chain per cluster job,
2. merge the finished chain runs into one combined run,
3. plot from the merged run.

The paper-oriented striatum retrograde submission script is:

```bash
/Users/gretarsson/Desktop/PrionNetworkModels/scripts/run_inferences.sh
```

That submits four single-chain jobs for each of:

- `DIFF`
- `DIFF-R`
- `DIFF-RF`

To sync all finished run folders back from the cluster to your local machine:

```bash
/Users/gretarsson/Desktop/PrionNetworkModels/scripts/sync_runs_from_cluster.sh
```

By default, that script pulls from:

- `alexanderc@cubic-login1:~/PrionNetworkModels/runs/`

and syncs into:

- `/Users/gretarsson/Desktop/PrionNetworkModels/runs/`

You can override the remote host or project path without editing the script, for example:

```bash
REMOTE_HOST=alexanderc@cubic-login5 \
REMOTE_PROJECT_DIR=~/PrionNetworkModels \
/Users/gretarsson/Desktop/PrionNetworkModels/scripts/sync_runs_from_cluster.sh
```

If the cluster jobs create runs like:

- `runs/striatum_DIFF-RF_RETRO_C1`
- `runs/striatum_DIFF-RF_RETRO_C2`
- `runs/striatum_DIFF-RF_RETRO_C3`
- `runs/striatum_DIFF-RF_RETRO_C4`

then you merge them with:

```bash
julia --project=/Users/gretarsson/Desktop/PrionNetworkModels \
  /Users/gretarsson/Desktop/PrionNetworkModels/scripts/merge_chains.jl \
  --prefix striatum_DIFF-RF_RETRO \
  --out-run-id striatum_DIFF-RF_RETRO
```

and then plot the merged run with:

```bash
julia --project=/Users/gretarsson/Desktop/PrionNetworkModels \
  /Users/gretarsson/Desktop/PrionNetworkModels/scripts/plot_run.jl \
  --run /Users/gretarsson/Desktop/PrionNetworkModels/runs/striatum_DIFF-RF_RETRO
```

## What Plots Are Implemented Now

Current plotting support:

- predicted vs observed
- retrodiction by region

These plots are generated from the run bundle. For retrodiction, the plotting code also reads `posterior.h5`, simulates a dense trajectory from the posterior mean parameters, and uses the inferred observation noise to show:

- the posterior mean trajectory,
- 50% and 90% noise bands,
- observed mean data points,
- a shared y-axis across all regional retrodiction plots in the same run.

The plotting code currently reads the canonical observations file from `spec.toml` and computes replicate summaries on the fly. That preserves the raw replicates in one place and leaves room to add future plot modes such as:

- plot the mean only
- plot all replicate observations
- plot mean with standard-error bars

Planned next:

- diagnostics plots from posterior chains
- out-of-sample plots
- `beta` vs `gamma` plots for `DIFF-RF`

## Repository Layout

- `src/`
  - core code
- `scripts/`
  - command-line entrypoints
- `configs/`
  - example and paper configs
- `data/`
  - example and curated paper inputs
- `runs/`
  - fit outputs
- `docs/`
  - design docs

For more design detail, see:

- [architecture.md](/Users/gretarsson/Desktop/PrionNetworkModels/docs/architecture.md)
- [run_format.md](/Users/gretarsson/Desktop/PrionNetworkModels/docs/run_format.md)
- [implementation_checklist.md](/Users/gretarsson/Desktop/PrionNetworkModels/docs/implementation_checklist.md)
