# Run Format

## Overview

Each inference run should produce a self-contained run bundle under `runs/<run_id>/`.

The run bundle is the contract between:

- inference,
- posterior summarization,
- diagnostics,
- plotting,
- downstream analysis.

## Required Files

- `spec.toml`
  - user-authored or resolved run specification
- `metadata.json`
  - fully resolved metadata for the run
- `posterior.h5`
  - posterior samples and numeric arrays
- `posterior_summary.csv`
  - summary statistics for each inferred parameter
- `diagnostics.json`
  - convergence and fit diagnostics

## Recommended Additional Files

- `predictions_train.csv`
- `predictions_full.csv`

Observed datasets should usually remain canonical under `data/` and be referenced from `spec.toml` rather than copied into every run bundle. This avoids duplicating the same observation matrices across many runs.

When the dataset lives inside the repository, those references should be stored as repo-relative paths so the same run bundle can move cleanly between environments that preserve the repo layout.

## Metadata Requirements

The metadata layer should include:

- model name
- transport mechanism
- parameter-sharing mode
- network source
- observation source
- region labels
- seed regions
- whether seed was inferred
- training timepoints
- held-out timepoints
- prior profile
- sampler settings

## Posterior HDF5 Layout

Initial proposal:

- `/chains/samples`
  - posterior draw matrix
- `/chains/parameter_names`
  - ordered parameter names
- `/chains/chain_ids`
  - chain indices
- `/model/state_names`
  - state vector labels
- `/data/region_labels`
  - region labels
- `/data/timepoints_train`
  - fitted timepoints
- `/data/timepoints_full`
  - full timepoints when available

This layout may change once the first concrete implementation is in place, but the goal is to preserve stable semantic meaning rather than opaque object dumps.
