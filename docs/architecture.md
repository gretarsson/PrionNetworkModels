# Architecture

## Purpose

`PrionNetworkModels` is intended to be a reusable modeling toolkit for prion-like propagation on networks, with a clean separation between:

- model definitions,
- inference configuration,
- run storage,
- posterior simulation,
- plotting and diagnostics,
- later biological downstream analyses.

## V1 Scope

The first version focuses on:

- `DIFF`
- `DIFF-R`
- `DIFF-RF`

with support for:

- retrograde transport,
- anterograde transport,
- bidirectional transport,
- euclidean transport,
- bilateral parameter sharing,
- multi-seed inference,
- structured holdout settings.

## Core Design Principle

The central design choice is to replace ad hoc serialized outputs with a structured run bundle.

Every inference run should be reproducible from a resolved run specification and should save enough metadata to:

- regenerate posterior simulations,
- identify training vs held-out observations,
- build standard plots without guessing file relationships,
- later attach downstream analyses such as transcriptomic correlations.

## Module Responsibilities

- `models.jl`
  - ODE definitions for `DIFF`, `DIFF-R`, and `DIFF-RF`
- `priors.jl`
  - prior profiles and prior construction
- `data_io.jl`
  - pathology dataset loading and validation
- `networks.jl`
  - connectivity loading and transport-mode transforms
- `parameter_sharing.jl`
  - independent vs bilateral pairing logic
- `problem_builder.jl`
  - ODEProblem construction from model + data + configuration
- `inference.jl`
  - MCMC orchestration
- `run_bundle.jl`
  - run bundle read/write logic
- `posterior.jl`
  - posterior summaries and parameter extraction
- `diagnostics.jl`
  - `Rhat`, ESS, log-likelihood, and convergence summaries
- `prediction.jl`
  - posterior simulation, retrodiction, held-out prediction
- `plotting.jl`
  - standard figure generation

## Planned Run Bundle

Each run will live in its own folder under `runs/`.

Suggested contents:

```text
runs/<run_id>/
  spec.toml
  metadata.json
  posterior.h5
  posterior_summary.csv
  diagnostics.json
  predictions_train.csv
  predictions_full.csv
```

The observed datasets themselves should usually stay in `data/` and be referenced from `spec.toml` rather than copied into each run directory.

## Storage Strategy

- Configuration: `TOML`
- Metadata and structured summaries: `JSON`
- Tabular summaries: `CSV`
- Posterior draws and dense numeric arrays: `HDF5`

This keeps metadata inspectable while avoiding the brittleness of Julia-specific serialized blobs.
