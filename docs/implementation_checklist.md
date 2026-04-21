# Implementation Checklist

## Phase 1: Skeleton and Contracts

- create new standalone repository
- define folder layout
- define run bundle contract
- define v1 model names and public terminology
- define configuration schema

## Phase 2: Core Port

- port `DIFF`
- port `DIFF-R`
- port `DIFF-RF`
- port transport mechanisms:
  - retrograde
  - anterograde
  - bidirectional
  - euclidean
- port bilateral parameter-sharing logic
- port dataset loading and region alignment

## Phase 3: Inference Pipeline

- implement config parsing
- implement resolved run spec generation
- implement ODE problem builder
- implement MCMC runner
- implement per-chain output handling
- implement chain merge workflow
- implement HDF5 posterior writing

## Phase 4: Core Analysis

- posterior summary extraction
- predicted vs observed
- retrodiction by region
- run diagnostics
- `beta` vs `gamma` plot for `DIFF-RF`
- out-of-sample evaluation

## Phase 5: Dataset Packaging

- identify processed striatal fitting-ready dataset
- identify processed hippocampal fitting-ready dataset
- copy only required files into `data/paper/`
- document provenance of each kept file

## Phase 6: Validation

- choose one reference striatal workflow
- confirm new repo reproduces core plots
- choose one reference hippocampal workflow
- compare outputs qualitatively and quantitatively
- document known differences from legacy repo

## Phase 7: Future Extensions

- transcriptomic correlation workflow
- vulnerability-axis analysis
- GSEA pipeline hooks
- additional dynamical processes
- package-quality tests and CI
