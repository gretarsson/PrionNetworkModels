# Hippocampus Paper Inputs

This folder contains the fitting-ready hippocampal pSyn dataset ported from the legacy `synuclein_spread` workflow.

- `observations.csv` comes from `synuclein_spread/data/hippocampal/hippocampal_syn_only.csv`.
- `network.csv` is the filtered labeled connectome used by the legacy hippocampus jobs. It is identical to `data/paper/striatum/network.csv` in this repository and was copied here so the hippocampus bundle is self-contained.

The legacy hippocampus jobs used seed indices `[53, 55, 56]`, corresponding to:

- `53`: `iCA1`
- `55`: `iCA3`
- `56`: `iDG`

The current paper configs fit raw replicate observations by default, matching the striatum configs in this repository.
