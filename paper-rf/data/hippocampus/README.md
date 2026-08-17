# Hippocampus Paper Inputs

This folder contains the fitting-ready hippocampal pSyn dataset used for the manuscript analyses.

- `observations.csv` contains the hippocampal PFF injection pathology observations.
- `observations_right_ipsi.csv` is the manuscript analysis table, with the right hemisphere mapped to ipsilateral labels.
- `network.csv` is the filtered labeled connectome used for the hippocampus jobs. It is identical to `paper-rf/data/striatum/network.csv` in this repository and was copied here so the hippocampus bundle is self-contained.

The hippocampus jobs use seed indices `[53, 55, 56]`, corresponding to:

- `53`: `iCA1`
- `55`: `iCA3`
- `56`: `iDG`

The current paper configs fit raw replicate observations by default, matching the striatum configs in this repository.
