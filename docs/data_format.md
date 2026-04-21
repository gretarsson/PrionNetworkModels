# Data Format

## Network Input

Expected network input:

- square CSV
- first column contains row labels
- remaining columns contain weights
- column headers define canonical region order

## Observation Input

Expected pathology input:

- first column identifies sample or replicate
- second column contains timepoints
- remaining columns contain region-level pathology values
- missing values are allowed

## Region Alignment

Observation data should be aligned to the network region order.

If a region exists in the network but not in the pathology table, the loader may fill that region with missing values if this behavior is explicitly supported by the dataset-loading layer.

## Bilateral Labeling

The current legacy convention uses hemisphere-prefixed labels such as:

- `iCA1`
- `cCA1`

V1 should preserve compatibility with that convention while isolating the pairing logic inside `parameter_sharing.jl`.
