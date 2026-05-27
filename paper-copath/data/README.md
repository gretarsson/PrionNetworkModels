# Synuclein/Tau/A-beta Co-pathology Data

This directory contains co-pathology paper inputs copied from the legacy
`synuclein_spread/data/syn_tau_abeta` archive.

The `LOCAL-RF` configs currently use the four longitudinal tau/synuclein
pathology tables:

- `syn_pathology_app.csv`
- `syn_pathology_mapt.csv`
- `tau_pathology_app.csv`
- `tau_pathology_mapt.csv`

Each table uses the standard observation format:

- column 1: mouse/sample id
- column 2: months post injection (`mpi`)
- remaining columns: regional pathology values

The local model is fit on the 412-region connectome subset using
`network.csv` as the canonical region index. `LOCAL-RF` does not use network
coupling, but the current run machinery still uses a network/index file to
define region order and to subset the wider 564-column pathology tables.

These inputs are intentionally kept under `paper-copath/data/`, separate from
the rise-and-fall paper inputs in `paper-rf/data/`.

The A-beta inputs are regional wide tables copied from the same archive:

- `ab40_pathology_mpff.csv`
- `ab42_pathology_mpff.csv`
- `ab40_pathology_adphf.csv`
- `ab42_pathology_adphf.csv`
- `ab40_pathology_control.csv`
- `ab42_pathology_control.csv`

The source A-beta tables contain APP/MAPTApp KI injected treatments (`mPFF`,
`AD PHF`) and non-injected controls (`none`). They do not contain a separate
MAPT amyloid condition, so downstream APP minus MAPT parameter shifts are
compared against treatment minus control A-beta burden.

The deterministic `LOCAL-RF` initial-condition test uses `u0 = 3.364e-5`,
the all-region mean pathology at the first striatum timepoint from
`paper-rf/data/striatum/observations.csv`.
