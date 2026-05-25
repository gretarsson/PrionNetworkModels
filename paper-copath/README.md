# Co-pathology Paper Workflow

This directory contains the project-specific inputs and cluster wrapper for the
synuclein/tau/A-beta co-pathology analyses.

Current `LOCAL-RF` configs:

- `configs/syn_app_local_rf.toml`
- `configs/syn_app_local_rf_fixed_u0.toml`
- `configs/syn_mapt_local_rf_fixed_u0.toml`
- `configs/tau_app_local_rf_fixed_u0.toml`
- `configs/tau_mapt_local_rf_fixed_u0.toml`
- `configs/syn_mapt_local_rf.toml`
- `configs/tau_app_local_rf.toml`
- `configs/tau_mapt_local_rf.toml`

The `*_fixed_u0.toml` config uses deterministic regional initial conditions for
`LOCAL-RF`:

```toml
[seeding]
infer_local_u0 = false
```

When `infer_local_u0 = false`, `LOCAL-RF` sets every regional initial condition
to `local_u0_value`. The default value is `3.364e-5`, the all-region mean
pathology at the first striatum timepoint used for the exploratory scale check.
Override it with:

```toml
local_u0_value = 1.0e-4
```

Omit `infer_local_u0` or set it to `true` to infer regional initial conditions.

To submit four chains for each dataset on CUBIC:

```bash
bash paper-copath/run_local_rf_inferences.sh
```

To submit four chains for each dataset with deterministic local initial
conditions:

```bash
bash paper-copath/run_copath_inferences.sh
```

By default this uses the `all` partition with a two-day wall time. Override
without editing the script:

```bash
LOCAL_RF_CHAINS=2 SLURM_TIME=1-00:00:00 bash paper-copath/run_local_rf_inferences.sh
```
