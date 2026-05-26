# Co-pathology Paper Workflow

This directory contains the project-specific inputs and cluster wrapper for the
synuclein/tau/A-beta co-pathology analyses.

Current `LOCAL-RF` configs:

- `configs/copath_syn_app_inferu0.toml`
- `configs/copath_syn_mapt_inferu0.toml`
- `configs/copath_tau_app_inferu0.toml`
- `configs/copath_tau_mapt_inferu0.toml`
- `configs/copath_syn_app_detu0.toml`
- `configs/copath_syn_mapt_detu0.toml`
- `configs/copath_tau_app_detu0.toml`
- `configs/copath_tau_mapt_detu0.toml`

The `*_inferu0.toml` configs infer one initial condition per region. The
`*_detu0.toml` configs use deterministic regional initial conditions for
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

To submit four inferred-initial-condition chains for each copath dataset on
CUBIC:

```bash
bash paper-copath/run_copath_inferences.sh
```

To submit deterministic-initial-condition chains instead:

```bash
COPATH_U0_MODE=detu0 bash paper-copath/run_copath_inferences.sh
```

By default this uses the `all` partition with a two-day wall time. Override
without editing the script:

```bash
LOCAL_RF_CHAINS=2 SLURM_TIME=1-00:00:00 bash paper-copath/run_copath_inferences.sh
```
