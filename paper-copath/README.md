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

## Independent region-wise RF fits

For the exploratory independent regional model, submit all four co-pathology
datasets with one script:

```bash
bash paper-copath/run_region_rf_copath.sh
```

This submits one SLURM array for each dataset: `syn_app`, `syn_mapt`, `tau_app`,
and `tau_mapt`. Each array task fits one brain region with its own `alpha`,
`beta`, `gamma`, `u0`, and `sigma`. The default array is `1-412%40`, so CUBIC can
run up to 40 regions at once while each region still runs four chains. The
default priors used by this exploratory script are:

```text
alpha ~ Normal+(0, 1.0)
beta  ~ Normal(0, 1)
gamma ~ Normal+(0, 0.1)
u0    ~ Normal+(0, 0.01)
sigma ~ LogNormal(0, 1)
```

The default ODE `maxiters` for these array jobs is `50000`. The script uses all
replicate observations by default; set `REGION_RF_MEAN_DATA=1` only if you want
to fit to timepoint means instead.

Outputs are stored under:

```text
runs/region_rf/copath_<dataset>/
```

Each region gets its own directory containing `posterior.h5`,
`posterior_summary.csv`, `diagnostics.csv`, `predictions_train.csv`, and fit/trace
plots. After the arrays finish, collect all four datasets into tidy summary
tables with:

```bash
bash paper-copath/collect_region_rf_copath.sh
```

This writes `region_rf_summary.csv` and
`region_rf_posterior_summary_long.csv` in each dataset output directory.

To build the assembled REGION-RF plots for all four datasets, run:

```bash
bash paper-copath/plot_region_rf_copath.sh
```

This writes `predictions_train.csv`, `plots/predicted_vs_observed.*`,
`plots/diagnostics/`, and `plots/retrodiction/` into each
`runs/region_rf/copath_<dataset>/` folder.
