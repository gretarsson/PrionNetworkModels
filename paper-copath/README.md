# Co-pathology Paper Workflow

This directory contains the project-specific inputs and cluster wrapper for the
synuclein/tau/A-beta co-pathology analyses.

Current `LOCAL-RF` configs:

- `configs/syn_app_local_rf.toml`
- `configs/syn_mapt_local_rf.toml`
- `configs/tau_app_local_rf.toml`
- `configs/tau_mapt_local_rf.toml`

To submit four chains for each dataset on CUBIC:

```bash
bash paper-copath/run_local_rf_inferences.sh
```

By default this uses the `all` partition with a two-day wall time. Override
without editing the script:

```bash
LOCAL_RF_CHAINS=2 SLURM_TIME=1-00:00:00 bash paper-copath/run_local_rf_inferences.sh
```
