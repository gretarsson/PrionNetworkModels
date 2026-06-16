#!/usr/bin/env python3
"""Recompute REGION-RF summaries after dropping the worst chain in non-converged fits.

This is a diagnostic sensitivity analysis. It does not modify the original
REGION-RF run directories.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import norm


DATASETS = ["syn_app", "syn_mapt", "tau_app", "tau_mapt"]
OBSERVATION_FILES = {
    "syn_app": "syn_pathology_app.csv",
    "syn_mapt": "syn_pathology_mapt.csv",
    "tau_app": "tau_pathology_app.csv",
    "tau_mapt": "tau_pathology_mapt.csv",
}
MAIN_PARAMETERS = ["alpha", "beta", "gamma"]
ALL_PARAMETERS = ["alpha", "beta", "gamma", "u0", "sigma"]
RHAT_CUTOFF = 1.05


def decode_names(values) -> list[str]:
    names = []
    for value in values:
        if isinstance(value, bytes):
            names.append(value.decode())
        else:
            names.append(str(value))
    return names


def load_h5(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    with h5py.File(path, "r") as h5:
        samples = np.asarray(h5["chains/samples"])
        chain_ids = np.asarray(h5["chains/chain_ids"], dtype=int)
        parameter_names = decode_names(h5["chains/parameter_names"][()])
    if samples.shape[0] == len(parameter_names):
        samples = samples.T
    if samples.shape[0] != len(chain_ids):
        raise ValueError(f"Sample/chain length mismatch in {path}: {samples.shape}, {chain_ids.shape}")
    return samples.astype(float), chain_ids, parameter_names


def split_rhat(values: np.ndarray, chain_ids: np.ndarray) -> float:
    chains = []
    for chain_id in sorted(np.unique(chain_ids)):
        chain_values = values[chain_ids == chain_id]
        chain_values = chain_values[np.isfinite(chain_values)]
        half = len(chain_values) // 2
        if half < 2:
            continue
        chains.append(chain_values[:half])
        chains.append(chain_values[-half:])
    if len(chains) < 2:
        return np.nan
    n = min(len(c) for c in chains)
    arr = np.vstack([c[:n] for c in chains])
    within = np.var(arr, axis=1, ddof=1)
    w = float(np.mean(within))
    b = float(n * np.var(np.mean(arr, axis=1), ddof=1))
    if w == 0:
        return 1.0 if b == 0 else np.inf
    var_plus = ((n - 1) / n) * w + b / n
    return float(np.sqrt(var_plus / w))


def summarize_samples(samples: np.ndarray, chain_ids: np.ndarray, parameter_names: list[str]) -> pd.DataFrame:
    rows = []
    for idx, parameter in enumerate(parameter_names):
        values = samples[:, idx]
        rows.append(
            {
                "parameter": parameter,
                "mean": float(np.nanmean(values)),
                "std": float(np.nanstd(values, ddof=1)),
                "mcse": float(np.nanstd(values, ddof=1) / np.sqrt(np.isfinite(values).sum())),
                "ess_bulk": np.nan,
                "ess_tail": np.nan,
                "rhat": split_rhat(values, chain_ids),
                "ess_per_sec": np.nan,
            }
        )
    return pd.DataFrame(rows)


def region_solution(timepoints: np.ndarray, alpha: float, beta: float, gamma: float, u0: float) -> np.ndarray | None:
    if not all(np.isfinite([alpha, beta, gamma, u0])):
        return None
    if max(timepoints) <= 0:
        return np.full_like(timepoints, u0, dtype=float)

    def rhs(_t, state):
        x, y = state
        return [alpha * x * (beta - y - x), gamma * x]

    try:
        sol = solve_ivp(
            rhs,
            (0.0, float(max(timepoints))),
            [max(u0, 0.0), 0.0],
            t_eval=np.asarray(timepoints, dtype=float),
            method="RK45",
            rtol=1e-6,
            atol=1e-8,
        )
    except Exception:
        return None
    if not sol.success or sol.y.shape[1] != len(timepoints):
        return None
    pred = sol.y[0]
    if not np.all(np.isfinite(pred)):
        return None
    return pred


def chain_loglik(
    samples: np.ndarray,
    chain_ids: np.ndarray,
    parameter_names: list[str],
    chain_id: int,
    obs_timepoints: np.ndarray,
    obs_values: np.ndarray,
) -> float:
    idx = {name: i for i, name in enumerate(parameter_names)}
    chain_samples = samples[chain_ids == chain_id, :]
    means = {name: float(np.nanmean(chain_samples[:, idx[name]])) for name in ALL_PARAMETERS}
    unique_times = np.array(sorted(pd.unique(obs_timepoints)), dtype=float)
    pred = region_solution(unique_times, means["alpha"], means["beta"], means["gamma"], means["u0"])
    if pred is None:
        return -np.inf
    pred_by_time = dict(zip(unique_times, pred))
    mu = np.array([pred_by_time[float(t)] for t in obs_timepoints], dtype=float)
    sigma = max(means["sigma"], 1e-12)
    return float(np.sum(norm.logpdf(obs_values, loc=mu, scale=sigma)))


def load_region_observations(observations: pd.DataFrame, region: str) -> tuple[np.ndarray, np.ndarray]:
    time_col = observations.columns[1]
    if region not in observations.columns:
        raise KeyError(f"{region} not present in observations")
    sub = observations[[time_col, region]].copy()
    sub[region] = pd.to_numeric(sub[region], errors="coerce")
    sub = sub[np.isfinite(sub[region].to_numpy(dtype=float))]
    return sub[time_col].to_numpy(dtype=float), sub[region].to_numpy(dtype=float)


def process_dataset(project_root: Path, dataset: str, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = project_root / "runs" / "region_rf" / f"copath_{dataset}"
    posterior = pd.read_csv(root / "region_rf_posterior_summary_long.csv")
    diagnostics = pd.read_csv(root / "region_rf_summary.csv")
    observations = pd.read_csv(project_root / "paper-copath" / "data" / OBSERVATION_FILES[dataset])

    main = posterior[posterior["parameter"].isin(MAIN_PARAMETERS)]
    region_max = main.groupby("region_index", as_index=False)["rhat"].max().rename(columns={"rhat": "original_main_max_rhat"})
    nonconverged = set(region_max.loc[region_max["original_main_max_rhat"] > RHAT_CUTOFF, "region_index"])

    adjusted_rows = []
    decision_rows = []
    for _, diag in diagnostics.sort_values("region_index").iterrows():
        region_index = int(diag["region_index"])
        region = str(diag["region"])
        run_id = str(diag["run_id"])
        run_dir = root / run_id
        original_rows = posterior[posterior["region_index"] == region_index].copy()
        original_main_max = float(region_max.loc[region_max["region_index"] == region_index, "original_main_max_rhat"].iloc[0])

        dropped_chain = np.nan
        chain_logliks = {}
        use_original = region_index not in nonconverged
        adjusted = original_rows
        adjusted_main_max = original_main_max

        if not use_original:
            samples, chain_ids, parameter_names = load_h5(run_dir / "posterior.h5")
            obs_t, obs_y = load_region_observations(observations, region)
            for chain_id in sorted(np.unique(chain_ids)):
                chain_logliks[int(chain_id)] = chain_loglik(samples, chain_ids, parameter_names, int(chain_id), obs_t, obs_y)
            dropped_chain = min(chain_logliks, key=chain_logliks.get)
            keep = chain_ids != dropped_chain
            adjusted = summarize_samples(samples[keep, :], chain_ids[keep], parameter_names)
            adjusted["run_id"] = run_id
            adjusted["region_index"] = region_index
            adjusted["region"] = region
            adjusted["rank"] = int(diag["rank"])
            adjusted_main_max = float(adjusted.loc[adjusted["parameter"].isin(MAIN_PARAMETERS), "rhat"].max())

        adjusted_rows.append(adjusted)
        decision_rows.append(
            {
                "dataset": dataset,
                "run_id": run_id,
                "region_index": region_index,
                "region": region,
                "rank": int(diag["rank"]),
                "original_main_max_rhat": original_main_max,
                "adjusted_main_max_rhat": adjusted_main_max,
                "was_nonconverged": region_index in nonconverged,
                "dropped_chain": dropped_chain,
                **{f"chain_{chain_id}_loglik": value for chain_id, value in chain_logliks.items()},
            }
        )

    adjusted_posterior = pd.concat(adjusted_rows, ignore_index=True)
    decisions = pd.DataFrame(decision_rows)
    adjusted_posterior.to_csv(out_dir / f"{dataset}_region_rf_posterior_summary_long.csv", index=False)
    decisions.to_csv(out_dir / f"{dataset}_drop_worst_chain_decisions.csv", index=False)
    return adjusted_posterior, decisions


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    out_dir = project_root / "paper-copath" / "results" / "region_rf_drop_worst_chain"
    out_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for dataset in DATASETS:
        adjusted, decisions = process_dataset(project_root, dataset, out_dir)
        main = adjusted[adjusted["parameter"].isin(MAIN_PARAMETERS)]
        region_max = main.groupby("region_index")["rhat"].max()
        summaries.append(
            {
                "dataset": dataset,
                "regions_adjusted": int(decisions["was_nonconverged"].sum()),
                "regions_passing_after": int((region_max <= RHAT_CUTOFF).sum()),
                "regions_total": int(region_max.shape[0]),
                "max_rhat_after": float(main["rhat"].max()),
            }
        )
    pd.DataFrame(summaries).to_csv(out_dir / "drop_worst_chain_dataset_summary.csv", index=False)
    print(out_dir)


if __name__ == "__main__":
    main()
