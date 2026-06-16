#!/usr/bin/env python3
"""Prepare REGION-RF-like summary tables from iterative chain-drop outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


DATASETS = ["syn_app", "syn_mapt", "tau_app", "tau_mapt"]
PARAMETERS = ["alpha", "beta", "gamma", "u0", "sigma"]


def prepare_dataset(project_root: Path, dataset: str, adjusted_dir: Path) -> None:
    original_root = project_root / "runs" / "region_rf" / f"copath_{dataset}"
    original_summary = pd.read_csv(original_root / "region_rf_summary.csv")
    adjusted_long = pd.read_csv(adjusted_dir / f"{dataset}_region_rf_posterior_summary_long.csv")
    decisions = pd.read_csv(adjusted_dir / f"{dataset}_iterative_chain_drop_decisions.csv")

    means = (
        adjusted_long[adjusted_long["parameter"].isin(PARAMETERS)]
        .pivot_table(index="region_index", columns="parameter", values="mean", aggfunc="first")
        .reset_index()
    )
    out = original_summary.drop(columns=[p for p in PARAMETERS if p in original_summary.columns])
    out = out.merge(means, on="region_index", how="left")
    out = out.merge(
        decisions[
            [
                "region_index",
                "original_main_max_rhat",
                "final_main_max_rhat",
                "status",
                "n_dropped",
                "dropped_chains",
                "retained_chains",
            ]
        ],
        on="region_index",
        how="left",
    )
    out.to_csv(adjusted_dir / f"{dataset}_region_rf_summary.csv", index=False)


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    adjusted_dir = project_root / "paper-copath" / "results" / "region_rf_iterative_drop_low_likelihood_chains"
    for dataset in DATASETS:
        prepare_dataset(project_root, dataset, adjusted_dir)
    print(adjusted_dir)


if __name__ == "__main__":
    main()
