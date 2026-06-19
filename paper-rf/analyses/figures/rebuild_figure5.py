#!/usr/bin/env python3
"""Regenerate Figure 5 scientific panels into `paper-rf/figures/Figure5/`."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from panel_export import ROOT, write_missing_requirements


RUNS = {
    "DIFF-R": ROOT / "runs" / "striatum_DIFF-R_RETRO_T-1",
    "DIFF-RF": ROOT / "runs" / "striatum_DIFF-RF_RETRO_T-1",
}
OBS_PATH = ROOT / "paper-rf" / "data" / "striatum" / "observations.csv"
TRAIN_TIMES = [0.1, 0.2, 0.3, 0.5, 1.0, 3.0, 6.0]
HELDOUT_TIME = 9.0
BLUE = "#1f62b5"
RED = "#d55a44"


def load_observed() -> pd.DataFrame:
    raw = pd.read_csv(OBS_PATH, na_values=["NA"])
    long = raw.melt(id_vars=["mouse", "mpi"], var_name="region", value_name="observed")
    long["mpi"] = long["mpi"].astype(float)
    return long.dropna(subset=["observed"])


def load_predictions(run_dir: Path) -> pd.DataFrame:
    pred = pd.read_csv(run_dir / "predictions_train.csv")
    long = pred.melt(id_vars=["region"], var_name="mpi", value_name="predicted")
    long["mpi"] = long["mpi"].astype(float)
    return long


def load_dense_predictions(run_dir: Path) -> pd.DataFrame:
    pred = pd.read_csv(run_dir / "predictions_mode_dense.csv")
    long = pred.melt(id_vars=["region"], var_name="mpi", value_name="predicted")
    long["mpi"] = long["mpi"].astype(float)
    return long


def r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    mask = np.isfinite(observed) & np.isfinite(predicted)
    observed = observed[mask]
    predicted = predicted[mask]
    if len(observed) < 2:
        return np.nan
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - observed.mean()) ** 2)
    return 1 - ss_res / ss_tot


def merged_points(observed: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    obs_mean = observed.groupby(["region", "mpi"], as_index=False)["observed"].mean()
    return obs_mean.merge(predictions, on=["region", "mpi"], how="inner")


def plot_predicted_observed(out_base: Path, observed: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(9.4, 2.35), sharex=False, sharey=False)
    panel_specs = [
        ("DIFF-R", "in sample", TRAIN_TIMES, "o", BLUE),
        ("DIFF-R", "out of sample", [HELDOUT_TIME], "^", RED),
        ("DIFF-RF", "in sample", TRAIN_TIMES, "o", BLUE),
        ("DIFF-RF", "out of sample", [HELDOUT_TIME], "^", RED),
    ]
    for ax, (model, title, times, marker, color) in zip(axes, panel_specs):
        pred = load_predictions(RUNS[model])
        pts = merged_points(observed[observed["mpi"].isin(times)], pred)
        r2 = r_squared(pts["observed"].to_numpy(), pts["predicted"].to_numpy())
        vmax = max(pts["observed"].max(), pts["predicted"].max())
        lim = max(0.02, vmax * 1.05)
        ax.plot([0, lim], [0, lim], color="0.7", lw=3.0, zorder=0)
        ax.scatter(
            pts["observed"],
            pts["predicted"],
            s=22 if marker == "o" else 34,
            marker=marker,
            color=color,
            alpha=0.72,
            edgecolor="none",
        )
        ax.text(0.04 * lim, 0.90 * lim, f"R^2 = {r2:.3f}", fontsize=9)
        ax.set_title(title, fontsize=11, pad=2)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Predicted", fontsize=10)
    axes[2].set_ylabel("Predicted", fontsize=10)
    for ax in axes:
        ax.set_xlabel("Observed", fontsize=10)
    axes[0].text(0.5, 1.20, "DIFF-R", transform=axes[0].transAxes, fontsize=13, fontstyle="italic")
    axes[2].text(0.5, 1.20, "DIFF-RF", transform=axes[2].transAxes, fontsize=13, fontstyle="italic")
    fig.tight_layout(pad=0.4, w_pad=0.8)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def top_regions(observed: pd.DataFrame, n: int = 4) -> list[str]:
    peaks = observed.groupby(["region", "mpi"])["observed"].mean().groupby("region").max()
    return peaks.sort_values(ascending=False).head(n).index.tolist()


def plot_trajectories(out_base: Path, observed: pd.DataFrame) -> None:
    regions = top_regions(observed, 4)
    summary = (
        observed.groupby(["region", "mpi"])["observed"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    fig, axes = plt.subplots(2, 4, figsize=(9.4, 3.6), sharex=True)
    for row, model in enumerate(["DIFF-R", "DIFF-RF"]):
        pred = load_dense_predictions(RUNS[model])
        for col, region in enumerate(regions):
            ax = axes[row, col]
            pdata = pred[pred["region"] == region].sort_values("mpi")
            rdata = summary[summary["region"] == region].sort_values("mpi")
            train = rdata[rdata["mpi"].isin(TRAIN_TIMES)]
            held = rdata[rdata["mpi"] == HELDOUT_TIME]
            ax.plot(pdata["mpi"], pdata["predicted"], color="black", lw=2.5)
            ax.errorbar(
                train["mpi"],
                train["mean"],
                yerr=train["std"].fillna(0.0),
                fmt="o",
                color=BLUE,
                ecolor="#78a5de",
                elinewidth=2.0,
                capsize=3,
                markersize=6.5,
            )
            ax.errorbar(
                held["mpi"],
                held["mean"],
                yerr=held["std"].fillna(0.0),
                fmt="^",
                color=RED,
                ecolor="#e59b8d",
                elinewidth=2.0,
                capsize=3,
                markersize=7.0,
            )
            ax.fill_between(
                pdata["mpi"],
                pdata["predicted"] * 0.88,
                pdata["predicted"] * 1.12,
                color="0.85",
                alpha=0.65,
                linewidth=0,
                zorder=0,
            )
            if row == 0:
                ax.set_title(region, fontsize=10, fontstyle="italic", pad=2)
            if col == 0:
                ax.set_ylabel(f"{model}\nalpha-synuclein\npathology", fontsize=9)
            ax.set_xlim(0, 9.2)
            ymax = max(pdata["predicted"].max(), rdata["mean"].max() + rdata["std"].fillna(0.0).max())
            ax.set_ylim(0, ymax * 1.18 if ymax > 0 else 1)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel("Time (months)", fontsize=9)
    fig.tight_layout(pad=0.45, w_pad=0.65, h_pad=0.45)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    out_dir = ROOT / "paper-rf" / "figures" / "Figure5"
    out_dir.mkdir(parents=True, exist_ok=True)
    required = [
        *(run / "predictions_train.csv" for run in RUNS.values()),
        *(run / "predictions_mode_dense.csv" for run in RUNS.values()),
        OBS_PATH,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        write_missing_requirements(out_dir / "missing_requirements.md", missing)
        raise FileNotFoundError("Missing inputs for Figure 5:\n" + "\n".join(missing))
    observed = load_observed()
    plot_predicted_observed(out_dir / "A_heldout_predicted_observed", observed)
    plot_trajectories(out_dir / "B_heldout_regional_trajectories", observed)
    write_missing_requirements(out_dir / "missing_requirements.md", [])
    print(out_dir)


if __name__ == "__main__":
    main()
