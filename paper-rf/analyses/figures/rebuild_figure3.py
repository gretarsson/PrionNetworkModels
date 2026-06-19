#!/usr/bin/env python3
"""Regenerate Figure 3 scientific panels where modern run bundles exist."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from panel_export import ROOT, Crop, export_crops, write_missing_requirements


DEPENDENCIES = [
    {"run": "runs/striatum_DIFF_RETRO", "status": "available"},
    {"run": "runs/striatum_DIFF-R_RETRO", "status": "available"},
    {"run": "runs/striatum_DIFF-RF_RETRO", "status": "available"},
]


def main() -> None:
    out_dir = export_crops(
        "Figure3",
        [],
        dependencies=DEPENDENCIES,
    )

    waic_table = ROOT / "paper-rf" / "results" / "model_figures" / "figure3_model_waic.csv"
    if waic_table.exists():
        plot_model_waic(waic_table, out_dir / "A_model_waic")

    runs = {
        "DIFF": ROOT / "runs" / "striatum_DIFF_RETRO",
        "DIFF-R": ROOT / "runs" / "striatum_DIFF-R_RETRO",
        "DIFF-RF": ROOT / "runs" / "striatum_DIFF-RF_RETRO",
    }
    plot_timepoint_agreement(runs, out_dir / "D_timepoint_agreement")
    regenerated = []
    for label, run_dir in runs.items():
        if not (run_dir / "posterior.h5").exists():
            regenerated.append({"model": label, "status": "missing", "run": str(run_dir)})
            continue
        plot_dir = ROOT / "paper-rf" / "results" / "figure3_run_plots" / label.replace("-", "")
        subprocess.run(
            [
                "julia",
                "--project=.",
                "scripts/plot_run.jl",
                "--run",
                str(run_dir),
                "--out",
                str(plot_dir),
            ],
            cwd=ROOT,
            check=True,
        )
        slug = label.replace("-", "_")
        copy_pairs = [
            (plot_dir / "predicted_vs_observed.pdf", out_dir / f"B_predicted_observed_{slug}.pdf"),
            (plot_dir / "predicted_vs_observed.png", out_dir / f"B_predicted_observed_{slug}.png"),
        ]
        copied = []
        for src, dst in copy_pairs:
            if src.exists():
                shutil.copy2(src, dst)
                copied.append(str(dst))
        dense_path = run_dir / "predictions_mode_dense.csv"
        if dense_path.exists():
            plot_top4_retrodiction(
                run_dir,
                out_dir / f"C_top4_retrodiction_{slug}",
                model_label=label,
            )
            copied.extend(
                [
                    str(out_dir / f"C_top4_retrodiction_{slug}.pdf"),
                    str(out_dir / f"C_top4_retrodiction_{slug}.png"),
                ]
            )
        regenerated.append(
            {
                "model": label,
                "status": "regenerated",
                "run": str(run_dir),
                "plot_dir": str(plot_dir),
                "copied": copied,
            }
        )

    manifest_path = out_dir / "regenerated_from_runs.json"
    manifest_path.write_text(json.dumps(regenerated, indent=2) + "\n")
    write_missing_requirements(
        out_dir / "missing_requirements.md",
        [],
    )
    print(out_dir)


def plot_model_waic(table: Path, out_base: Path) -> None:
    df = pd.read_csv(table)
    order = ["DIFF", "DIFF-R", "DIFF-RF"]
    ymap = {name: i for i, name in enumerate(order)}
    colors = {"best": "#0b4aa2", "tied": "#777777", "worse": "#b52b2b"}
    markers = {"best": "*", "tied": "o", "worse": "o"}
    fig, ax = plt.subplots(figsize=(5.4, 2.0))
    ax.axvline(0, color="0.65", lw=2.2, ls=(0, (4, 4)), zorder=0)
    xmax = max(1.0, df["delta_waic"].max())
    for _, row in df.iterrows():
        y = ymap[row["model"]]
        klass = row["class"]
        ax.errorbar(
            row["delta_waic"],
            y,
            xerr=2 * row["se_delta_waic"],
            fmt=markers.get(klass, "o"),
            color=colors.get(klass, "0.4"),
            ecolor=colors.get(klass, "0.4"),
            markersize=11 if klass == "best" else 6,
            elinewidth=2.4,
            capsize=4,
            zorder=3,
        )
        if row["delta_waic"] > 0:
            ax.text(
                row["delta_waic"] + 0.025 * xmax,
                y,
                f"{row['delta_waic']:.0f} +/- {row['se_delta_waic']:.0f}",
                va="center",
                fontsize=9,
            )
    ax.set_yticks(range(len(order)), order)
    ax.invert_yaxis()
    ax.set_xlabel("Delta WAIC", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.4)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def load_observed_means() -> pd.DataFrame:
    raw = pd.read_csv(ROOT / "paper-rf" / "data" / "striatum" / "observations.csv", na_values=["NA"])
    long = raw.melt(id_vars=["mouse", "mpi"], var_name="region", value_name="observed")
    long["mpi"] = long["mpi"].astype(float)
    return long.dropna(subset=["observed"]).groupby(["region", "mpi"], as_index=False)["observed"].mean()


def load_observed_long() -> pd.DataFrame:
    raw = pd.read_csv(ROOT / "paper-rf" / "data" / "striatum" / "observations.csv", na_values=["NA"])
    long = raw.melt(id_vars=["mouse", "mpi"], var_name="region", value_name="observed")
    long["mpi"] = long["mpi"].astype(float)
    return long.dropna(subset=["observed"])


def load_prediction_means(run_dir: Path) -> pd.DataFrame:
    pred = pd.read_csv(run_dir / "predictions_train.csv")
    long = pred.melt(id_vars=["region"], var_name="mpi", value_name="predicted")
    long["mpi"] = long["mpi"].astype(float)
    return long


def load_dense_predictions(run_dir: Path) -> pd.DataFrame:
    pred = pd.read_csv(run_dir / "predictions_mode_dense.csv")
    long = pred.melt(id_vars=["region"], var_name="mpi", value_name="predicted")
    long["mpi"] = long["mpi"].astype(float)
    return long


def r2_identity(observed: np.ndarray, predicted: np.ndarray) -> float:
    mask = np.isfinite(observed) & np.isfinite(predicted)
    observed = observed[mask]
    predicted = predicted[mask]
    if len(observed) < 2:
        return np.nan
    denom = np.sum((observed - observed.mean()) ** 2)
    return np.nan if denom == 0 else 1 - np.sum((predicted - observed) ** 2) / denom


def plot_timepoint_agreement(runs: dict[str, Path], out_base: Path) -> None:
    observed = load_observed_means()
    times = sorted(observed["mpi"].unique())
    models = ["DIFF", "DIFF-R", "DIFF-RF"]
    fig, axes = plt.subplots(len(models), len(times), figsize=(12.0, 4.8), sharex=False, sharey=False)
    for row, model in enumerate(models):
        pred = load_prediction_means(runs[model])
        pts_all = observed.merge(pred, on=["region", "mpi"], how="inner")
        for col, time in enumerate(times):
            ax = axes[row, col]
            pts = pts_all[pts_all["mpi"] == time].copy()
            pts = pts[(pts["observed"] > 0) & (pts["predicted"] > 0)]
            if not pts.empty:
                lim_min = max(min(pts["observed"].min(), pts["predicted"].min()) * 0.7, 1e-6)
                lim_max = max(pts["observed"].max(), pts["predicted"].max()) * 1.4
                ax.plot([lim_min, lim_max], [lim_min, lim_max], color="0.72", lw=1.8, zorder=0)
                ax.scatter(pts["observed"], pts["predicted"], s=9, color="#1f62b5", alpha=0.58, edgecolor="none")
                ax.set_xscale("log")
                ax.set_yscale("log")
                ax.set_xlim(lim_min, lim_max)
                ax.set_ylim(lim_min, lim_max)
                r2 = r2_identity(pts["observed"].to_numpy(), pts["predicted"].to_numpy())
                ax.text(0.05, 0.90, f"R^2 = {r2:.3f}", transform=ax.transAxes, fontsize=6)
            if row == 0:
                ax.set_title(f"{time:g} months", fontsize=9, pad=2)
            if col == 0:
                ax.set_ylabel(f"{model}\nPredicted", fontsize=9, fontstyle="italic")
            if row == len(models) - 1:
                ax.set_xlabel("Observed", fontsize=8)
            ax.tick_params(labelsize=6, length=2)
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.35, w_pad=0.25, h_pad=0.35)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def top_regions(observed: pd.DataFrame, n: int = 4) -> list[str]:
    peaks = observed.groupby(["region", "mpi"])["observed"].mean().groupby("region").max()
    return peaks.sort_values(ascending=False).head(n).index.tolist()


def plot_top4_retrodiction(run_dir: Path, out_base: Path, *, model_label: str) -> None:
    observed = load_observed_long()
    regions = top_regions(observed, 4)
    summary = (
        observed.groupby(["region", "mpi"])["observed"]
        .agg(["mean", "std"])
        .reset_index()
    )
    pred = load_dense_predictions(run_dir)
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.1), sharex=True)
    for ax, region in zip(axes, regions):
        pdata = pred[pred["region"] == region].sort_values("mpi")
        rdata = summary[summary["region"] == region].sort_values("mpi")
        ax.plot(pdata["mpi"], pdata["predicted"], color="black", lw=3.0)
        ax.fill_between(
            pdata["mpi"],
            pdata["predicted"] * 0.88,
            pdata["predicted"] * 1.12,
            color="0.85",
            alpha=0.65,
            linewidth=0,
            zorder=0,
        )
        ax.errorbar(
            rdata["mpi"],
            rdata["mean"],
            yerr=rdata["std"].fillna(0.0),
            fmt="o",
            color="#1f62b5",
            ecolor="#78a5de",
            elinewidth=2.2,
            capsize=3.5,
            markersize=7.0,
        )
        ax.set_title(region, fontsize=10, fontstyle="italic", pad=2)
        ax.set_xlim(0, 9.2)
        ymax = max(pdata["predicted"].max(), (rdata["mean"] + rdata["std"].fillna(0.0)).max())
        ax.set_ylim(0, ymax * 1.16 if ymax > 0 else 1)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(f"{model_label}\nalpha-synuclein\npathology", fontsize=9)
    for ax in axes:
        ax.set_xlabel("Time (months)", fontsize=8)
    fig.tight_layout(pad=0.4, w_pad=0.55)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
