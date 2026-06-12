#!/usr/bin/env python3
"""Plot chain-wise PC1 directions for no-seed DIFF-RF biological axes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FILTER_ORDER = ["all", "beta_positive", "updated"]
FILTER_LABELS = {
    "all": "All regions",
    "beta_positive": "Beta-positive",
    "updated": "Updated",
}

CHAIN_COLORS = {
    "C1": "#1f77b4",
    "C2": "#2ca02c",
    "C3": "#d62728",
    "C4": "#9467bd",
    "C5": "#ff7f0e",
    "C6": "#17becf",
    "C7": "#8c564b",
    "C8": "#e377c2",
}


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 300,
        "savefig.dpi": 300,
    })


def arrow(ax, x: float, y: float, color: str, alpha: float, linestyle: str, linewidth: float) -> None:
    ax.annotate(
        "",
        xy=(x, y),
        xytext=(0, 0),
        arrowprops={
            "arrowstyle": "-|>",
            "color": color,
            "alpha": alpha,
            "lw": linewidth,
            "linestyle": linestyle,
            "mutation_scale": 14,
            "shrinkA": 0,
            "shrinkB": 0,
        },
    )


def plot_panel(ax, df: pd.DataFrame, level: str, scale_by_loglik: bool) -> None:
    theta = np.linspace(0, 2 * np.pi, 256)
    ax.plot(np.cos(theta), np.sin(theta), color="#D0D0D0", lw=0.9, zorder=0)
    ax.axhline(0, color="#B8B8B8", lw=0.8, zorder=0)
    ax.axvline(0, color="#B8B8B8", lw=0.8, zorder=0)

    level_df = df[df["filter_level"] == level].copy()
    for chain in sorted(level_df["chain"].unique()):
        color = CHAIN_COLORS.get(chain, "#4C4C4C")
        chain_df = level_df[level_df["chain"] == chain]

        for _, row in chain_df[chain_df["dataset"] == "striatum"].iterrows():
            weight = float(row.get("loglik_weight", 1.0)) if scale_by_loglik else 1.0
            x = float(row["pc1_loading_beta"]) * weight
            y = float(row["pc1_loading_gamma"]) * weight
            arrow(ax, x, y, color=color, alpha=0.30, linestyle=(0, (2, 2)), linewidth=2.0)
            ax.scatter([x], [y], s=18, marker="s", color=color, alpha=0.45, edgecolor="none", zorder=3)

        for _, row in chain_df[chain_df["dataset"] == "hippocampus"].iterrows():
            weight = float(row.get("loglik_weight", 1.0)) if scale_by_loglik else 1.0
            x = float(row["pc1_loading_beta"]) * weight
            y = float(row["pc1_loading_gamma"]) * weight
            arrow(ax, x, y, color=color, alpha=0.86, linestyle="-", linewidth=2.2)
            ax.scatter([x], [y], s=24, marker="o", color=color, alpha=0.95, edgecolor="white", linewidth=0.4, zorder=4)
            label_radius = max(weight, 0.12)
            raw_x = float(row["pc1_loading_beta"]) * label_radius
            raw_y = float(row["pc1_loading_gamma"]) * label_radius
            ax.text(raw_x * 1.06, raw_y * 1.06, chain, color=color, ha="center", va="center", fontsize=9)

    ax.set_title(FILTER_LABELS[level])
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("PC1 loading on beta")
    ax.grid(True, color="#E8E8E8", lw=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        default="paper-rf/results/sidequests/no_seed_diff_rf/chainwise/chainwise_axis_summary.csv",
        help="chainwise_axis_summary.csv",
    )
    parser.add_argument(
        "--out",
        default="paper-rf/figures/no_seed_diff_rf/chainwise/pc1_direction_by_filter",
        help="Output path without extension.",
    )
    parser.add_argument(
        "--weights",
        help="Optional CSV with chain,dataset,loglik_weight columns. Arrows are scaled by loglik_weight.",
    )
    args = parser.parse_args()

    setup_style()
    df = pd.read_csv(args.summary)
    scale_by_loglik = args.weights is not None
    if args.weights is not None:
        weights = pd.read_csv(args.weights)
        keep_cols = ["chain", "dataset", "loglik_all", "loglik_weight"]
        df = df.merge(weights[keep_cols], on=["chain", "dataset"], how="left")
        df["loglik_weight"] = df["loglik_weight"].fillna(1.0)

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8), sharey=True)
    for ax, level in zip(axes, FILTER_ORDER):
        plot_panel(ax, df, level, scale_by_loglik)
    axes[0].set_ylabel("PC1 loading on gamma")

    chains = [chain for chain in CHAIN_COLORS if chain in set(df["chain"])]
    chain_handles = [
        plt.Line2D([0], [0], color=color, lw=2.4, label=chain)
        for chain, color in CHAIN_COLORS.items()
        if chain in chains
    ]
    dataset_handles = [
        plt.Line2D([0], [0], color="#444444", lw=2.2, marker="o", label="Hippocampus"),
        plt.Line2D([0], [0], color="#444444", lw=2.0, linestyle=(0, (2, 2)), marker="s", alpha=0.45, label="Striatum"),
    ]
    fig.legend(
        handles=chain_handles + dataset_handles,
        loc="lower center",
        ncol=min(10, len(chain_handles) + len(dataset_handles)),
        frameon=False,
        bbox_to_anchor=(0.5, -0.03),
        handlelength=2.2,
        columnspacing=1.4,
    )
    fig.tight_layout(rect=[0, 0.08, 1, 1])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
