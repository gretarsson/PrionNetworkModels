#!/usr/bin/env python3
"""Regenerate Figure 4 scientific panels into `paper-rf/figures/Figure4/`."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from panel_export import ROOT, Crop, export_crops, write_missing_requirements


DEPENDENCIES = [
    {"run": "runs/striatum_DIFF-RF_RETRO", "status": "available"},
    {"run": "runs/striatum_DIFF-RF_RETRO_connectivity_nulls", "status": "missing"},
    {"run": "runs/striatum_DIFF-RF_RETRO_seed_nulls", "status": "missing"},
]


def main() -> None:
    out_dir = ROOT / "paper-rf" / "figures" / "Figure4"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = ROOT / "paper-rf" / "results" / "model_figures" / "figure4_null_waic.csv"
    if table.exists():
        plot_null_waic(table, out_dir / "null_waic")
        write_missing_requirements(out_dir / "missing_requirements.md", [])
    else:
        out_dir = export_crops(
            "Figure4",
            [
                Crop(
                    "null_waic",
                    "nulls.png",
                    820,
                    0,
                    1730,
                    1183,
                    "WAIC distributions for connectivity and seeding null models.",
                    "Temporary fallback until figure4_null_waic.csv is generated.",
                ),
            ],
            dependencies=DEPENDENCIES,
        )
        write_missing_requirements(
            out_dir / "missing_requirements.md",
            [
                "Run paper-rf/analyses/figures/import_synuclein_spread_runs.jl once, then paper-rf/analyses/figures/export_run_bundle_waic.jl.",
            ],
        )
    print(out_dir)


def plot_null_waic(table: Path, out_base: Path) -> None:
    df = pd.read_csv(table)
    order = ["connectivity", "seeding"]
    labels = ["connectivity null", "seeding null"]
    rng = np.random.default_rng(123456)
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    null_color = "0.35"
    empirical_color = "#0b4aa2"

    data = [
        df[(df["null_type"] == kind) & (~df["is_empirical"]) & (df["waic"] < -58000)]["waic"].to_numpy()
        for kind in order
    ]
    ax.boxplot(
        data,
        positions=[1, 2],
        widths=0.35,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": null_color, "edgecolor": null_color, "alpha": 0.8},
        medianprops={"color": "black", "linewidth": 2.5},
        whiskerprops={"color": "black", "linewidth": 2.2},
        capprops={"color": "black", "linewidth": 2.2},
    )
    for i, values in enumerate(data, start=1):
        jitter = rng.normal(0, 0.045, len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, s=28, color="black", alpha=0.6, zorder=2)
    for i, kind in enumerate(order, start=1):
        empirical = df[(df["null_type"] == kind) & (df["is_empirical"])]["waic"].to_numpy()
        if len(empirical):
            ax.scatter([i], [empirical[0]], s=85, color=empirical_color, zorder=5)
    ax.set_xticks([1, 2], labels)
    ax.set_ylabel("WAIC", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.5)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
