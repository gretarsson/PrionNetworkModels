#!/usr/bin/env python3
"""Regenerate Figure 2 scientific panels into `paper-rf/figures/Figure2/`."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from panel_export import ROOT, Crop, export_crops, write_missing_requirements


DEPENDENCIES = [
    {"run": "runs/striatum_DIFF_RETRO and transport variants", "status": "missing"},
    {"run": "runs/striatum_DIFF-R_RETRO and transport variants", "status": "missing"},
    {"run": "runs/striatum_DIFF-RF_RETRO", "status": "partial"},
]


def main() -> None:
    out_dir = ROOT / "paper-rf" / "figures" / "Figure2"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = ROOT / "paper-rf" / "results" / "model_figures" / "figure2_transport_waic.csv"
    if table.exists():
        plot_transport_waic(table, out_dir / "transport_waic")
        write_missing_requirements(out_dir / "missing_requirements.md", [])
    else:
        out_dir = export_crops(
            "Figure2",
            [
                Crop(
                    "transport_waic",
                    "transport.png",
                    1200,
                    0,
                    1567,
                    1633,
                    "WAIC comparison for Euclidean, anterograde, retrograde, and bidirectional transport.",
                    "Temporary fallback until figure2_transport_waic.csv is generated.",
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


def plot_transport_waic(table: Path, out_base: Path) -> None:
    df = pd.read_csv(table)
    panels = ["DIFF", "DIFF-R", "DIFF-RF"]
    order = {
        "DIFF": ["anterograde", "euclidean", "retrograde", "bidirectional"],
        "DIFF-R": ["anterograde", "euclidean", "retrograde", "bidirectional"],
        "DIFF-RF": ["anterograde", "euclidean", "bidirectional", "retrograde"],
    }
    colors = {"best": "#0b4aa2", "tied": "#777777", "worse": "#b52b2b"}
    markers = {"best": "*", "tied": "o", "worse": "o"}

    fig, axes = plt.subplots(3, 1, figsize=(6.4, 5.6), sharex=False)
    for ax, panel in zip(axes, panels):
        sub = df[df["panel"] == panel].copy()
        ylabels = order[panel]
        ymap = {name: i for i, name in enumerate(ylabels)}
        ax.axvline(0, color="0.65", lw=2.2, ls=(0, (4, 4)), zorder=0)
        for _, row in sub.iterrows():
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
                    row["delta_waic"] + 0.02 * max(1, sub["delta_waic"].max()),
                    y,
                    f"{row['delta_waic']:.0f} +/- {row['se_delta_waic']:.0f}",
                    va="center",
                    fontsize=9,
                )
        ax.set_title(panel, loc="left", fontsize=13, pad=2)
        ax.set_yticks(range(len(ylabels)), ylabels)
        ax.invert_yaxis()
        ax.tick_params(labelsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Delta WAIC", fontsize=12)
    fig.tight_layout(pad=0.4)
    for ext in ("pdf", "png"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
