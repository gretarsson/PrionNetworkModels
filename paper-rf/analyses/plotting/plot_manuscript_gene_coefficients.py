#!/usr/bin/env python3
"""Create manuscript-style gene coefficient PCA panels.

This script is intentionally narrow: it makes the gene coefficient panels used
for the rise-fall manuscript from an already-generated paper-rf filter level.
The default filter is beta_positive, which corresponds to the updated primary
analysis choice for the manuscript refresh.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


DATASETS = ("striatum", "hippocampus")


@dataclass(frozen=True)
class GeneCoefficientData:
    name: str
    coefficients: pd.DataFrame
    correlations: pd.DataFrame
    pca_summary: pd.DataFrame
    pc1: np.ndarray
    pc2: np.ndarray
    pc1_var: float
    pc2_var: float
    n_regions: int
    n_genes: int
    center: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        default="paper-rf/results/filtering",
        help="Root containing filter-level paper-rf results.",
    )
    parser.add_argument(
        "--filter-level",
        default="beta_positive",
        help="Filter level to plot, e.g. all, beta_positive, or updated.",
    )
    parser.add_argument(
        "--out-dir",
        default="paper-rf/figures/manuscript/gene_coefficients_beta_positive",
        help="Directory for generated figures.",
    )
    parser.add_argument(
        "--stats-dir",
        default="paper-rf/results/manuscript/gene_coefficients_beta_positive",
        help="Directory for generated summary CSVs.",
    )
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.labelsize": 22,
            "axes.titlesize": 15,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 17,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#d9d9d9",
            "grid.linewidth": 0.45,
            "grid.alpha": 0.55,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_dataset(results_root: Path, filter_level: str, dataset: str) -> GeneCoefficientData:
    data_dir = results_root / filter_level / "transcriptomics" / dataset
    coefficients = pd.read_csv(data_dir / "gene_parameter_coefficients.csv")
    correlations = pd.read_csv(data_dir / "gene_eta_correlations.csv")
    pca_summary = pd.read_csv(data_dir / "pca_summary.csv")

    pc1_row = pca_summary.loc[pca_summary["component"] == "PC1"].iloc[0]
    pc2_row = pca_summary.loc[pca_summary["component"] == "PC2"].iloc[0]
    pc1 = np.array([pc1_row["loading_beta"], pc1_row["loading_gamma"]], dtype=float)
    pc2 = np.array([pc2_row["loading_beta"], pc2_row["loading_gamma"]], dtype=float)

    # Keep the same orientation convention as the old manuscript plot: the
    # dominant fall-aligned PC1 points upward in gamma.
    if pc1[1] < 0:
        pc1 = -pc1
    if np.linalg.det(np.vstack([pc1, pc2])) < 0:
        pc2 = -pc2

    return GeneCoefficientData(
        name=dataset,
        coefficients=coefficients,
        correlations=correlations,
        pca_summary=pca_summary,
        pc1=pc1,
        pc2=pc2,
        pc1_var=float(pc1_row["explained_variance_ratio"]),
        pc2_var=float(pc2_row["explained_variance_ratio"]),
        n_regions=int(pc1_row["n_regions"]),
        n_genes=int(pc1_row["n_genes"]),
        center=coefficients[["coef_beta", "coef_gamma"]].mean(axis=0).to_numpy(dtype=float),
    )


def merged_plot_table(data: GeneCoefficientData) -> pd.DataFrame:
    return data.coefficients.merge(
        data.correlations[["gene", "r"]],
        on="gene",
        how="left",
        validate="one_to_one",
    )


def pca_segment(center: np.ndarray, vec: np.ndarray, x: np.ndarray, y: np.ndarray, scale: float) -> tuple[np.ndarray, np.ndarray]:
    xrange = float(np.nanmax(x) - np.nanmin(x))
    yrange = float(np.nanmax(y) - np.nanmin(y))
    length = scale * min(
        xrange / max(abs(vec[0]), np.finfo(float).eps),
        yrange / max(abs(vec[1]), np.finfo(float).eps),
    )
    return center - length * vec, center + length * vec


def axis_limits(values: np.ndarray, pad_frac: float = 0.06) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    pad = pad_frac * max(hi - lo, np.finfo(float).eps)
    return lo - pad, hi + pad


def draw_panel(
    ax: plt.Axes,
    data: GeneCoefficientData,
    *,
    title: str | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    show_ylabel: bool = True,
    add_colorbar: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
) -> object:
    table = merged_plot_table(data)
    x = table["coef_beta"].to_numpy(dtype=float)
    y = table["coef_gamma"].to_numpy(dtype=float)
    color = table["r"].to_numpy(dtype=float)
    color_values = color[np.isfinite(color)]
    if vmin is None:
        vmin = float(np.min(color_values))
    if vmax is None:
        vmax = float(np.max(color_values))

    scatter = ax.scatter(
        x,
        y,
        c=color,
        cmap="coolwarm",
        vmin=vmin,
        vmax=vmax,
        s=72,
        alpha=0.60,
        linewidths=0.45,
        edgecolors="#4d4d4d",
        rasterized=True,
    )

    p2a, p2b = pca_segment(data.center, data.pc2, x, y, 0.26)
    p1a, p1b = pca_segment(data.center, data.pc1, x, y, 0.34)
    ax.plot([p2a[0], p2b[0]], [p2a[1], p2b[1]], color="#8c8c8c", lw=9, solid_capstyle="round")
    ax.plot([p1a[0], p1b[0]], [p1a[1], p1b[1]], color="black", lw=9, solid_capstyle="round")

    ax.set_xlabel(r"gene coefficient for $z(\beta)$")
    if show_ylabel:
        ax.set_ylabel(r"gene coefficient for $z(\gamma)$")
    if title:
        ax.set_title(title, loc="left", fontweight="bold", pad=6)
    ax.set_xlim(xlim or axis_limits(x))
    ax.set_ylim(ylim or axis_limits(y))

    handles = [
        Line2D([0], [0], color="#8c8c8c", lw=7, solid_capstyle="round", label="PC2"),
        Line2D([0], [0], color="black", lw=7, solid_capstyle="round", label="PC1"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, handlelength=1.35, handletextpad=0.15)

    ax.text(
        0.98,
        0.97,
        f"PC1 {100 * data.pc1_var:.1f}%\nPC2 {100 * data.pc2_var:.1f}%\nregions n = {data.n_regions}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "none", "alpha": 0.72},
    )

    if add_colorbar:
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.055, pad=0.03)
        cbar.set_label(r"corr(gene, $\eta$)", rotation=90, labelpad=15, fontsize=22)
        cbar.ax.tick_params(labelsize=15)
    return scatter


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    fig.savefig(out_dir / f"{stem}.png", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_stats(datasets: list[GeneCoefficientData], stats_dir: Path) -> None:
    rows = []
    for data in datasets:
        rows.append(
            {
                "dataset": data.name,
                "n_regions": data.n_regions,
                "n_genes": data.n_genes,
                "pc1_loading_beta": data.pc1[0],
                "pc1_loading_gamma": data.pc1[1],
                "pc2_loading_beta": data.pc2[0],
                "pc2_loading_gamma": data.pc2[1],
                "pc1_explained_variance": data.pc1_var,
                "pc2_explained_variance": data.pc2_var,
            }
        )
    pd.DataFrame(rows).to_csv(stats_dir / "gene_coefficient_pca_stats.csv", index=False)

    if len(datasets) == 2:
        a, b = datasets
        cosine = float(np.dot(a.pc1, b.pc1) / (np.linalg.norm(a.pc1) * np.linalg.norm(b.pc1)))
        cosine = max(-1.0, min(1.0, cosine))
        angle = math.degrees(math.acos(cosine))
        corr_table = a.correlations[["gene", "r"]].merge(
            b.correlations[["gene", "r"]],
            on="gene",
            suffixes=("_striatum", "_hippocampus"),
        )
        gene_r = float(corr_table["r_striatum"].corr(corr_table["r_hippocampus"]))
        pd.DataFrame(
            [
                {
                    "dataset_a": a.name,
                    "dataset_b": b.name,
                    "pc1_cosine_similarity": cosine,
                    "pc1_angle_degrees": angle,
                    "gene_eta_correlation_pearson": gene_r,
                    "n_shared_genes": len(corr_table),
                }
            ]
        ).to_csv(stats_dir / "striatum_hippocampus_pc1_comparison.csv", index=False)


def plot_pc1_comparison(datasets: list[GeneCoefficientData], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 4.5), constrained_layout=True)
    colors = {"striatum": "#1f77b4", "hippocampus": "#d62728"}
    for data in datasets:
        vec = data.pc1 / np.linalg.norm(data.pc1)
        ax.arrow(
            0,
            0,
            vec[0],
            vec[1],
            width=0.010,
            head_width=0.055,
            length_includes_head=True,
            color=colors.get(data.name, "black"),
            label=f"{data.name.capitalize()} ({100 * data.pc1_var:.1f}%, n={data.n_regions})",
        )
    ax.axhline(0, color="#bdbdbd", lw=0.8)
    ax.axvline(0, color="#bdbdbd", lw=0.8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.15, 1.05)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"coefficient for $z(\beta)$")
    ax.set_ylabel(r"coefficient for $z(\gamma)$")
    ax.legend(frameon=False, loc="lower left", fontsize=10, handlelength=1.8)
    save_figure(fig, out_dir, "striatum_hippocampus_pc1_directions")


def main() -> None:
    args = parse_args()
    setup_style()

    results_root = Path(args.results_root)
    out_dir = Path(args.out_dir)
    stats_dir = Path(args.stats_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    datasets = [load_dataset(results_root, args.filter_level, dataset) for dataset in DATASETS]

    for data in datasets:
        fig, ax = plt.subplots(figsize=(6.9, 5.4))
        draw_panel(ax, data, title=data.name.capitalize())
        save_figure(fig, out_dir, f"{data.name}_gene_coefficient_pca")

    all_corrs = []
    for data in datasets:
        all_corrs.extend(merged_plot_table(data)["r"].dropna().to_numpy(dtype=float))
    shared_abs = float(np.nanmax(np.abs(all_corrs)))

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), constrained_layout=True)
    scatter = None
    for ax, data in zip(axes, datasets):
        scatter = draw_panel(
            ax,
            data,
            title=data.name.capitalize(),
            show_ylabel=(data.name == "striatum"),
            add_colorbar=False,
            vmin=-shared_abs,
            vmax=shared_abs,
        )
    assert scatter is not None
    cbar = fig.colorbar(scatter, ax=axes, fraction=0.035, pad=0.02)
    cbar.set_label(r"corr(gene, $\eta$)", rotation=90, labelpad=13, fontsize=16)
    cbar.ax.tick_params(labelsize=12)
    save_figure(fig, out_dir, "striatum_hippocampus_gene_coefficient_pca")

    plot_pc1_comparison(datasets, out_dir)
    write_stats(datasets, stats_dir)

    print(f"Wrote figures to {out_dir}")
    print(f"Wrote stats to {stats_dir}")


if __name__ == "__main__":
    main()
