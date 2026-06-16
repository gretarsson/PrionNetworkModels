#!/usr/bin/env python3
"""Positive MAPT baseline r/beta values vs APP-MAPT parameter shifts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PARAMETERS = ["r", "beta"]
PARAMETER_LABELS = {"r": r"$r=\alpha\beta$", "beta": "beta"}
PROTEINS = [("syn", "Synuclein", "#2563eb"), ("tau", "Tau", "#dc2626")]


def parameter_values(df: pd.DataFrame, parameter: str, condition: str) -> pd.Series:
    if parameter == "r":
        return pd.to_numeric(df[f"alpha_{condition}"], errors="coerce") * pd.to_numeric(
            df[f"beta_{condition}"], errors="coerce"
        )
    return pd.to_numeric(df[f"{parameter}_{condition}"], errors="coerce")


def p_text(p: float) -> str:
    if not np.isfinite(p):
        return "p=NA"
    if p < 1e-4:
        return "p<1e-4"
    if p < 1e-3:
        return f"p={p:.1e}"
    return f"p={p:.3f}"


def bh_fdr(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    out = np.full_like(p, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return pd.Series(out, index=values.index)
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[finite])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.clip(ranked, 0, 1)
    return pd.Series(out, index=values.index)


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def load_data(comparison_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    stats_rows = []
    for protein, protein_label, color in PROTEINS:
        df = pd.read_csv(comparison_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
        df = df[df["active_any"].astype(bool)].copy()
        for parameter in PARAMETERS:
            mapt = parameter_values(df, parameter, "mapt")
            app = parameter_values(df, parameter, "app")
            sub = df[["region_index", "region"]].copy()
            sub["mapt"] = mapt
            sub["app"] = app
            sub["shift"] = app - mapt
            sub = sub[np.isfinite(sub["mapt"]) & np.isfinite(sub["shift"]) & (sub["mapt"] > 0)].copy()
            pearson = stats.pearsonr(sub["mapt"], sub["shift"])
            spearman = stats.spearmanr(sub["mapt"], sub["shift"])
            slope = stats.linregress(sub["mapt"], sub["shift"])
            stats_rows.append(
                {
                    "protein": protein,
                    "protein_label": protein_label,
                    "parameter": parameter,
                    "n_positive_mapt": len(sub),
                    "pearson_r": pearson.statistic,
                    "pearson_p": pearson.pvalue,
                    "spearman_rho": spearman.statistic,
                    "spearman_p": spearman.pvalue,
                    "slope": slope.slope,
                    "slope_p": slope.pvalue,
                }
            )
            sub["protein"] = protein
            sub["protein_label"] = protein_label
            sub["parameter"] = parameter
            sub["color"] = color
            rows.append(sub)
    stats_df = pd.DataFrame(stats_rows)
    stats_df["spearman_p_fdr"] = bh_fdr(stats_df["spearman_p"])
    return pd.concat(rows, ignore_index=True), stats_df


def plot_panel(rows: pd.DataFrame, stats_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.2), sharex=False, sharey=False)
    for row_idx, (protein, protein_label, color) in enumerate(PROTEINS):
        for col_idx, parameter in enumerate(PARAMETERS):
            ax = axes[row_idx, col_idx]
            sub = rows[(rows["protein"] == protein) & (rows["parameter"] == parameter)]
            stat = stats_df[(stats_df["protein"] == protein) & (stats_df["parameter"] == parameter)].iloc[0]
            x = sub["mapt"].to_numpy(dtype=float)
            y = sub["shift"].to_numpy(dtype=float)
            ax.scatter(x, y, s=22, color=color, alpha=0.68, linewidth=0)
            ax.axhline(0, color="0.55", lw=0.9, ls=":")
            fit = stats.linregress(x, y)
            xx = np.linspace(np.nanmin(x), np.nanmax(x), 150)
            ax.plot(xx, fit.intercept + fit.slope * xx, color="0.15", lw=1.4)
            ax.text(
                0.04,
                0.96,
                rf"$\rho$={stat['spearman_rho']:.2f}"
                + f"\n{p_text(stat['spearman_p'])}\nn={int(stat['n_positive_mapt'])}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.7,
                color="0.25",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.70, "pad": 1.3},
            )
            if row_idx == 0:
                ax.set_title(PARAMETER_LABELS[parameter])
            if col_idx == 0:
                ax.set_ylabel(f"{protein_label}\nAPP - MAPT")
            ax.set_xlabel(f"MAPT {PARAMETER_LABELS[parameter]} > 0")
            style_axis(ax)
    fig.suptitle("Positive MAPT baseline values vs APP-associated r/beta shifts", y=0.98)
    fig.text(
        0.5,
        0.03,
        "Pathology-active regions; one-chain retained fits included; regions filtered to positive MAPT baseline value.",
        ha="center",
        fontsize=9,
        color="0.35",
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    comparison_dir = project_root / "paper-copath" / "results" / "region_rf_condition_comparison"
    result_dir = project_root / "paper-copath" / "results" / "region_rf_condition_comparison"
    figure_dir = project_root / "paper-copath" / "figures" / "collaborator_update"
    rows, stats_df = load_data(comparison_dir)
    rows.to_csv(result_dir / "positive_mapt_baseline_vs_shift_regions.csv", index=False)
    stats_df.to_csv(result_dir / "positive_mapt_baseline_vs_shift_summary.csv", index=False)
    out_path = figure_dir / "11_positive_mapt_baseline_vs_r_beta_shift"
    plot_panel(rows, stats_df, out_path)
    print(out_path)


if __name__ == "__main__":
    main()
