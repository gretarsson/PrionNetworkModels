#!/usr/bin/env python3
"""Standalone APP effect plot for alpha * beta REGION-RF composite."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PROTEINS = [("syn", "Synuclein", "#2563eb"), ("tau", "Tau", "#dc2626")]
RHAT_CUTOFF = 1.05


def p_text(p: float) -> str:
    if p < 1e-4:
        return "p < 1e-4"
    if p < 1e-3:
        return f"p = {p:.1e}"
    return f"p = {p:.3f}"


def load_product_data(comparison_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    stats_rows = []
    for protein, label, color in PROTEINS:
        df = pd.read_csv(comparison_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
        mask = df["active_any"].astype(bool)
        for parameter in ["alpha", "beta"]:
            mask &= pd.to_numeric(df[f"{parameter}_rhat_app"], errors="coerce") <= RHAT_CUTOFF
            mask &= pd.to_numeric(df[f"{parameter}_rhat_mapt"], errors="coerce") <= RHAT_CUTOFF
        sub = df.loc[mask, ["region_index", "region", "alpha_app", "beta_app", "alpha_mapt", "beta_mapt"]].copy()
        for col in ["alpha_app", "beta_app", "alpha_mapt", "beta_mapt"]:
            sub[col] = pd.to_numeric(sub[col], errors="coerce")
        sub = sub.dropna()
        sub["app_alpha_beta"] = sub["alpha_app"] * sub["beta_app"]
        sub["mapt_alpha_beta"] = sub["alpha_mapt"] * sub["beta_mapt"]
        sub["diff"] = sub["app_alpha_beta"] - sub["mapt_alpha_beta"]
        sub = sub[np.isfinite(sub["app_alpha_beta"]) & np.isfinite(sub["mapt_alpha_beta"])]
        ttest = stats.ttest_rel(sub["app_alpha_beta"], sub["mapt_alpha_beta"])
        n = len(sub)
        mean_diff = sub["diff"].mean()
        sd_diff = sub["diff"].std(ddof=1)
        se = sd_diff / np.sqrt(n)
        ci = stats.t.ppf(0.975, n - 1) * se
        stats_rows.append(
            {
                "protein": protein,
                "protein_label": label,
                "n": n,
                "mapt_mean_alpha_beta": sub["mapt_alpha_beta"].mean(),
                "app_mean_alpha_beta": sub["app_alpha_beta"].mean(),
                "mean_diff_app_minus_mapt": mean_diff,
                "median_diff_app_minus_mapt": sub["diff"].median(),
                "sd_diff": sd_diff,
                "ci_low": mean_diff - ci,
                "ci_high": mean_diff + ci,
                "standardized_mean_diff": mean_diff / sd_diff if sd_diff > 0 else np.nan,
                "frac_app_greater": (sub["diff"] > 0).mean(),
                "paired_t": ttest.statistic,
                "paired_t_p": ttest.pvalue,
            }
        )
        for _, row in sub.iterrows():
            rows.append(
                {
                    "protein": protein,
                    "protein_label": label,
                    "color": color,
                    "region_index": row["region_index"],
                    "region": row["region"],
                    "MAPT": row["mapt_alpha_beta"],
                    "APP": row["app_alpha_beta"],
                    "diff": row["diff"],
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(stats_rows)


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_paired_panel(ax, sub: pd.DataFrame, stat: pd.Series, color: str) -> None:
    rng = np.random.default_rng(11)
    x0 = rng.normal(0, 0.018, len(sub))
    x1 = 1 + rng.normal(0, 0.018, len(sub))
    for a, b, y0, y1 in zip(x0, x1, sub["MAPT"], sub["APP"]):
        ax.plot([a, b], [y0, y1], color="0.72", lw=0.45, alpha=0.34, zorder=1)
    ax.scatter(x0, sub["MAPT"], s=12, color="0.45", alpha=0.4, linewidth=0, zorder=2)
    ax.scatter(x1, sub["APP"], s=12, color=color, alpha=0.4, linewidth=0, zorder=2)
    box = ax.boxplot(
        [sub["MAPT"], sub["APP"]],
        positions=[0, 1],
        widths=0.44,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "lw": 1.2},
    )
    box["boxes"][0].set_facecolor("0.86")
    box["boxes"][0].set_alpha(0.75)
    box["boxes"][1].set_facecolor(color)
    box["boxes"][1].set_alpha(0.25)
    ax.scatter(
        [0, 1],
        [sub["MAPT"].mean(), sub["APP"].mean()],
        marker="D",
        s=40,
        color=["0.15", color],
        edgecolor="white",
        linewidth=0.7,
        zorder=4,
    )
    ymin = min(sub["MAPT"].min(), sub["APP"].min())
    ymax = max(sub["MAPT"].max(), sub["APP"].max())
    yrange = ymax - ymin if ymax > ymin else 1.0
    bracket_y = ymax + 0.12 * yrange
    tick_y = ymax + 0.05 * yrange
    ax.plot([0, 0, 1, 1], [tick_y, bracket_y, bracket_y, tick_y], color="black", lw=1.0)
    ax.text(0.5, bracket_y + 0.035 * yrange, p_text(stat["paired_t_p"]), ha="center", va="bottom", fontsize=9)
    ax.text(
        0.03,
        0.94,
        f"n={int(stat['n'])}\nΔ={stat['mean_diff_app_minus_mapt']:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="0.25",
    )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MAPT", "APP"])
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(ymin - 0.08 * yrange, bracket_y + 0.20 * yrange)
    ax.set_ylabel(r"$\alpha \times \beta$")
    style_axis(ax)


def draw_difference_panel(ax, sub: pd.DataFrame, stat: pd.Series, color: str, y: float) -> None:
    values = sub["diff"].to_numpy()
    parts = ax.violinplot([values], positions=[y], vert=False, widths=0.62, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_alpha(0.25)
        body.set_edgecolor("none")
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    ax.plot([q1, q3], [y, y], color=color, lw=5.5, solid_capstyle="round")
    ax.scatter([med], [y], color="white", edgecolor=color, s=42, zorder=4, linewidth=1.1)
    ax.errorbar(
        stat["mean_diff_app_minus_mapt"],
        y + 0.22,
        xerr=[
            [stat["mean_diff_app_minus_mapt"] - stat["ci_low"]],
            [stat["ci_high"] - stat["mean_diff_app_minus_mapt"]],
        ],
        fmt="o",
        color=color,
        capsize=3,
        ms=5,
        lw=1.5,
    )
    ax.text(
        0.98,
        y,
        f"{100 * stat['frac_app_greater']:.0f}% >0",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9,
        color="0.25",
    )


def plot_figure(pairs: pd.DataFrame, stats_df: pd.DataFrame, out_path: Path) -> None:
    fig = plt.figure(figsize=(11.8, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.85], hspace=0.45, wspace=0.30)
    for col, (protein, label, color) in enumerate(PROTEINS):
        ax = fig.add_subplot(gs[0, col])
        sub = pairs[pairs["protein"] == protein]
        stat = stats_df[stats_df["protein"] == protein].iloc[0]
        draw_paired_panel(ax, sub, stat, color)
        ax.set_title(label)

    ax_diff = fig.add_subplot(gs[1, 0])
    for y, (protein, label, color) in enumerate(PROTEINS):
        sub = pairs[pairs["protein"] == protein]
        stat = stats_df[stats_df["protein"] == protein].iloc[0]
        draw_difference_panel(ax_diff, sub, stat, color, y)
    ax_diff.axvline(0, color="black", ls="--", lw=1.0)
    ax_diff.set_yticks([0, 1])
    ax_diff.set_yticklabels([label for _, label, _ in PROTEINS])
    ax_diff.set_xlabel(r"APP - MAPT difference in $\alpha \times \beta$")
    ax_diff.set_title("Signed regional shift")
    style_axis(ax_diff)

    ax_eff = fig.add_subplot(gs[1, 1])
    y = np.arange(len(PROTEINS))
    colors = [color for _, _, color in PROTEINS]
    ax_eff.axvline(0, color="black", ls="--", lw=1.0)
    for yi, (_, stat), color in zip(y, stats_df.iterrows(), colors):
        dz = stat["standardized_mean_diff"]
        se_dz = 1 / np.sqrt(stat["n"])
        ax_eff.errorbar(dz, yi, xerr=1.96 * se_dz, fmt="o", color=color, capsize=3, ms=7, lw=1.7)
        p = stat["paired_t_p"]
        ax_eff.text(
            dz + 0.045,
            yi,
            p_text(p),
            va="center",
            ha="left",
            fontsize=9,
            color="0.25",
        )
    ax_eff.set_yticks(y)
    ax_eff.set_yticklabels([label for _, label, _ in PROTEINS])
    ax_eff.invert_yaxis()
    ax_eff.set_xlabel(r"standardized paired shift, mean(diff) / SD(diff)")
    ax_eff.set_title("Effect size")
    style_axis(ax_eff)

    fig.suptitle(r"APP effect on REGION-RF growth drive ($\alpha \times \beta$)", fontsize=15, y=0.98)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    comparison_dir = project_root / "paper-copath" / "results" / "region_rf_condition_comparison"
    figure_dir = project_root / "paper-copath" / "figures" / "region_rf_condition_comparison"
    pairs, stats_df = load_product_data(comparison_dir)
    stats_df.to_csv(comparison_dir / "app_alpha_beta_product_summary.csv", index=False)
    plot_figure(pairs, stats_df, figure_dir / "app_alpha_beta_product")
    print(figure_dir / "app_alpha_beta_product")


if __name__ == "__main__":
    main()
