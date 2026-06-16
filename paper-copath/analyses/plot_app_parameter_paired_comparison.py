#!/usr/bin/env python3
"""Paired APP vs MAPT group-comparison plots for REGION-RF parameters."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PARAMETERS = ["r", "beta", "gamma"]
PROTEINS = [("syn", "Synuclein", "#2563eb"), ("tau", "Tau", "#dc2626")]
RHAT_CUTOFF = 1.05


def parameter_values(df: pd.DataFrame, parameter: str, condition: str) -> pd.Series:
    if parameter == "r":
        return pd.to_numeric(df[f"alpha_{condition}"], errors="coerce") * pd.to_numeric(
            df[f"beta_{condition}"], errors="coerce"
        )
    return pd.to_numeric(df[f"{parameter}_{condition}"], errors="coerce")


def parameter_mask(df: pd.DataFrame, parameter: str, condition: str) -> pd.Series:
    if parameter == "r":
        return (
            pd.to_numeric(df[f"alpha_rhat_{condition}"], errors="coerce") <= RHAT_CUTOFF
        ) & (
            pd.to_numeric(df[f"beta_rhat_{condition}"], errors="coerce") <= RHAT_CUTOFF
        )
    return pd.to_numeric(df[f"{parameter}_rhat_{condition}"], errors="coerce") <= RHAT_CUTOFF


def display_parameter(parameter: str) -> str:
    return r"$r=\alpha\beta$" if parameter == "r" else parameter


def p_text(p: float) -> str:
    if p < 1e-4:
        return "p < 1e-4"
    if p < 1e-3:
        return f"p = {p:.1e}"
    return f"p = {p:.3f}"


def load_pairs(
    comparison_dir: Path,
    require_rhat: bool = True,
    mapt_reference_zscore: bool = False,
    require_active: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    stats_rows = []
    for protein, label, color in PROTEINS:
        df = pd.read_csv(comparison_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
        for parameter in PARAMETERS:
            mask = pd.Series(True, index=df.index)
            if require_active:
                mask &= df["active_any"].astype(bool)
            if require_rhat:
                mask &= parameter_mask(df, parameter, "app")
                mask &= parameter_mask(df, parameter, "mapt")
            sub = df.loc[mask, ["region_index", "region"]].copy()
            sub["app"] = parameter_values(df.loc[mask], parameter, "app").to_numpy()
            sub["mapt"] = parameter_values(df.loc[mask], parameter, "mapt").to_numpy()
            sub = sub[np.isfinite(sub["app"]) & np.isfinite(sub["mapt"])]
            mapt_reference_mean = sub["mapt"].mean()
            mapt_reference_sd = sub["mapt"].std(ddof=1)
            if mapt_reference_zscore:
                sub["app_raw"] = sub["app"]
                sub["mapt_raw"] = sub["mapt"]
                sub["app"] = (sub["app"] - mapt_reference_mean) / mapt_reference_sd
                sub["mapt"] = (sub["mapt"] - mapt_reference_mean) / mapt_reference_sd
            diff = sub["app"] - sub["mapt"]
            ttest = stats.ttest_rel(sub["app"], sub["mapt"])
            mean_diff = diff.mean()
            frac_app_greater = (diff > 0).mean()
            frac_follow_trend = frac_app_greater if mean_diff >= 0 else 1 - frac_app_greater
            stats_rows.append(
                {
                    "protein": protein,
                    "protein_label": label,
                    "parameter": parameter,
                    "n": len(sub),
                    "mapt_mean": sub["mapt"].mean(),
                    "app_mean": sub["app"].mean(),
                    "mean_diff_app_minus_mapt": mean_diff,
                    "mapt_reference_mean": mapt_reference_mean,
                    "mapt_reference_sd": mapt_reference_sd,
                    "frac_app_greater": frac_app_greater,
                    "frac_follow_trend": frac_follow_trend,
                    "paired_t": ttest.statistic,
                    "paired_t_p": ttest.pvalue,
                }
            )
            for _, row in sub.iterrows():
                rows.append(
                    {
                        "protein": protein,
                        "protein_label": label,
                        "parameter": parameter,
                        "color": color,
                        "region_index": row["region_index"],
                        "region": row["region"],
                        "MAPT": row["mapt"],
                        "APP": row["app"],
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(stats_rows)


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_panel(ax, sub: pd.DataFrame, stat: pd.Series, color: str) -> None:
    rng = np.random.default_rng(7)
    mapt_x = 0 + rng.normal(0, 0.018, len(sub))
    app_x = 1 + rng.normal(0, 0.018, len(sub))
    for x0, x1, y0, y1 in zip(mapt_x, app_x, sub["MAPT"], sub["APP"]):
        ax.plot([x0, x1], [y0, y1], color="0.72", lw=0.45, alpha=0.34, zorder=1)
    ax.scatter(mapt_x, sub["MAPT"], s=11, color="0.45", alpha=0.38, linewidth=0, zorder=2)
    ax.scatter(app_x, sub["APP"], s=11, color=color, alpha=0.38, linewidth=0, zorder=2)

    box = ax.boxplot(
        [sub["MAPT"], sub["APP"]],
        positions=[0, 1],
        widths=0.44,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "lw": 1.2},
        boxprops={"lw": 1.0},
        whiskerprops={"lw": 1.0},
        capprops={"lw": 1.0},
    )
    box["boxes"][0].set_facecolor("0.86")
    box["boxes"][0].set_alpha(0.75)
    box["boxes"][1].set_facecolor(color)
    box["boxes"][1].set_alpha(0.25)

    means = [sub["MAPT"].mean(), sub["APP"].mean()]
    ax.scatter([0, 1], means, marker="D", s=38, color=["0.15", color], edgecolor="white", linewidth=0.7, zorder=4)

    ymin = min(sub["MAPT"].min(), sub["APP"].min())
    ymax = max(sub["MAPT"].max(), sub["APP"].max())
    yrange = ymax - ymin if ymax > ymin else 1.0
    bracket_y = ymax + 0.10 * yrange
    tick_y = ymax + 0.04 * yrange
    ax.plot([0, 0, 1, 1], [tick_y, bracket_y, bracket_y, tick_y], color="black", lw=1.0)
    ax.text(0.5, bracket_y + 0.03 * yrange, p_text(stat["paired_t_p"]), ha="center", va="bottom", fontsize=8.5)

    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(ymin - 0.08 * yrange, bracket_y + 0.18 * yrange)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MAPT", "APP"])
    style_axis(ax)


def plot_figure(
    pairs: pd.DataFrame,
    stats_df: pd.DataFrame,
    out_path: Path,
    mapt_reference_zscore: bool = False,
    show_n: bool = True,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(11.8, 7.4), sharex=False)
    for row_idx, (protein, label, color) in enumerate(PROTEINS):
        for col_idx, parameter in enumerate(PARAMETERS):
            ax = axes[row_idx, col_idx]
            sub = pairs[(pairs["protein"] == protein) & (pairs["parameter"] == parameter)]
            stat = stats_df[(stats_df["protein"] == protein) & (stats_df["parameter"] == parameter)].iloc[0]
            diff_label = "Δz" if mapt_reference_zscore else "Δ"
            draw_panel(ax, sub, stat, color)
            if row_idx == 0:
                ax.set_title(display_parameter(parameter))
            if col_idx == 0:
                ax.set_ylabel(label)
            label_text = (
                f"{diff_label}={stat['mean_diff_app_minus_mapt']:.3g}\n"
                f"{100 * stat['frac_follow_trend']:.0f}% trend"
            )
            if show_n:
                label_text = f"n={int(stat['n'])}\n{label_text}"
            ax.text(
                0.03,
                0.985,
                label_text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.5,
                color="0.25",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.70, "pad": 1.2},
            )
    if mapt_reference_zscore:
        fig.suptitle("Paired APP vs MAPT REGION-RF parameter comparison, MAPT-referenced z scale", fontsize=15, y=0.98)
    else:
        fig.suptitle("Paired APP vs MAPT REGION-RF parameter comparison", fontsize=15, y=0.98)
    fig.supxlabel("Condition", y=0.04)
    if mapt_reference_zscore:
        fig.supylabel("MAPT-referenced z score", x=0.01)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--include-one-chain-survivors",
        action="store_true",
        help="Include active regions retained after iterative chain removal even when R-hat is NA because only one chain remains.",
    )
    parser.add_argument(
        "--mapt-reference-zscore",
        action="store_true",
        help="Plot values on a MAPT-referenced z scale: (value - mean(MAPT)) / SD(MAPT), per protein and parameter.",
    )
    parser.add_argument(
        "--hide-n",
        action="store_true",
        help="Do not print sample sizes inside each panel.",
    )
    parser.add_argument(
        "--include-inactive-regions",
        action="store_true",
        help="Include every region in the paired table instead of filtering to active_any regions.",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[2]
    comparison_dir = project_root / "paper-copath" / "results" / "region_rf_condition_comparison"
    figure_dir = project_root / "paper-copath" / "figures" / "region_rf_condition_comparison"
    pairs, stats_df = load_pairs(
        comparison_dir,
        require_rhat=not args.include_one_chain_survivors,
        mapt_reference_zscore=args.mapt_reference_zscore,
        require_active=not args.include_inactive_regions,
    )
    suffix = "_pathological_regions" if args.include_one_chain_survivors else ""
    if args.include_inactive_regions:
        suffix = "_all_412_regions" if args.include_one_chain_survivors else "_all_regions_unfiltered"
    if args.mapt_reference_zscore:
        suffix += "_mapt_zscore"
    stats_df.to_csv(comparison_dir / f"app_parameter_paired_ttests{suffix}.csv", index=False)
    plot_figure(
        pairs,
        stats_df,
        figure_dir / f"app_parameter_paired_comparison{suffix}",
        mapt_reference_zscore=args.mapt_reference_zscore,
        show_n=not args.hide_n,
    )
    print(figure_dir / f"app_parameter_paired_comparison{suffix}")


if __name__ == "__main__":
    main()
