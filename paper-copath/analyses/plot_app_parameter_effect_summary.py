#!/usr/bin/env python3
"""Summarize APP effects on REGION-RF parameters as signed shifts."""

from __future__ import annotations

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


def load_shift_rows(comparison_dir: Path) -> pd.DataFrame:
    rows = []
    for protein, label, color in PROTEINS:
        df = pd.read_csv(comparison_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
        for parameter in PARAMETERS:
            mask = df["active_any"].astype(bool)
            mask &= parameter_mask(df, parameter, "app")
            mask &= parameter_mask(df, parameter, "mapt")
            app = parameter_values(df.loc[mask], parameter, "app")
            mapt = parameter_values(df.loc[mask], parameter, "mapt")
            diff = app - mapt
            diff = diff[np.isfinite(diff)]
            n = len(diff)
            mean = float(diff.mean())
            sd = float(diff.std(ddof=1))
            se = sd / np.sqrt(n)
            ci = stats.t.ppf(0.975, n - 1) * se
            ttest = stats.ttest_1samp(diff, popmean=0.0)
            rows.append(
                {
                    "protein": protein,
                    "protein_label": label,
                    "parameter": parameter,
                    "color": color,
                    "n": n,
                    "mean_diff": mean,
                    "median_diff": float(diff.median()),
                    "sd_diff": sd,
                    "ci_low": mean - ci,
                    "ci_high": mean + ci,
                    "standardized_mean_diff": mean / sd if sd > 0 else np.nan,
                    "frac_app_greater": float((diff > 0).mean()),
                    "paired_t": float(ttest.statistic),
                    "paired_t_p": float(ttest.pvalue),
                    "values": diff.to_numpy(),
                }
            )
    return pd.DataFrame(rows)


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_shift_distributions(summary: pd.DataFrame, out_path: Path) -> None:
    fig = plt.figure(figsize=(14.4, 8.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.05, 1.0], hspace=0.42, wspace=0.40)

    for col, parameter in enumerate(PARAMETERS):
        ax = fig.add_subplot(gs[0, col])
        sub = summary[summary["parameter"] == parameter]
        positions = np.arange(len(PROTEINS))
        for pos, (_, row) in zip(positions, sub.iterrows()):
            values = row["values"]
            parts = ax.violinplot(
                [values],
                positions=[pos],
                vert=False,
                widths=0.72,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body in parts["bodies"]:
                body.set_facecolor(row["color"])
                body.set_alpha(0.28)
                body.set_edgecolor("none")
            q1, med, q3 = np.percentile(values, [25, 50, 75])
            ax.plot([q1, q3], [pos, pos], color=row["color"], lw=6, alpha=0.85, solid_capstyle="round")
            ax.scatter([med], [pos], color="white", edgecolor=row["color"], s=46, zorder=4, linewidth=1.2)
            ax.errorbar(
                row["mean_diff"],
                pos + 0.24,
                xerr=[[row["mean_diff"] - row["ci_low"]], [row["ci_high"] - row["mean_diff"]]],
                fmt="o",
                color=row["color"],
                capsize=3,
                ms=5,
                lw=1.5,
            )
            ax.text(
                0.98,
                pos,
                f"{100 * row['frac_app_greater']:.0f}% >0",
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=9,
                color="0.25",
            )
        ax.axvline(0, color="black", ls="--", lw=1.1)
        ax.set_yticks(positions)
        ax.set_yticklabels([label for _, label, _ in PROTEINS])
        ax.set_xlabel("APP - MAPT")
        ax.set_title(display_parameter(parameter))
        style_axis(ax)

    ax = fig.add_subplot(gs[1, :3])
    y = np.arange(len(summary))
    labels = [f"{row.protein_label} {display_parameter(row.parameter)}" for row in summary.itertuples()]
    short_labels = [
        f"{'Syn' if row.protein == 'syn' else 'Tau'} {display_parameter(row.parameter)}"
        for row in summary.itertuples()
    ]
    colors = summary["color"].to_list()
    ax.axvline(0, color="black", ls="--", lw=1.1)
    for yi, (_, row) in zip(y, summary.iterrows()):
        dz = row["standardized_mean_diff"]
        n = row["n"]
        se_dz = 1 / np.sqrt(n)
        ax.errorbar(
            dz,
            yi,
            xerr=1.96 * se_dz,
            fmt="o",
            color=row["color"],
            ecolor=row["color"],
            capsize=3,
            ms=7,
            lw=1.7,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("standardized paired shift, mean(APP - MAPT) / SD(diff)")
    ax.set_title("Direction and relative size of APP effect")
    ax.set_xlim(-0.72, 0.68)
    style_axis(ax)

    ax_t = fig.add_subplot(gs[:, 3])
    signed_logp = []
    for _, row in summary.iterrows():
        p = max(row["paired_t_p"], np.nextafter(0, 1))
        signed_logp.append(np.sign(row["mean_diff"]) * -np.log10(p))
    y_t = np.arange(len(summary))
    ax_t.axvline(0, color="black", lw=1.0)
    ax_t.axvline(-np.log10(0.05), color="0.55", ls="--", lw=0.9)
    ax_t.axvline(np.log10(0.05), color="0.55", ls="--", lw=0.9)
    ax_t.barh(y_t, signed_logp, color=colors, alpha=0.85)
    for yi, (_, row), val in zip(y_t, summary.iterrows(), signed_logp):
        p = row["paired_t_p"]
        trend_frac = row["frac_app_greater"] if row["mean_diff"] >= 0 else 1 - row["frac_app_greater"]
        p_text = "p<1e-4" if p < 1e-4 else f"p={p:.2g}"
        label = f"{p_text}\n{100 * trend_frac:.0f}% trend"
        ax_t.text(
            val + (0.25 if val >= 0 else -0.25),
            yi,
            label,
            ha="left" if val >= 0 else "right",
            va="center",
            fontsize=8.5,
            color="0.25",
        )
    ax_t.set_yticks(y_t)
    ax_t.set_yticklabels(short_labels)
    ax_t.invert_yaxis()
    ax_t.set_xlabel("signed -log10 paired t-test p")
    ax_t.set_title("Paired t-test\nAPP - MAPT vs 0")
    style_axis(ax_t)

    fig.suptitle("APP effect on REGION-RF parameters", fontsize=15, y=0.98)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    comparison_dir = project_root / "paper-copath" / "results" / "region_rf_condition_comparison"
    figure_dir = project_root / "paper-copath" / "figures" / "region_rf_condition_comparison"
    summary = load_shift_rows(comparison_dir)
    table = summary.drop(columns=["values", "color"])
    table.to_csv(comparison_dir / "app_parameter_effect_summary.csv", index=False)
    plot_shift_distributions(summary, figure_dir / "app_parameter_effect_summary")
    print(figure_dir / "app_parameter_effect_summary")


if __name__ == "__main__":
    main()
