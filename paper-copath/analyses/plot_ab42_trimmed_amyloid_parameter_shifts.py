#!/usr/bin/env python3
"""AB42 sensitivity panel after trimming extreme amyloid-load values."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PARAMETERS = ["r", "beta", "gamma"]
PARAMETER_LABELS = {"r": r"$r=\alpha\beta$", "beta": "beta", "gamma": "gamma"}
PROTEIN_LABELS = {"syn": "Synuclein", "tau": "Tau"}
PROTEIN_COLORS = {"syn": "#2563eb", "tau": "#dc2626"}


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def bh_fdr(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full_like(p, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[finite])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.clip(ranked, 0, 1)
    return out


def correlation_stats(x: np.ndarray, y: np.ndarray) -> dict[str, float | int]:
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 3:
        return {
            "n": int(finite.sum()),
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
        }
    pr = stats.pearsonr(x[finite], y[finite])
    sr = stats.spearmanr(x[finite], y[finite])
    return {
        "n": int(finite.sum()),
        "pearson_r": float(pr.statistic),
        "pearson_p": float(pr.pvalue),
        "spearman_r": float(sr.statistic),
        "spearman_p": float(sr.pvalue),
    }


def p_text(p: float) -> str:
    if not np.isfinite(p):
        return "p=NA"
    if p < 1e-4:
        return "p<1e-4"
    if p < 1e-3:
        return f"p={p:.1e}"
    return f"p={p:.2g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", default="paper-copath/results/amyloid_sensitivity")
    parser.add_argument("--figure-dir", default="paper-copath/figures/collaborator_update")
    parser.add_argument("--trim-quantile", type=float, default=0.025)
    args = parser.parse_args()

    result_dir = Path(args.result_dir)
    figure_dir = Path(args.figure_dir)
    amy_col = "ab42_treatment_mean_prelimval"
    trim_q = args.trim_quantile

    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.6), sharex=False)
    for row_idx, protein in enumerate(["syn", "tau"]):
        df = pd.read_csv(result_dir / f"{protein}_amyloid_sensitivity_region_table.csv")
        df = df[df["active_rhat"].astype(bool)].copy()
        x_all = pd.to_numeric(df[amy_col], errors="coerce")
        finite_x = x_all[np.isfinite(x_all)]
        low = float(finite_x.quantile(trim_q))
        high = float(finite_x.quantile(1 - trim_q))
        df["ab42_trimmed_included"] = np.isfinite(x_all) & (x_all >= low) & (x_all <= high)

        for col_idx, parameter in enumerate(PARAMETERS):
            ax = axes[row_idx, col_idx]
            y_col = f"{parameter}_diff_app_minus_mapt"
            sub = df[df["ab42_trimmed_included"]].copy()
            x = pd.to_numeric(sub[amy_col], errors="coerce").to_numpy(dtype=float)
            y = pd.to_numeric(sub[y_col], errors="coerce").to_numpy(dtype=float)
            st = correlation_stats(x, y)
            rows.append(
                {
                    "protein": protein,
                    "parameter": parameter,
                    "amyloid": "ab42",
                    "trim_quantile_each_tail": trim_q,
                    "ab42_low_cutoff": low,
                    "ab42_high_cutoff": high,
                    "n_before_trim_finite_ab42": int(finite_x.shape[0]),
                    "n_after_trim": st["n"],
                    **st,
                }
            )

            ax.scatter(
                x,
                y,
                s=28,
                alpha=0.72,
                color=PROTEIN_COLORS[protein],
                edgecolor="white",
                linewidth=0.35,
            )
            finite = np.isfinite(x) & np.isfinite(y)
            if finite.sum() >= 3:
                slope, intercept, *_ = stats.linregress(x[finite], y[finite])
                xx = np.linspace(np.nanmin(x[finite]), np.nanmax(x[finite]), 100)
                ax.plot(xx, intercept + slope * xx, color="black", lw=1.8)
            ax.axhline(0, color="0.75", lw=0.8)
            ax.text(
                0.04,
                0.96,
                f"r={st['pearson_r']:.2f}, {p_text(st['pearson_p'])}\n"
                rf"$\rho$={st['spearman_r']:.2f}, {p_text(st['spearman_p'])}"
                f"\nn={st['n']}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.8,
            )
            ax.set_title(f"{PROTEIN_LABELS[protein]}: delta {PARAMETER_LABELS[parameter]}")
            ax.set_xlabel("AB42 amyloid load")
            ax.set_ylabel("APP - MAPT")
            style_axis(ax)

    stats_df = pd.DataFrame(rows)
    stats_df["pearson_p_fdr"] = bh_fdr(stats_df["pearson_p"].to_numpy())
    stats_df["spearman_p_fdr"] = bh_fdr(stats_df["spearman_p"].to_numpy())
    result_dir.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(result_dir / "ab42_parameter_shift_correlations_trimmed.csv", index=False)

    pct = 100 * trim_q
    fig.suptitle(f"AB42 load vs APP-MAPT parameter shifts after trimming outer {pct:g}% tails", y=0.98)
    fig.text(
        0.5,
        0.02,
        "Pathology-active regions; one-chain retained fits included; AB42 trimmed within each protein.",
        ha="center",
        fontsize=9,
        color="0.35",
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.09, top=0.90, wspace=0.34, hspace=0.52)

    figure_dir.mkdir(parents=True, exist_ok=True)
    out_path = figure_dir / "07b_ab42_amyloid_vs_parameter_shifts_trimmed"
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
