#!/usr/bin/env python3
"""Test whether MAPT gamma strata modify APP-associated gamma shifts."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture


PROTEINS = [("syn", "Synuclein", "#2563eb"), ("tau", "Tau", "#dc2626")]
RHAT_CUTOFF = 1.05
RANDOM_SEED = 13


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


def p_text(p: float) -> str:
    if not np.isfinite(p):
        return "p = NA"
    if p < 1e-4:
        return "p < 1e-4"
    if p < 1e-3:
        return f"p = {p:.1e}"
    return f"p = {p:.3f}"


def load_gamma_pairs(
    comparison_dir: Path,
    protein: str,
    require_rhat: bool = True,
    require_active: bool = True,
) -> pd.DataFrame:
    df = pd.read_csv(comparison_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
    mask = pd.Series(True, index=df.index)
    if require_active:
        mask &= df["active_any"].astype(bool)
    if require_rhat:
        mask &= pd.to_numeric(df["gamma_rhat_app"], errors="coerce") <= RHAT_CUTOFF
        mask &= pd.to_numeric(df["gamma_rhat_mapt"], errors="coerce") <= RHAT_CUTOFF
    out = df.loc[mask, ["region_index", "region", "gamma_mapt", "gamma_app"]].copy()
    out["gamma_mapt"] = pd.to_numeric(out["gamma_mapt"], errors="coerce")
    out["gamma_app"] = pd.to_numeric(out["gamma_app"], errors="coerce")
    out = out.dropna()
    out["gamma_diff_app_minus_mapt"] = out["gamma_app"] - out["gamma_mapt"]
    return out


def assign_mapt_gamma_strata(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = df[["gamma_mapt"]].to_numpy()
    bic_rows = []
    models: dict[int, GaussianMixture] = {}
    for k in [1, 2, 3]:
        gm = GaussianMixture(n_components=k, n_init=100, random_state=RANDOM_SEED)
        gm.fit(x)
        models[k] = gm
        bic_rows.append(
            {
                "components": k,
                "bic": gm.bic(x),
                "aic": gm.aic(x),
                "component_means": ";".join(f"{v:.6g}" for v in sorted(gm.means_.ravel())),
            }
        )

    gm2 = models[2]
    labels = gm2.predict(x)
    means = gm2.means_.ravel()
    order = np.argsort(means)
    label_name = {order[0]: "low MAPT gamma", order[1]: "high MAPT gamma"}
    label_rank = {order[0]: 0, order[1]: 1}

    out = df.copy()
    out["mapt_gamma_stratum"] = [label_name[label] for label in labels]
    out["mapt_gamma_stratum_rank"] = [label_rank[label] for label in labels]
    out["mapt_gamma_stratum_probability"] = gm2.predict_proba(x).max(axis=1)
    out["mapt_gamma_stratum_component_mean"] = [means[label] for label in labels]
    return out, pd.DataFrame(bic_rows)


def summarize_strata(df: pd.DataFrame, protein: str, label: str) -> list[dict[str, float | int | str]]:
    rows = []
    for stratum in ["low MAPT gamma", "high MAPT gamma"]:
        sub = df[df["mapt_gamma_stratum"] == stratum]
        diff = sub["gamma_diff_app_minus_mapt"]
        paired = stats.ttest_rel(sub["gamma_app"], sub["gamma_mapt"])
        n_app_greater = int((diff > 0).sum())
        binom = stats.binomtest(n_app_greater, n=len(sub), p=0.5, alternative="two-sided")
        rows.append(
            {
                "protein": protein,
                "protein_label": label,
                "comparison": "within_stratum_paired_shift",
                "mapt_gamma_stratum": stratum,
                "n": len(sub),
                "gamma_mapt_mean": sub["gamma_mapt"].mean(),
                "gamma_app_mean": sub["gamma_app"].mean(),
                "mean_diff_app_minus_mapt": diff.mean(),
                "median_diff_app_minus_mapt": diff.median(),
                "frac_app_greater": n_app_greater / len(sub),
                "paired_t": paired.statistic,
                "p_value": paired.pvalue,
                "binomial_sign_p": binom.pvalue,
            }
        )

    low = df[df["mapt_gamma_stratum"] == "low MAPT gamma"]["gamma_diff_app_minus_mapt"]
    high = df[df["mapt_gamma_stratum"] == "high MAPT gamma"]["gamma_diff_app_minus_mapt"]
    welch = stats.ttest_ind(high, low, equal_var=False)
    spearman = stats.spearmanr(df["gamma_mapt"], df["gamma_diff_app_minus_mapt"])
    rows.append(
        {
            "protein": protein,
            "protein_label": label,
            "comparison": "high_minus_low_shift_difference",
            "mapt_gamma_stratum": "high - low",
            "n": len(df),
            "gamma_mapt_mean": df["gamma_mapt"].mean(),
            "gamma_app_mean": df["gamma_app"].mean(),
            "mean_diff_app_minus_mapt": high.mean() - low.mean(),
            "median_diff_app_minus_mapt": high.median() - low.median(),
            "frac_app_greater": np.nan,
            "paired_t": welch.statistic,
            "p_value": welch.pvalue,
            "binomial_sign_p": np.nan,
        }
    )
    rows.append(
        {
            "protein": protein,
            "protein_label": label,
            "comparison": "continuous_mapt_gamma_vs_shift_spearman",
            "mapt_gamma_stratum": "continuous",
            "n": len(df),
            "gamma_mapt_mean": df["gamma_mapt"].mean(),
            "gamma_app_mean": df["gamma_app"].mean(),
            "mean_diff_app_minus_mapt": spearman.statistic,
            "median_diff_app_minus_mapt": np.nan,
            "frac_app_greater": np.nan,
            "paired_t": np.nan,
            "p_value": spearman.pvalue,
            "binomial_sign_p": np.nan,
        }
    )
    return rows


def component_density(model: GaussianMixture, xs: np.ndarray) -> list[np.ndarray]:
    densities = []
    component_order = np.argsort(model.means_.ravel())
    for idx in component_order:
        mean = model.means_.ravel()[idx]
        var = model.covariances_.ravel()[idx]
        weight = model.weights_[idx]
        sd = np.sqrt(var)
        densities.append(weight * stats.norm.pdf(xs, mean, sd))
    return densities


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_results(
    region_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    bic_df: pd.DataFrame,
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.2))
    stratum_colors = {"low MAPT gamma": "#0f766e", "high MAPT gamma": "#b45309"}

    for row_idx, (protein, label, color) in enumerate(PROTEINS):
        sub = region_df[region_df["protein"] == protein]
        x = sub[["gamma_mapt"]].to_numpy()
        gm2 = GaussianMixture(n_components=2, n_init=100, random_state=RANDOM_SEED).fit(x)
        xs = np.linspace(sub["gamma_mapt"].min(), sub["gamma_mapt"].max(), 400)

        ax = axes[row_idx, 0]
        ax.hist(sub["gamma_mapt"], bins=28, density=True, color="0.82", edgecolor="white")
        for density, stratum in zip(component_density(gm2, xs), ["low MAPT gamma", "high MAPT gamma"]):
            ax.plot(xs, density, color=stratum_colors[stratum], lw=2.0)
        ax.plot(xs, np.exp(gm2.score_samples(xs.reshape(-1, 1))), color="0.15", lw=1.2, ls="--")
        best_k = bic_df[bic_df["protein"] == protein].sort_values("bic").iloc[0]["components"]
        delta_bic_2_vs_1 = (
            bic_df[(bic_df["protein"] == protein) & (bic_df["components"] == 1)]["bic"].iloc[0]
            - bic_df[(bic_df["protein"] == protein) & (bic_df["components"] == 2)]["bic"].iloc[0]
        )
        ax.set_title(f"{label}: MAPT gamma distribution")
        ax.set_xlabel("MAPT gamma")
        ax.set_ylabel("Density")
        ax.text(
            0.03,
            0.96,
            f"2-state split shown\nbest BIC: {int(best_k)} components\nΔBIC 1-2 = {delta_bic_2_vs_1:.1f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            color="0.25",
        )
        style_axis(ax)

        ax = axes[row_idx, 1]
        for stratum, ss in sub.groupby("mapt_gamma_stratum"):
            ax.scatter(
                ss["gamma_mapt"],
                ss["gamma_diff_app_minus_mapt"],
                s=20,
                alpha=0.72,
                linewidth=0,
                color=stratum_colors[stratum],
                label=stratum.replace(" MAPT gamma", ""),
            )
        slope = stats.linregress(sub["gamma_mapt"], sub["gamma_diff_app_minus_mapt"])
        line_y = slope.intercept + slope.slope * xs
        ax.plot(xs, line_y, color="0.2", lw=1.3)
        ax.axhline(0, color="0.35", lw=1.0, ls=":")
        cont = summary_df[
            (summary_df["protein"] == protein)
            & (summary_df["comparison"] == "continuous_mapt_gamma_vs_shift_spearman")
        ].iloc[0]
        ax.set_title("Continuous trend")
        ax.set_xlabel("MAPT gamma")
        ax.set_ylabel("APP - MAPT gamma")
        ax.text(
            0.03,
            0.96,
            f"Spearman rho={cont['mean_diff_app_minus_mapt']:.2f}\n{p_text(cont['p_value'])}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            color="0.25",
        )
        ax.legend(frameon=False, fontsize=8, loc="lower left")
        style_axis(ax)

        ax = axes[row_idx, 2]
        positions = [0, 1]
        strata = ["low MAPT gamma", "high MAPT gamma"]
        rng = np.random.default_rng(RANDOM_SEED)
        values = [sub[sub["mapt_gamma_stratum"] == stratum]["gamma_diff_app_minus_mapt"] for stratum in strata]
        box = ax.boxplot(
            values,
            positions=positions,
            widths=0.48,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "lw": 1.1},
            boxprops={"lw": 1.0},
            whiskerprops={"lw": 1.0},
            capprops={"lw": 1.0},
        )
        for patch, stratum in zip(box["boxes"], strata):
            patch.set_facecolor(stratum_colors[stratum])
            patch.set_alpha(0.22)
        for pos, stratum in zip(positions, strata):
            ss = sub[sub["mapt_gamma_stratum"] == stratum]
            jitter = rng.normal(0, 0.055, len(ss))
            ax.scatter(
                pos + jitter,
                ss["gamma_diff_app_minus_mapt"],
                s=16,
                alpha=0.52,
                linewidth=0,
                color=stratum_colors[stratum],
            )
            stat = summary_df[
                (summary_df["protein"] == protein)
                & (summary_df["comparison"] == "within_stratum_paired_shift")
                & (summary_df["mapt_gamma_stratum"] == stratum)
            ].iloc[0]
            trend = stat["frac_app_greater"] if stat["mean_diff_app_minus_mapt"] >= 0 else 1 - stat["frac_app_greater"]
            ax.text(
                0.04 + 0.50 * pos,
                0.96,
                f"n={int(stat['n'])}\nΔ={stat['mean_diff_app_minus_mapt']:.3g}\n{100 * trend:.0f}% trend\n{p_text(stat['p_value'])}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.2,
                color="0.25",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
            )
        interaction = summary_df[
            (summary_df["protein"] == protein)
            & (summary_df["comparison"] == "high_minus_low_shift_difference")
        ].iloc[0]
        ax.axhline(0, color="0.35", lw=1.0, ls=":")
        ax.set_title(f"Shift by MAPT gamma state\nhigh-low shift {p_text(interaction['p_value'])}")
        ax.set_xticks(positions)
        ax.set_xticklabels(["Low", "High"])
        ax.set_ylabel("APP - MAPT gamma")
        style_axis(ax)

    fig.suptitle("Does baseline MAPT gamma state modify APP-associated gamma shifts?", fontsize=15, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
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
        "--include-inactive-regions",
        action="store_true",
        help="Include every region in the paired table instead of filtering to active_any regions.",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[2]
    comparison_dir = project_root / "paper-copath" / "results" / "region_rf_condition_comparison"
    figure_dir = project_root / "paper-copath" / "figures" / "region_rf_condition_comparison"
    suffix = "_pathological_regions" if args.include_one_chain_survivors else ""
    if args.include_inactive_regions:
        suffix = "_all_412_regions" if args.include_one_chain_survivors else "_all_regions_unfiltered"

    all_regions = []
    all_bics = []
    all_summary = []
    for protein, label, _ in PROTEINS:
        pairs = load_gamma_pairs(
            comparison_dir,
            protein,
            require_rhat=not args.include_one_chain_survivors,
            require_active=not args.include_inactive_regions,
        )
        pairs, bic = assign_mapt_gamma_strata(pairs)
        pairs.insert(0, "protein_label", label)
        pairs.insert(0, "protein", protein)
        bic.insert(0, "protein_label", label)
        bic.insert(0, "protein", protein)
        all_regions.append(pairs)
        all_bics.append(bic)
        all_summary.extend(summarize_strata(pairs, protein, label))

    region_df = pd.concat(all_regions, ignore_index=True)
    bic_df = pd.concat(all_bics, ignore_index=True)
    summary_df = pd.DataFrame(all_summary)
    summary_df["p_value_fdr"] = bh_fdr(summary_df["p_value"])
    summary_df["binomial_sign_p_fdr"] = bh_fdr(summary_df["binomial_sign_p"])

    region_df.to_csv(comparison_dir / f"gamma_mapt_stratified_app_effect_regions{suffix}.csv", index=False)
    bic_df.to_csv(comparison_dir / f"gamma_mapt_strata_gmm_bic{suffix}.csv", index=False)
    summary_df.to_csv(comparison_dir / f"gamma_mapt_stratified_app_effect_summary{suffix}.csv", index=False)
    plot_results(region_df, summary_df, bic_df, figure_dir / f"gamma_mapt_stratified_app_effect{suffix}")
    print(figure_dir / f"gamma_mapt_stratified_app_effect{suffix}")


if __name__ == "__main__":
    main()
