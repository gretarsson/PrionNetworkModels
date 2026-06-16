#!/usr/bin/env python3
"""Combined convergence summary for co-pathology REGION-RF fits."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATASETS = [
    ("syn_app", "Syn APP", "#2563eb"),
    ("syn_mapt", "Syn MAPT", "#60a5fa"),
    ("tau_app", "Tau APP", "#dc2626"),
    ("tau_mapt", "Tau MAPT", "#f97316"),
]
MAIN_PARAMETERS = ["alpha", "beta", "gamma"]
ALL_PARAMETERS = ["alpha", "beta", "gamma", "u0", "sigma"]
RHAT_CUTOFF = 1.05


def load_dataset(dataset: str, project_root: Path, adjusted_dir: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = project_root / "runs" / "region_rf" / f"copath_{dataset}"
    adjusted_path = None if adjusted_dir is None else adjusted_dir / f"{dataset}_region_rf_posterior_summary_long.csv"
    posterior_path = adjusted_path if adjusted_path is not None and adjusted_path.exists() else root / "region_rf_posterior_summary_long.csv"
    posterior = pd.read_csv(posterior_path)
    summary = pd.read_csv(root / "region_rf_summary.csv")
    active = summary[["region_index", "observed_peak_mean"]].copy()
    active["active"] = pd.to_numeric(active["observed_peak_mean"], errors="coerce").fillna(0) > 0
    posterior = posterior.merge(active[["region_index", "active"]], on="region_index", how="left")
    posterior["active"] = posterior["active"].fillna(False)
    return posterior, summary


def region_max_rhat(posterior: pd.DataFrame, parameters: list[str]) -> pd.DataFrame:
    sub = posterior[posterior["parameter"].isin(parameters)].copy()
    region = (
        sub.groupby(["region_index", "region"], as_index=False)
        .agg(max_rhat=("rhat", "max"), active=("active", "max"))
        .sort_values("region_index")
    )
    return region


def summarize_dataset(dataset: str, label: str, posterior: pd.DataFrame) -> dict[str, float | int | str]:
    rows: dict[str, float | int | str] = {"dataset": dataset, "label": label}
    for key, parameters in [("main", MAIN_PARAMETERS), ("all", ALL_PARAMETERS)]:
        sub = posterior[posterior["parameter"].isin(parameters)]
        region = region_max_rhat(posterior, parameters)
        active_region = region[region["active"]]
        rows[f"{key}_n_parameter_rows"] = len(sub)
        rows[f"{key}_max_rhat"] = float(sub["rhat"].max())
        rows[f"{key}_mean_rhat"] = float(sub["rhat"].mean())
        rows[f"{key}_parameter_rows_gt_1_01"] = int((sub["rhat"] > 1.01).sum())
        rows[f"{key}_parameter_rows_gt_1_05"] = int((sub["rhat"] > RHAT_CUTOFF).sum())
        rows[f"{key}_regions_total"] = len(region)
        rows[f"{key}_regions_pass"] = int((region["max_rhat"] <= RHAT_CUTOFF).sum())
        rows[f"{key}_active_regions_total"] = len(active_region)
        rows[f"{key}_active_regions_pass"] = int((active_region["max_rhat"] <= RHAT_CUTOFF).sum())
    for parameter in MAIN_PARAMETERS:
        p = posterior[posterior["parameter"] == parameter]
        rows[f"{parameter}_rows_gt_1_05"] = int((p["rhat"] > RHAT_CUTOFF).sum())
        rows[f"{parameter}_max_rhat"] = float(p["rhat"].max())
    return rows


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.sort(np.asarray(values, dtype=float))
    y = np.arange(1, len(values) + 1) / len(values)
    return values, y


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_figure(project_root: Path, out_dir: Path, adjusted_dir: Path | None = None, title_suffix: str = "") -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    loaded = []
    summary_rows = []
    parameter_rows = []
    for dataset, label, color in DATASETS:
        posterior, _ = load_dataset(dataset, project_root, adjusted_dir)
        loaded.append((dataset, label, color, posterior, region_max_rhat(posterior, MAIN_PARAMETERS)))
        summary_rows.append(summarize_dataset(dataset, label, posterior))
        for parameter in MAIN_PARAMETERS:
            p = posterior[posterior["parameter"] == parameter]
            parameter_rows.append(
                {
                    "dataset": dataset,
                    "label": label,
                    "parameter": parameter,
                    "pass_fraction": float((p["rhat"] <= RHAT_CUTOFF).mean()),
                    "max_rhat": float(p["rhat"].max()),
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "region_rf_convergence_summary.csv", index=False)
    pd.DataFrame(parameter_rows).to_csv(out_dir / "region_rf_parameter_convergence_summary.csv", index=False)

    fig = plt.figure(figsize=(12.8, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.48, wspace=0.42)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_ecdf = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[1, 0])
    ax_rank = fig.add_subplot(gs[1, 1])

    labels = [label for _, label, _ in DATASETS]
    x = np.arange(len(labels))
    all_pass = summary["main_regions_pass"] / summary["main_regions_total"]
    active_pass = summary["main_active_regions_pass"] / summary["main_active_regions_total"]
    width = 0.34
    ax_bar.bar(x - width / 2, all_pass, width, color="#94a3b8", label="all regions")
    ax_bar.bar(x + width / 2, active_pass, width, color="#0f766e", label="pathology-active")
    ax_bar.axhline(1.0, color="0.75", lw=0.8)
    ax_bar.set_ylim(0.80, 1.08)
    ax_bar.set_ylabel(f"fraction with max R-hat <= {RHAT_CUTOFF}")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels, rotation=25, ha="right")
    ax_bar.set_title("Regional convergence")
    ax_bar.legend(frameon=False, loc="lower left")
    for i, row in summary.iterrows():
        ax_bar.text(
            i - width / 2,
            min(all_pass.iloc[i] + 0.008, 1.035),
            f"{int(row['main_regions_pass'])}/{int(row['main_regions_total'])}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )
        ax_bar.text(
            i + width / 2,
            min(active_pass.iloc[i] + 0.008, 1.035),
            f"{int(row['main_active_regions_pass'])}/{int(row['main_active_regions_total'])}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )
    style_axis(ax_bar)

    for dataset, label, color, _, region in loaded:
        vals, y = ecdf(region["max_rhat"].to_numpy())
        ax_ecdf.plot(vals, y, color=color, lw=2.0, label=label)
    ax_ecdf.axvline(RHAT_CUTOFF, color="black", ls="--", lw=1.1)
    ax_ecdf.set_xlim(0.995, 1.25)
    ax_ecdf.set_ylim(0, 1.01)
    ax_ecdf.set_xlabel("regional max R-hat across alpha, beta, gamma")
    ax_ecdf.set_ylabel("cumulative fraction of regions")
    ax_ecdf.set_title("Distribution of regional R-hat")
    ax_ecdf.legend(frameon=False, loc="lower right")
    style_axis(ax_ecdf)

    heat = pd.DataFrame(parameter_rows).pivot(index="label", columns="parameter", values="pass_fraction")
    heat = heat.loc[labels, MAIN_PARAMETERS]
    im = ax_heat.imshow(heat.to_numpy(), vmin=0.85, vmax=1.0, cmap="viridis", aspect="auto")
    ax_heat.set_xticks(np.arange(len(MAIN_PARAMETERS)))
    ax_heat.set_xticklabels(MAIN_PARAMETERS)
    ax_heat.set_yticks(np.arange(len(labels)))
    ax_heat.set_yticklabels(labels)
    ax_heat.set_title(f"Parameter-row pass fraction, R-hat <= {RHAT_CUTOFF}")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            ax_heat.text(j, i, f"{100 * heat.iloc[i, j]:.1f}%", ha="center", va="center", fontsize=9, color="white")
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.040, pad=0.04)
    cbar.set_label("pass fraction")

    for dataset, label, color, _, region in loaded:
        rank = np.arange(1, len(region) + 1)
        ordered = region.sort_values("max_rhat", ascending=False).reset_index(drop=True)
        rank = np.arange(1, len(ordered) + 1)
        ax_rank.scatter(rank, ordered["max_rhat"], s=9, alpha=0.65, color=color, label=label, linewidth=0)
    ax_rank.axhline(RHAT_CUTOFF, color="black", ls="--", lw=1.1)
    ax_rank.set_yscale("log")
    ax_rank.set_xlabel("regions ordered by max R-hat")
    ax_rank.set_ylabel("regional max R-hat", labelpad=10)
    ax_rank.set_title("Tail of non-converged regional fits")
    ax_rank.legend(frameon=False, loc="upper right", ncol=2, fontsize=8)
    style_axis(ax_rank)

    title = "Co-pathology REGION-RF convergence diagnostics"
    if title_suffix:
        title = f"{title} ({title_suffix})"
    fig.suptitle(title, y=0.98, fontsize=15)
    fig.savefig(out_dir / "region_rf_convergence_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "region_rf_convergence_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adjusted-dir", default=None, help="Optional directory with <dataset>_region_rf_posterior_summary_long.csv files.")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--title-suffix", default="")
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[2]
    out_dir = Path(args.out_dir) if args.out_dir else project_root / "paper-copath" / "figures" / "region_rf_convergence"
    adjusted_dir = Path(args.adjusted_dir) if args.adjusted_dir else None
    make_figure(project_root, out_dir, adjusted_dir, args.title_suffix)
    print(out_dir)


if __name__ == "__main__":
    main()
