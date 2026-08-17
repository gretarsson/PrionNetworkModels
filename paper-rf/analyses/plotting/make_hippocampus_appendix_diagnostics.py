#!/usr/bin/env python3
"""Build standalone hippocampus diagnostic panels for the appendix."""

from pathlib import Path
import csv
import json
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
ACCEPTED_RUN = ROOT / "runs" / "hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C4"
ALL_CHAIN_RUN = ROOT / "runs" / "hippocampus_DIFF-RF_RETRO_striatum-global-priors_C1_C2_C3_C4"
OUT = ROOT / "paper-rf" / "figures" / "hippocampus_appendix"
BLUE = "#3f5fa8"


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#e6e6e6", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=11, width=0.8, length=3)


def seed_label_map():
    metadata = json.loads((ACCEPTED_RUN / "metadata.json").read_text())
    seed_indices = metadata["spec"]["seeding"]["seed_indices"]
    network_path = ROOT / metadata["spec"]["data"]["network"]
    with network_path.open() as f:
        labels = next(csv.reader(f))[1:]
    out = {}
    for seed_number, region_index in enumerate(seed_indices, start=1):
        out[f"seed_values[{seed_number}]"] = f"seed {labels[region_index - 1]}"
    return out


def pretty_parameter(name, seeds):
    if name == "rho":
        return r"$\rho$"
    if name == "alpha":
        return r"$\alpha$"
    if name == "sigma":
        return r"$\sigma$"
    if name in seeds:
        return seeds[name]
    return name


def plot_rhat_panel(rhat: pd.DataFrame, family: str, out_name: str, title: str):
    sub = rhat.loc[rhat["family"].eq(family)].copy()
    seeds = seed_label_map()
    if family in {"beta", "gamma"}:
        sub["rank"] = sub["parameter"].str.extract(r"\[(\d+)\]").astype(int)
        sub = sub.sort_values("rank")
        x = np.arange(1, len(sub) + 1)
        xlabel = "Parameter index"
        width = 5.0
    else:
        order = ["rho", "alpha", "seed_values[1]", "seed_values[2]", "seed_values[3]", "sigma"]
        sub["order"] = sub["parameter"].map({name: i for i, name in enumerate(order)})
        sub = sub.sort_values("order")
        x = np.arange(1, len(sub) + 1)
        xlabel = ""
        width = 5.4

    fig, ax = plt.subplots(figsize=(width, 3.6))
    ax.scatter(
        x,
        sub["rhat"],
        s=50 if family == "global" else 24,
        color=BLUE,
        alpha=0.82,
        edgecolor="#262626",
        linewidth=0.45,
        zorder=3,
    )
    for val, col, lw in ((1.00, "#8c8c8c", 1.0), (1.01, "#238b45", 1.2), (1.05, "#fdae61", 1.4), (1.10, "#ef3b2c", 1.4)):
        ax.axhline(val, color=col, linewidth=lw, linestyle=(0, (6, 4)), zorder=1)
    ax.set_ylabel(r"$\hat{R}$", fontsize=14, labelpad=9)
    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_title(title, fontsize=14, pad=9)
    ax.set_ylim(0.998, 1.105)

    if family == "global":
        ax.set_xticks(x)
        ax.set_xticklabels([pretty_parameter(v, seeds) for v in sub["parameter"]], rotation=35, ha="right")
        ax.set_xlim(0.45, len(sub) + 0.55)
    else:
        ax.set_xlim(0, len(sub) + 1)

    style_axes(ax)
    fig.subplots_adjust(left=0.16, right=0.99, top=0.86, bottom=0.28 if family == "global" else 0.17)
    fig.savefig(OUT / out_name, bbox_inches="tight")
    plt.close(fig)


def plot_loglik():
    metrics = pd.read_csv(ALL_CHAIN_RUN / "plots" / "diagnostics" / "chain_fit_metrics.csv")
    retained = metrics["chain"].isin([3, 4])
    colors = np.where(retained, "#2f6f9f", "#9a9a9a")

    fig, ax = plt.subplots(figsize=(5.1, 3.6))
    ax.scatter(metrics["chain"], metrics["loglik_all"], s=72, color=colors, edgecolor="black", linewidth=0.5, zorder=3)
    ax.plot(metrics["chain"], metrics["loglik_all"], color="#c8c8c8", linewidth=1.0, zorder=1)
    ax.set_xlabel("Source chain", fontsize=12)
    ax.set_ylabel("Log-likelihood", fontsize=12)
    ax.set_xticks(metrics["chain"])
    ax.margins(x=0.22)

    yspan = metrics["loglik_all"].max() - metrics["loglik_all"].min()
    yspan = yspan if yspan > 0 else 1
    for _, row in metrics.iterrows():
        label = "retained" if row["chain"] in [3, 4] else "excluded"
        ax.text(
            row["chain"] + 0.16,
            row["loglik_all"] + 0.035 * yspan,
            label,
            fontsize=10,
            color="#2f6f9f" if label == "retained" else "#707070",
            ha="left",
            va="bottom",
        )

    style_axes(ax)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.96, bottom=0.18)
    fig.savefig(OUT / "log_likelihood_all_observations.pdf", bbox_inches="tight")
    plt.close(fig)


def copy_fit_panels():
    shutil.copy2(
        ACCEPTED_RUN / "plots" / "predicted_vs_observed.pdf",
        OUT / "predicted_vs_observed.pdf",
    )
    shutil.copy2(
        ACCEPTED_RUN / "plots" / "retrodiction" / "top_pathology_panels" / "top_observed_pathology_1_to_4.pdf",
        OUT / "top_4_regions.pdf",
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    copy_fit_panels()
    rhat = pd.read_csv(ACCEPTED_RUN / "plots" / "diagnostics" / "rhat_summary.csv")
    plot_rhat_panel(rhat, "global", "global_parameters.pdf", "Global parameters")
    plot_rhat_panel(rhat, "beta", "beta_parameters.pdf", r"$\beta$ parameters")
    plot_rhat_panel(rhat, "gamma", "gamma_parameters.pdf", r"$\gamma$ parameters")
    plot_loglik()


if __name__ == "__main__":
    main()
