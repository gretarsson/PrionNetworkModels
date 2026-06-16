#!/usr/bin/env python3
"""Simple cell-type comparison panel for APP-MAPT parameter shifts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUTCOMES = ["r_diff_app_minus_mapt", "beta_diff_app_minus_mapt", "gamma_diff_app_minus_mapt"]
OUTCOME_LABELS = {
    "r_diff_app_minus_mapt": r"$\Delta r$",
    "beta_diff_app_minus_mapt": r"$\Delta \beta$",
    "gamma_diff_app_minus_mapt": r"$\Delta \gamma$",
}
PROTEIN_LABELS = {"syn": "Synuclein", "tau": "Tau"}
CELL_ORDER = ["Dopa", "Nora", "Sero", "Hist", "Chol", "GABA", "GABA-Glyc", "Glut", "Glut-GABA", "Unknown"]


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    result_dir = project_root / "paper-copath" / "results" / "amyloid_sensitivity"
    figure_dir = project_root / "paper-copath" / "figures" / "amyloid_sensitivity"

    corr = pd.read_csv(result_dir / "celltype_delta_axis_correlations.csv")
    corr = corr[corr["outcome"].isin(OUTCOMES)].copy()
    corr["cell_type"] = pd.Categorical(corr["cell_type"], categories=CELL_ORDER, ordered=True)
    corr = corr.dropna(subset=["cell_type"])

    summary = corr[
        ["protein", "outcome", "cell_type", "n", "spearman_r", "spearman_p", "p_perm", "p_fdr"]
    ].sort_values(["protein", "cell_type", "outcome"])
    summary.to_csv(result_dir / "celltype_parameter_shift_panel_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.6), sharey=True)
    image = None
    for ax, protein in zip(axes, ["syn", "tau"]):
        sub = corr[corr["protein"] == protein]
        mat = (
            sub.pivot(index="cell_type", columns="outcome", values="spearman_r")
            .reindex(index=CELL_ORDER, columns=OUTCOMES)
            .to_numpy(dtype=float)
        )
        q = (
            sub.pivot(index="cell_type", columns="outcome", values="p_fdr")
            .reindex(index=CELL_ORDER, columns=OUTCOMES)
            .to_numpy(dtype=float)
        )
        image = ax.imshow(mat, cmap="coolwarm", vmin=-0.35, vmax=0.35, aspect="auto")
        ax.set_title(PROTEIN_LABELS[protein])
        ax.set_xticks(np.arange(len(OUTCOMES)))
        ax.set_xticklabels([OUTCOME_LABELS[o] for o in OUTCOMES])
        ax.set_yticks(np.arange(len(CELL_ORDER)))
        ax.set_yticklabels(CELL_ORDER)
        ax.tick_params(length=0)
        for y in range(mat.shape[0]):
            for x in range(mat.shape[1]):
                if not np.isfinite(mat[y, x]):
                    continue
                text = f"{mat[y, x]:.2f}"
                color = "white" if abs(mat[y, x]) > 0.22 else "0.18"
                ax.text(x, y, text, ha="center", va="center", fontsize=8.2, color=color)
                if np.isfinite(q[y, x]) and q[y, x] < 0.10:
                    marker = "*" if q[y, x] < 0.05 else "."
                    ax.text(x + 0.33, y - 0.28, marker, ha="center", va="center", fontsize=10, color="black")
        style_axis(ax)

    fig.subplots_adjust(left=0.10, right=0.86, top=0.82, bottom=0.18, wspace=0.08)
    cax = fig.add_axes([0.89, 0.22, 0.018, 0.54])
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_label("Spearman rho")
    fig.suptitle("Cell-type composition vs APP-MAPT parameter shifts", y=0.96)
    fig.text(
        0.48,
        0.06,
        "Pathology-active regions; one-chain retained fits included. * permutation FDR < 0.05, . FDR < 0.10",
        ha="center",
        fontsize=9,
        color="0.35",
    )

    figure_dir.mkdir(parents=True, exist_ok=True)
    out_path = figure_dir / "celltype_parameter_shift_panel"
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
