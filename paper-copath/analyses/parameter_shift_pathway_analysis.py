#!/usr/bin/env python3
"""Direct gene/pathway associations with APP-MAPT parameter shifts."""

from __future__ import annotations

import argparse
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PARAMETERS = ["r", "beta", "gamma"]
PARAMETER_LABELS = {"r": r"$\Delta r$", "beta": r"$\Delta \beta$", "gamma": r"$\Delta \gamma$"}
PROTEIN_LABELS = {"syn": "Synuclein", "tau": "Tau"}
PROTEIN_COLORS = {"syn": "#2563eb", "tau": "#dc2626"}


def strip_hemi(region: str) -> str:
    region = str(region).strip()
    return region[1:] if len(region) > 1 and region[0] in {"i", "c"} else region


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


def load_expression(path: Path) -> tuple[pd.DataFrame, list[str]]:
    expr = pd.read_csv(path)
    expr = expr.rename(columns={expr.columns[0]: "region_base"})
    gene_cols = [c for c in expr.columns if c != "region_base"]
    return expr, gene_cols


def load_parameter_table(comparison_dir: Path, protein: str) -> pd.DataFrame:
    df = pd.read_csv(comparison_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
    df["region_base"] = df["region"].map(strip_hemi)
    df["r_app"] = pd.to_numeric(df["alpha_app"], errors="coerce") * pd.to_numeric(df["beta_app"], errors="coerce")
    df["r_mapt"] = pd.to_numeric(df["alpha_mapt"], errors="coerce") * pd.to_numeric(df["beta_mapt"], errors="coerce")
    for parameter in PARAMETERS:
        diff_col = f"{parameter}_diff_app_minus_mapt"
        if diff_col not in df.columns:
            df[diff_col] = pd.to_numeric(df[f"{parameter}_app"], errors="coerce") - pd.to_numeric(
                df[f"{parameter}_mapt"], errors="coerce"
            )
    return df[df["active_any"].astype(bool)].copy()


def gene_correlations(merged: pd.DataFrame, gene_cols: list[str], outcome: str) -> pd.DataFrame:
    y = pd.to_numeric(merged[outcome], errors="coerce").to_numpy(dtype=float)
    rows = []
    for gene in gene_cols:
        x = pd.to_numeric(merged[gene], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 8:
            continue
        pr = stats.pearsonr(x[mask], y[mask])
        sr = stats.spearmanr(x[mask], y[mask])
        rows.append(
            {
                "gene": gene,
                "n_regions": int(mask.sum()),
                "pearson_r": float(pr.statistic),
                "pearson_p": float(pr.pvalue),
                "spearman_r": float(sr.statistic),
                "spearman_p": float(sr.pvalue),
            }
        )
    out = pd.DataFrame(rows).sort_values("pearson_r", ascending=False)
    out["pearson_p_fdr"] = bh_fdr(out["pearson_p"].to_numpy())
    out["spearman_p_fdr"] = bh_fdr(out["spearman_p"].to_numpy())
    return out


def run_gsea(ranked: pd.DataFrame, gmt_path: Path, out_dir: Path, seed: int) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    rnk_path = out_dir / "ranked_genes.rnk"
    ranked[["gene", "pearson_r"]].dropna().drop_duplicates("gene").to_csv(
        rnk_path, sep="\t", header=False, index=False
    )
    pre_res = gp.prerank(
        rnk=str(rnk_path),
        gene_sets=str(gmt_path),
        min_size=10,
        max_size=500,
        permutation_num=1000,
        outdir=str(out_dir / "gseapy"),
        seed=seed,
        verbose=False,
    )
    return pre_res.res2d.copy()


def clean_term(term: str) -> str:
    term = str(term)
    term = term.replace("KEGG_2019_Mouse__", "")
    term = term.replace("_", " ")
    return term


def truncate(text: str, max_len: int = 42) -> str:
    return text if len(text) <= max_len else f"{text[: max_len - 1]}..."


def select_terms(gsea: pd.DataFrame, max_terms: int = 14) -> list[str]:
    sig = gsea.copy()
    sig["abs_nes"] = sig["NES"].abs()
    sig = sig.sort_values(["FDR q-val", "abs_nes"], ascending=[True, False])
    terms = []
    for term in sig["Term"]:
        if term not in terms:
            terms.append(term)
        if len(terms) >= max_terms:
            break
    return terms


def plot_panel(gsea: pd.DataFrame, figure_dir: Path) -> None:
    terms = select_terms(gsea)
    plot_df = gsea[gsea["Term"].isin(terms)].copy()
    plot_df["term_clean"] = plot_df["Term"].map(clean_term)
    plot_df["outcome_label"] = plot_df["parameter"].map(PARAMETER_LABELS)
    y_labels = [truncate(clean_term(t)) for t in terms[::-1]]
    term_to_y = {term: i for i, term in enumerate(terms[::-1])}

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), sharey=True)
    norm = plt.Normalize(-3.0, 3.0)
    cmap = plt.get_cmap("coolwarm")
    for ax, protein in zip(axes, ["syn", "tau"]):
        sub = plot_df[plot_df["protein"] == protein].copy()
        for _, row in sub.iterrows():
            x = PARAMETERS.index(row["parameter"])
            y = term_to_y[row["Term"]]
            fdr = float(row["FDR q-val"])
            size = 28 + 210 * np.clip(-np.log10(max(fdr, 1e-6)) / 6, 0, 1)
            ax.scatter(
                x,
                y,
                s=size,
                c=[cmap(norm(float(row["NES"])))],
                edgecolor="white",
                linewidth=0.55,
            )
        ax.set_title(PROTEIN_LABELS[protein])
        ax.set_xticks(np.arange(len(PARAMETERS)))
        ax.set_xticklabels([PARAMETER_LABELS[p] for p in PARAMETERS])
        ax.set_xlim(-0.5, len(PARAMETERS) - 0.5)
        ax.set_yticks(np.arange(len(y_labels)))
        ax.set_yticklabels(y_labels)
        ax.grid(axis="y", color="0.90", lw=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(length=0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes, pad=0.02, fraction=0.04)
    cbar.set_label("GSEA NES")
    for fdr, x in [(0.25, 0.80), (0.05, 0.90), (0.01, 1.00)]:
        axes[1].scatter([], [], s=28 + 210 * np.clip(-np.log10(fdr) / 6, 0, 1), color="0.5", label=f"FDR {fdr:g}")
    axes[1].legend(frameon=False, title="size", loc="lower right", bbox_to_anchor=(1.42, 0.02), fontsize=8)
    fig.suptitle("Pathways associated with APP-MAPT parameter shifts", y=0.98)
    fig.text(0.5, 0.02, "Pathology-active regions; one-chain retained fits included; direct parameter correlations, no PCA", ha="center", fontsize=9, color="0.35")
    fig.subplots_adjust(left=0.23, right=0.86, bottom=0.13, top=0.88, wspace=0.08)
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_path = figure_dir / "parameter_shift_pathway_gsea_panel"
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_top10_by_analysis(gsea: pd.DataFrame, figure_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(22.0, 10.6))
    norm = plt.Normalize(-4.2, 4.2)
    cmap = plt.get_cmap("coolwarm")
    q_min = 1e-3
    q_max_for_size = 0.15

    def q_to_size(q: np.ndarray | float) -> np.ndarray | float:
        q = np.asarray(q, dtype=float)
        score = -np.log10(np.clip(q, q_min, q_max_for_size))
        lo = -np.log10(q_max_for_size)
        hi = -np.log10(q_min)
        scaled = np.clip((score - lo) / (hi - lo), 0, 1)
        return 35 + 240 * scaled

    for row_idx, protein in enumerate(["syn", "tau"]):
        for col_idx, parameter in enumerate(PARAMETERS):
            ax = axes[row_idx, col_idx]
            sub = gsea[(gsea["protein"] == protein) & (gsea["parameter"] == parameter)].copy()
            sub["abs_nes"] = sub["NES"].abs()
            sub = sub.sort_values(["FDR q-val", "abs_nes"], ascending=[True, False]).head(10)
            sub = sub.iloc[::-1].reset_index(drop=True)

            y = np.arange(len(sub))
            fdr = pd.to_numeric(sub["FDR q-val"], errors="coerce").to_numpy(dtype=float)
            sizes = q_to_size(fdr)
            ax.scatter(
                sub["NES"],
                y,
                s=sizes,
                c=sub["NES"],
                cmap=cmap,
                norm=norm,
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
            ax.axvline(0, color="0.55", lw=0.9)
            ax.grid(axis="y", color="0.92", lw=0.7)
            ax.set_xlim(-4.55, 5.05)
            ax.set_yticks(y)
            ax.set_yticklabels([truncate(clean_term(t), 38) for t in sub["Term"]], fontsize=8.2)
            ax.set_xlabel("GSEA NES")
            ax.set_title(f"{PROTEIN_LABELS[protein]} {PARAMETER_LABELS[parameter]}")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="y", length=0)

    fig.subplots_adjust(left=0.15, right=0.83, bottom=0.09, top=0.92, wspace=1.05, hspace=0.38)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = fig.add_axes([0.84, 0.20, 0.018, 0.62])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("GSEA NES")
    legend_q = [0.15, 0.05, 0.01, 0.001]
    handles = [
        axes[1, 2].scatter([], [], s=q_to_size(q), color="0.45", edgecolor="white", linewidth=0.6)
        for q in legend_q
    ]
    labels = ["q=0.15", "q=0.05", "q=0.01", "q<=0.001"]
    fig.legend(
        handles,
        labels,
        title="FDR q-value",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(0.875, 0.50),
        fontsize=9,
        title_fontsize=10,
    )
    fig.suptitle("Top pathways associated with APP-MAPT parameter shifts", y=0.985)
    fig.text(
        0.5,
        0.025,
        "Top 10 pathways selected separately within each protein and parameter. Pathology-active regions; one-chain retained fits included; no PCA.",
        ha="center",
        fontsize=9,
        color="0.35",
    )
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_path = figure_dir / "parameter_shift_pathway_gsea_top10_by_analysis"
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-dir", default="paper-copath/results/region_rf_condition_comparison")
    parser.add_argument("--expression", default="paper-rf/data/transcriptomics/avg_Pangea_exp.csv")
    parser.add_argument("--gene-sets", default="paper-rf/results/figure_6_7/all/striatum/enrichment/gseapy/gene_sets.gmt")
    parser.add_argument("--out-dir", default="paper-copath/results/parameter_shift_pathways")
    parser.add_argument("--figure-dir", default="paper-copath/figures/parameter_shift_pathways")
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    comparison_dir = Path(args.comparison_dir)
    expression, gene_cols = load_expression(Path(args.expression))
    out_dir = Path(args.out_dir)
    figure_dir = Path(args.figure_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_gsea = []
    for protein in ["syn", "tau"]:
        params = load_parameter_table(comparison_dir, protein)
        merged = params.merge(expression, on="region_base", how="inner")
        for parameter in PARAMETERS:
            outcome = f"{parameter}_diff_app_minus_mapt"
            corr = gene_correlations(merged, gene_cols, outcome)
            param_dir = out_dir / protein / parameter
            param_dir.mkdir(parents=True, exist_ok=True)
            corr.to_csv(param_dir / "gene_parameter_shift_correlations.csv", index=False)
            gsea = run_gsea(corr, Path(args.gene_sets), param_dir, args.seed)
            gsea["protein"] = protein
            gsea["parameter"] = parameter
            gsea["outcome"] = outcome
            gsea.to_csv(param_dir / "gsea_parameter_shift_all.csv", index=False)
            all_gsea.append(gsea)

    gsea_all = pd.concat(all_gsea, ignore_index=True)
    gsea_all.to_csv(out_dir / "gsea_parameter_shift_all_proteins.csv", index=False)
    plot_panel(gsea_all, figure_dir)
    plot_top10_by_analysis(gsea_all, figure_dir)
    print(out_dir)
    print(figure_dir)


if __name__ == "__main__":
    main()
