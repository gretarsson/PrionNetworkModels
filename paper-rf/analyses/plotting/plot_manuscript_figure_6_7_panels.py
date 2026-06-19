#!/usr/bin/env python3
"""Create manuscript-style Figure 6/7 panels from clean figure_6_7 outputs."""

from __future__ import annotations

import argparse
import math
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact


FILTERS = ["all", "beta_positive"]
DATASETS = [
    ("striatum", "Striatum", "#2b6cb0"),
    ("hippocampus_C1_C4", "Hippocampus C1/C4", "#c2410c"),
]
CELL_ORDER = ["frac_Dopa", "frac_Glut", "frac_GABA-Glyc", "frac_Sero", "frac_Nora", "frac_Hist", "frac_Chol", "frac_GABA", "frac_Glut-GABA", "frac_Unknown"]
CATEGORY_ORDER = ["Metabolism", "Protein homeostasis", "Synapse", "Neurodeg. disease", "Other"]
CATEGORY_COLORS = {
    "Metabolism": "#2ca24d",
    "Protein homeostasis": "#ff7f0e",
    "Synapse": "#8c61aa",
    "Neurodeg. disease": "#c94f54",
    "Other": "#858585",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#dddddd",
            "grid.linewidth": 0.45,
            "grid.alpha": 0.55,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def p_text(p: float) -> str:
    if not np.isfinite(p):
        return "p = NA"
    if p < 1e-4:
        return "p < 1e-4"
    return f"p = {p:.3g}"


def q_to_size(q: float) -> float:
    if not np.isfinite(q):
        return 18.0
    score = -math.log10(max(float(q), 1e-6))
    return float(np.clip(10 + 12 * score, 18, 82))


def category_for_term(term: str) -> str:
    s = str(term).lower()
    if any(k in s for k in ["parkinson", "alzheimer", "huntington", "prion", "amyotrophic", "als", "neurodegenerative"]):
        return "Neurodeg. disease"
    if any(k in s for k in ["synapse", "synaptic", "vesicle", "neurotransmitter", "long-term potentiation", "long-term depression"]):
        return "Synapse"
    if any(k in s for k in ["proteasome", "ubiquitin", "autophagy", "lysosome", "endoplasmic reticulum", "protein", "ribosome", "translation", "folding", "chaperone"]):
        return "Protein homeostasis"
    if any(k in s for k in ["metabolism", "oxidative phosphorylation", "mitochond", "glycolysis", "tca", "citrate cycle", "fatty acid", "lipid", "biosynthesis", "pyruvate", "thermogenesis"]):
        return "Metabolism"
    return "Other"


def load_gsea(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t").copy()
    fdr_col = "FDR q-val" if "FDR q-val" in df.columns else "FDR"
    df["FDR"] = pd.to_numeric(df[fdr_col], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df["TermName"] = df["Term"].astype(str).str.replace("_", " ", regex=False)
    df["Category"] = df["TermName"].map(category_for_term)
    return df.dropna(subset=["NES", "FDR"])


def load_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


def pca_segments(coefs: pd.DataFrame, pca: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = coefs[["coef_beta", "coef_gamma"]].mean(axis=0).to_numpy(dtype=float)
    pc1 = pca.loc[pca["component"] == "PC1", ["loading_beta", "loading_gamma"]].iloc[0].to_numpy(dtype=float)
    pc2 = pca.loc[pca["component"] == "PC2", ["loading_beta", "loading_gamma"]].iloc[0].to_numpy(dtype=float)
    if pc1[1] < 0:
        pc1 = -pc1
    if np.linalg.det(np.vstack([pc1, pc2])) < 0:
        pc2 = -pc2
    return center, pc1, pc2


def segment(center: np.ndarray, vec: np.ndarray, x: np.ndarray, y: np.ndarray, scale: float) -> tuple[np.ndarray, np.ndarray]:
    xrange = float(np.nanmax(x) - np.nanmin(x))
    yrange = float(np.nanmax(y) - np.nanmin(y))
    length = scale * min(
        xrange / max(abs(vec[0]), np.finfo(float).eps),
        yrange / max(abs(vec[1]), np.finfo(float).eps),
    )
    return center - length * vec, center + length * vec


def plot_panel_a(trans_dir: Path, out_dir: Path) -> None:
    coefs = pd.read_csv(trans_dir / "gene_parameter_coefficients.csv")
    corr = pd.read_csv(trans_dir / "gene_eta_correlations.csv")[["gene", "r"]]
    pca = pd.read_csv(trans_dir / "pca_summary.csv")
    pc1_row = pca.loc[pca["component"] == "PC1"].iloc[0]
    center, pc1, pc2 = pca_segments(coefs, pca)
    df = coefs.merge(corr, on="gene", how="left")
    x = df["coef_beta"].to_numpy(dtype=float)
    y = df["coef_gamma"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(5.3, 4.4))
    max_abs = float(np.nanmax(np.abs(df["r"])))
    sc = ax.scatter(x, y, c=df["r"], cmap="coolwarm", vmin=-max_abs, vmax=max_abs, s=25, alpha=0.58, edgecolors="#454545", linewidths=0.22, rasterized=True)
    p2a, p2b = segment(center, pc2, x, y, 0.23)
    p1a, p1b = segment(center, pc1, x, y, 0.36)
    ax.plot([p2a[0], p2b[0]], [p2a[1], p2b[1]], color="#8c8c8c", lw=5, solid_capstyle="round")
    ax.plot([p1a[0], p1b[0]], [p1a[1], p1b[1]], color="black", lw=5, solid_capstyle="round")
    ax.set_xlabel(r"gene coefficient for $z(\beta)$")
    ax.set_ylabel(r"gene coefficient for $z(\gamma)$")
    ax.set_title(rf"PC1 defines $\eta = {pc1[0]:.2f}z(\beta) + {pc1[1]:.2f}z(\gamma)$")
    ax.legend(
        handles=[
            Line2D([0], [0], color="#8c8c8c", lw=4, solid_capstyle="round", label="PC2"),
            Line2D([0], [0], color="black", lw=4, solid_capstyle="round", label="PC1"),
        ],
        frameon=False,
        loc="lower left",
        handlelength=1.7,
        handletextpad=0.25,
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.025)
    cbar.set_label(r"corr(gene, $\eta$)")
    ax.text(0.98, 0.97, f"PC1 {100 * pc1_row['explained_variance_ratio']:.1f}%", transform=ax.transAxes, ha="right", va="top", fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 3})
    save(fig, out_dir / "panel_A_gene_coefficient_pca")


def plot_panel_b(trans_dir: Path, out_dir: Path) -> None:
    df = pd.read_csv(trans_dir / "region_axis.csv")
    pca = pd.read_csv(trans_dir / "pca_summary.csv").iloc[0]
    x = df["beta"].to_numpy(dtype=float)
    y = df["gamma"].to_numpy(dtype=float)
    gx = np.linspace(np.nanmin(x), np.nanmax(x), 160)
    gy = np.linspace(np.nanmin(y), np.nanmax(y), 160)
    xx, yy = np.meshgrid(gx, gy)
    z_beta = (xx - np.nanmean(x)) / np.nanstd(x)
    z_gamma = (yy - np.nanmean(y)) / np.nanstd(y)
    eta_grid = pca["loading_beta"] * z_beta + pca["loading_gamma"] * z_gamma
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.contourf(xx, yy, eta_grid, levels=18, cmap="coolwarm", alpha=0.78)
    ax.contour(xx, yy, eta_grid, levels=6, colors="0.35", alpha=0.45, linewidths=1.5)
    ax.scatter(x, y, s=78, facecolor="white", edgecolor="black", linewidth=1.5, alpha=0.92)
    ax.set_xlabel(r"rise parameter $\beta$")
    ax.set_ylabel(r"fall parameter $\gamma$")
    cbar = fig.colorbar(im, ax=ax, pad=0.025)
    cbar.set_label(r"vulnerability axis $\eta$")
    save(fig, out_dir / "panel_B_beta_gamma_eta_plane")


def plot_panel_d(enrich_dir: Path, out_dir: Path) -> None:
    df = load_gsea(enrich_dir / "gsea_results_all.tsv")
    ranked = pd.read_csv(enrich_dir / "ranked_genes.rnk", sep="\t", header=None, names=["gene", "score"])
    ranked["score"] = pd.to_numeric(ranked["score"], errors="coerce")
    ranked = ranked.dropna(subset=["gene", "score"])
    ranked["gene_key"] = ranked["gene"].astype(str).str.upper()
    score_by_gene = ranked.drop_duplicates("gene_key").set_index("gene_key")["score"]
    gene_sets = {term: {gene.upper() for gene in genes} for term, genes in load_gmt(enrich_dir / "gseapy" / "gene_sets.gmt").items()}
    sig = df[df["FDR"] <= 0.05].copy()
    if sig.empty:
        return
    sig["absNES"] = sig["NES"].abs()
    top = sig.sort_values("absNES", ascending=False).head(10).iloc[::-1].copy()
    fig, ax = plt.subplots(figsize=(6.6, 4.9))
    scores = ranked["score"].to_numpy(dtype=float)
    violins = ax.violinplot(
        [scores for _ in range(len(top))],
        positions=np.arange(len(top)),
        orientation="horizontal",
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in violins["bodies"]:
        body.set_facecolor("#b5b5b5")
        body.set_edgecolor("none")
        body.set_alpha(0.42)
    rng = np.random.default_rng(7)
    for i, row in enumerate(top.itertuples(index=False)):
        color = CATEGORY_COLORS[row.Category]
        member_scores = score_by_gene.reindex(sorted(gene_sets.get(row.Term, set()))).dropna().to_numpy(dtype=float)
        if len(member_scores):
            jitter = rng.uniform(-0.12, 0.12, size=len(member_scores))
            ax.scatter(
                member_scores,
                i + jitter,
                s=q_to_size(float(row.FDR)),
                color=color,
                alpha=0.86,
                edgecolors="none",
                zorder=3,
                rasterized=True,
            )
        ax.text(1.04, i, f"{row.NES:+.2f}", va="center", fontsize=11, transform=ax.get_yaxis_transform())
    ax.axvline(0, color="0.25", lw=0.9)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(["\n".join(textwrap.wrap(t, width=34)) for t in top["TermName"]])
    ax.set_xlabel("r")
    lim = float(np.nanquantile(np.abs(scores), 0.995))
    ax.set_xlim(-lim * 1.08, lim * 1.08)
    ax.text(1.04, len(top) - 0.2, "NES", ha="left", va="bottom", fontsize=12, transform=ax.get_yaxis_transform())
    ax.text(1.15, 0.88, "FDR", transform=ax.transAxes, ha="left", va="center", fontsize=10, clip_on=False)
    for j, q in enumerate([0.05, 0.01, 0.001]):
        y_pos = 0.82 - 0.055 * j
        ax.scatter(
            [1.15],
            [y_pos],
            s=q_to_size(q),
            color="0.35",
            alpha=0.75,
            edgecolors="none",
            transform=ax.transAxes,
            clip_on=False,
            zorder=10,
        )
        ax.text(1.19, y_pos, f"q={q:g}", transform=ax.transAxes, ha="left", va="center", fontsize=8.5, clip_on=False)
    ax.legend(
        handles=[Patch(color=CATEGORY_COLORS[c], label=c) for c in CATEGORY_ORDER],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.28),
        ncol=3,
    )
    save(fig, out_dir / "panel_D_top_gsea_terms")


def category_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sig = df["FDR"] <= 0.05
    for category in CATEGORY_ORDER:
        in_cat = df["Category"] == category
        a = int((sig & in_cat).sum())
        b = int((sig & ~in_cat).sum())
        c = int((~sig & in_cat).sum())
        d = int((~sig & ~in_cat).sum())
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds_ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append((category, a, int(in_cat.sum()), math.log2(odds_ha), p))
    out = pd.DataFrame(rows, columns=["Category", "n_sig", "n_total", "log2_odds_ratio", "p"])
    out["q"] = bh_fdr(out["p"].to_numpy(dtype=float))
    return out


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def plot_panels_e_f(enrich_dir: Path, out_dir: Path) -> None:
    df = load_gsea(enrich_dir / "gsea_results_all.tsv")
    summary = category_summary(df)
    summary.to_csv(enrich_dir / "category_enrichment_fisher.csv", index=False)
    fig, ax = plt.subplots(figsize=(4.8, 3.35))
    plot_df = summary.iloc[::-1]
    ax.barh(plot_df["Category"], plot_df["log2_odds_ratio"], color=[CATEGORY_COLORS[c] for c in plot_df["Category"]])
    ax.axvline(0, color="0.25", lw=0.9)
    ax.set_xlabel(r"$\log_2$(odds ratio)")
    xmin = min(-2.0, float(plot_df["log2_odds_ratio"].min()) - 0.35)
    xmax = max(2.0, float(plot_df["log2_odds_ratio"].max()) + 1.0)
    ax.set_xlim(xmin, xmax)
    for i, row in enumerate(plot_df.itertuples(index=False)):
        xpos = row.log2_odds_ratio + 0.08 if row.log2_odds_ratio >= 0 else 0.08
        ax.text(xpos, i, f"BH q={row.q:.3g}", va="center", ha="left", fontsize=8.5)
    save(fig, out_dir / "panel_E_category_log2_odds")

    cats = [c for c in CATEGORY_ORDER if c in set(df["Category"])]
    data = [df.loc[df["Category"] == c, "NES"].to_numpy(dtype=float) for c in cats]
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    vp = ax.violinplot(data, showmedians=True, showextrema=False)
    for body, cat in zip(vp["bodies"], cats):
        body.set_facecolor(CATEGORY_COLORS[cat])
        body.set_edgecolor("none")
        body.set_alpha(0.88)
    if "cmedians" in vp:
        vp["cmedians"].set_color("black")
    ax.set_xticks(np.arange(1, len(cats) + 1))
    ax.set_xticklabels(cats, rotation=24, ha="right", fontsize=8.5)
    ax.set_ylabel("NES")
    save(fig, out_dir / "panel_F_category_nes_violin")


def plot_panel_g(cell_dir: Path, out_dir: Path) -> None:
    df = pd.read_csv(cell_dir / "eta_celltype_correlations.csv").set_index("cell_type")
    df = df.reindex([c for c in CELL_ORDER if c in df.index]).reset_index()
    labels = df["cell_type"].str.replace("frac_", "", regex=False)
    y = np.arange(len(df))[::-1]
    fig, ax = plt.subplots(figsize=(4.7, 3.9))
    ax.barh(y, df["spearman_rho"], color="#8c61aa")
    ax.axvline(0, color="0.25", lw=0.9)
    xmin = min(-0.24, float(df["spearman_rho"].min()) - 0.06)
    xmax = max(0.42, float(df["spearman_rho"].max()) + 0.12)
    ax.set_xlim(xmin, xmax)
    for yi, row in zip(y, df.itertuples(index=False)):
        ax.text(xmax - 0.01 * (xmax - xmin), yi, f"q={row.p_fdr:.3g}", va="center", ha="right", fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"Spearman $\rho(\eta,\mathrm{CLR(cell\ type)})$")
    save(fig, out_dir / "panel_G_cell_type_correlations")


def plot_panel_h(cell_dir: Path, out_dir: Path) -> None:
    joint = pd.read_csv(cell_dir / "eta_celltype_joint_table.csv")
    stat = pd.read_csv(cell_dir / "eta_monoaminergic_stats.csv").iloc[0]
    x = joint["monoaminergic_score"].to_numpy(dtype=float)
    y = joint["eta"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(4.4, 3.9))
    ax.scatter(x, y, s=58, color="#8c61aa", alpha=0.78, edgecolor="white", linewidth=0.4)
    if len(joint) >= 3:
        slope, intercept, *_ = stats.linregress(x, y)
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 120)
        ax.plot(xx, intercept + slope * xx, color="black", lw=2.0)
    ax.set_xlabel("Monoaminergic score")
    ax.set_ylabel(r"vulnerability axis $\eta$")
    ax.set_title(rf"Spearman $\rho$={stat['spearman_rho']:.3f}, {p_text(float(stat['p_perm']))}")
    save(fig, out_dir / "panel_H_monoaminergic_eta")


def scatter_compare(
    x: pd.Series,
    y: pd.Series,
    xlabel: str,
    ylabel: str,
    out: Path,
    color: str,
    ax: plt.Axes | None = None,
    title_prefix: str | None = None,
) -> None:
    mask = np.isfinite(x.to_numpy(dtype=float)) & np.isfinite(y.to_numpy(dtype=float))
    xv = x.to_numpy(dtype=float)[mask]
    yv = y.to_numpy(dtype=float)[mask]
    r, p = stats.spearmanr(xv, yv)
    lim = float(np.nanquantile(np.abs(np.concatenate([xv, yv])), 0.995))
    lim = max(lim, 0.05)
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(4.2, 4.0))
    else:
        fig = ax.figure
    ax.scatter(xv, yv, s=18, color=color, alpha=0.22, edgecolor="#202020", linewidth=0.18, rasterized=True)
    ax.plot([-lim, lim], [-lim, lim], color="0.45", lw=1.4, ls=(0, (6, 5)))
    slope, intercept, *_ = stats.linregress(xv, yv)
    xx = np.linspace(-lim, lim, 120)
    ax.plot(xx, intercept + slope * xx, color="black", lw=2.2)
    ax.axhline(0, color="0.55", lw=1.1, ls=(0, (6, 5)))
    ax.axvline(0, color="0.55", lw=1.1, ls=(0, (6, 5)))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    title = rf"Spearman $\rho$ = {r:.3f}, {p_text(p)}"
    if title_prefix:
        title = f"{title_prefix}\n{title}"
    ax.set_title(title)
    if own_fig:
        save(fig, out)


def plot_comparisons(results_root: Path, figures_root: Path, filter_key: str) -> None:
    sdir = results_root / filter_key / "striatum" / "transcriptomics"
    hdir = results_root / filter_key / "hippocampus_C1_C4" / "transcriptomics"
    out = figures_root / filter_key / "comparisons"
    s_coef = pd.read_csv(sdir / "gene_parameter_coefficients.csv")
    h_coef = pd.read_csv(hdir / "gene_parameter_coefficients.csv")
    coef = s_coef.merge(h_coef, on="gene", suffixes=("_striatum", "_hippocampus"))
    scatter_compare(
        coef["coef_gamma_striatum"],
        coef["coef_gamma_hippocampus"],
        r"Striatum gene coefficient for $z(\gamma)$",
        r"Hippocampus gene coefficient for $z(\gamma)$",
        out / "panel_A_gene_gamma_coefficient_comparison",
        "#b32f32",
    )
    scatter_compare(
        coef["coef_beta_striatum"],
        coef["coef_beta_hippocampus"],
        r"Striatum gene coefficient for $z(\beta)$",
        r"Hippocampus gene coefficient for $z(\beta)$",
        out / "panel_B_gene_beta_coefficient_comparison",
        "#225ea8",
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0), constrained_layout=True)
    scatter_compare(
        coef["coef_gamma_striatum"],
        coef["coef_gamma_hippocampus"],
        r"Striatum gene coefficient for $z(\gamma)$",
        r"Hippocampus gene coefficient for $z(\gamma)$",
        out / "_unused_gamma",
        "#b32f32",
        ax=axes[0],
        title_prefix="A",
    )
    scatter_compare(
        coef["coef_beta_striatum"],
        coef["coef_beta_hippocampus"],
        r"Striatum gene coefficient for $z(\beta)$",
        r"Hippocampus gene coefficient for $z(\beta)$",
        out / "_unused_beta",
        "#225ea8",
        ax=axes[1],
        title_prefix="B",
    )
    save(fig, out / "panel_AB_gene_beta_gamma_coefficient_comparison")
    s_eta = pd.read_csv(sdir / "gene_eta_correlations.csv")[["gene", "r"]]
    h_eta = pd.read_csv(hdir / "gene_eta_correlations.csv")[["gene", "r"]]
    eta = s_eta.merge(h_eta, on="gene", suffixes=("_striatum", "_hippocampus"))
    scatter_compare(
        eta["r_striatum"],
        eta["r_hippocampus"],
        r"Striatum gene correlation with $\eta$",
        r"Hippocampus gene correlation with $\eta$",
        out / "panel_D_gene_eta_correlation_comparison",
        "#4a3f8f",
    )


def plot_dataset(results_root: Path, figures_root: Path, filter_key: str, dataset: str, omit_regional_panels: bool) -> None:
    base = results_root / filter_key / dataset
    out = figures_root / filter_key / dataset
    plot_panel_a(base / "transcriptomics", out)
    if not omit_regional_panels:
        plot_panel_b(base / "transcriptomics", out)
    plot_panel_d(base / "enrichment", out)
    plot_panels_e_f(base / "enrichment", out)
    plot_panel_g(base / "cell_types", out)
    plot_panel_h(base / "cell_types", out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", default="paper-rf/results/figure_6_7")
    parser.add_argument("--figures-root", default="paper-rf/figures/figure_6_7")
    parser.add_argument("--omit-regional-panels", action="store_true")
    args = parser.parse_args()
    setup_style()
    results_root = Path(args.results_root)
    figures_root = Path(args.figures_root)
    for filter_key in FILTERS:
        for dataset, _, _ in DATASETS:
            plot_dataset(results_root, figures_root, filter_key, dataset, args.omit_regional_panels)
        plot_comparisons(results_root, figures_root, filter_key)
    print(figures_root)


if __name__ == "__main__":
    main()
