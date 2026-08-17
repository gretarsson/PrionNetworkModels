#!/usr/bin/env python3
"""Tighten linked Figure 6/7 biological panels for manual Illustrator layout."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import stats

from plot_manuscript_figure_6_7_panels import (
    CATEGORY_COLORS,
    CATEGORY_ORDER,
    CELL_ORDER,
    bh_fdr,
    load_gmt,
    load_gsea,
    p_text,
    pca_segments,
    segment,
    setup_style,
)


ROOT = Path(__file__).resolve().parents[3]
ENRICH_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "striatum" / "enrichment"
CELL_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "striatum" / "cell_types"
OUT_DIR = ROOT / "paper-rf" / "figures" / "Figure6"
HIPPO_ENRICH_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "hippocampus_C1_C4" / "enrichment"
HIPPO_CELL_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "hippocampus_C1_C4" / "cell_types"
STRIATUM_TRANS_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "striatum" / "transcriptomics"
HIPPO_TRANS_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "hippocampus_C1_C4" / "transcriptomics"
PC1_DIR_CSV = ROOT / "paper-rf" / "figures" / "pooled_z" / "all" / "pc1_direction" / "pc1_direction_comparison.csv"
FIG7_OUT_DIR = ROOT / "paper-rf" / "figures" / "Figure7"
FIG7_EXTRA_OUT_DIR = ROOT / "paper-rf" / "figures" / "figure7"

TERM_LABELS = {
    "Parkinson disease": "Parkinson\ndisease",
    "Oxidative phosphorylation": "Oxidative\nphosphorylation",
    "Alzheimer disease": "Alzheimer\ndisease",
    "Citrate cycle (TCA cycle)": "TCA cycle",
    "Huntington disease": "Huntington\ndisease",
    "Synaptic vesicle cycle": "Synaptic\nvesicle cycle",
    "Terpenoid backbone biosynthesis": "Terpenoid\nbiosynthesis",
    "Pyruvate metabolism": "Pyruvate\nmetabolism",
    "Steroid biosynthesis": "Steroid\nbiosynthesis",
    "Starch and sucrose metabolism": "Starch/sucrose\nmetabolism",
    "Protein processing in endoplasmic reticulum": "Protein processing\nin ER",
    "Non-alcoholic fatty liver disease (NAFLD)": "NAFLD",
    "Valine, leucine and isoleucine degradation": "Valine, leucine and\nisoleucine degradation",
    "Biosynthesis of unsaturated fatty acids": "Biosynthesis of\nunsaturated fatty acids",
    "Primary bile acid biosynthesis": "Primary bile acid\nbiosynthesis",
}

CATEGORY_LABELS = {
    "Metabolism": "Metab.",
    "Protein homeostasis": "Protein\nhomeost.",
    "Synapse": "Synapse",
    "Neurodeg. disease": "Neurodeg.",
    "Other": "Other",
}


def q_to_size_panel_b(q: float) -> float:
    if not np.isfinite(q):
        return 18.0
    score = -math.log(max(float(q), 1e-8))
    return float(np.clip(8 + 8.5 * score, 20, 74))


def save_pdf(fig: plt.Figure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def q_stars(q: float) -> str:
    if not np.isfinite(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def q_text(q: float) -> str:
    if not np.isfinite(q):
        return "q=n/a"
    if q < 0.001:
        return "q<0.001"
    return f"q={q:.3g}"


def p_stars(p: float) -> str:
    return q_stars(p)


def add_bold_title_stars(ax: plt.Axes, stars: str) -> None:
    if not stars:
        return
    fig = ax.figure
    fig.canvas.draw()
    title = ax.title
    renderer = fig.canvas.get_renderer()
    title_box = title.get_window_extent(renderer=renderer)
    ax_box = ax.get_window_extent(renderer=renderer)
    x = (title_box.x1 - ax_box.x0) / ax_box.width + 0.004
    y = (title_box.y0 + 0.5 * title_box.height - ax_box.y0) / ax_box.height
    ax.text(
        x,
        y,
        f"({stars})",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        clip_on=False,
    )


def category_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sig = df["FDR"] <= 0.05
    for category in CATEGORY_ORDER:
        in_cat = df["Category"] == category
        a = int((sig & in_cat).sum())
        b = int((sig & ~in_cat).sum())
        c = int((~sig & in_cat).sum())
        d = int((~sig & ~in_cat).sum())
        odds_ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        _, p = __import__("scipy.stats", fromlist=["fisher_exact"]).fisher_exact(
            [[a, b], [c, d]], alternative="greater"
        )
        rows.append((category, a, int(in_cat.sum()), math.log2(odds_ha), p))
    out = pd.DataFrame(rows, columns=["Category", "n_sig", "n_total", "log2_odds_ratio", "p"])
    out["q"] = bh_fdr(out["p"].to_numpy(dtype=float))
    return out


def top_gsea_terms(enrich_dir: Path, n: int = 10) -> list[str]:
    df = load_gsea(enrich_dir / "gsea_results_all.tsv")
    sig = df[df["FDR"] <= 0.05].copy()
    sig["absNES"] = sig["NES"].abs()
    return sig.sort_values("absNES", ascending=False).head(n)["Term"].tolist()


def replicated_gsea_terms(striatum_dir: Path, hippocampus_dir: Path, n: int = 10) -> list[str]:
    striatum = load_gsea(striatum_dir / "gsea_results_all.tsv").set_index("Term")
    hippocampus = load_gsea(hippocampus_dir / "gsea_results_all.tsv").set_index("Term")
    common = striatum.index.intersection(hippocampus.index)
    replicated = pd.DataFrame(
        {
            "Term": common,
            "Category": hippocampus.loc[common, "Category"].to_numpy(),
            "striatum_NES": striatum.loc[common, "NES"].to_numpy(dtype=float),
            "striatum_q": striatum.loc[common, "FDR"].to_numpy(dtype=float),
            "hippocampus_NES": hippocampus.loc[common, "NES"].to_numpy(dtype=float),
            "hippocampus_q": hippocampus.loc[common, "FDR"].to_numpy(dtype=float),
        }
    )
    replicated = replicated[
        (replicated["striatum_q"] <= 0.05)
        & (replicated["hippocampus_q"] <= 0.05)
        & ((replicated["striatum_NES"] * replicated["hippocampus_NES"]) > 0)
    ].copy()
    replicated = replicated.sort_values("hippocampus_NES", ascending=False)
    replicated.to_csv(hippocampus_dir / "replicated_gsea_terms.csv", index=False)
    return replicated.head(n)["Term"].tolist()


def plot_gsea_top(enrich_dir: Path, out_file: Path, terms: list[str] | None = None) -> None:
    df = load_gsea(enrich_dir / "gsea_results_all.tsv")
    ranked = pd.read_csv(enrich_dir / "ranked_genes.rnk", sep="\t", header=None, names=["gene", "score"])
    ranked["score"] = pd.to_numeric(ranked["score"], errors="coerce")
    ranked = ranked.dropna(subset=["gene", "score"])
    ranked["gene_key"] = ranked["gene"].astype(str).str.upper()
    score_by_gene = ranked.drop_duplicates("gene_key").set_index("gene_key")["score"]
    gene_sets = {
        term: {gene.upper() for gene in genes}
        for term, genes in load_gmt(enrich_dir / "gseapy" / "gene_sets.gmt").items()
    }

    if terms is None:
        terms = top_gsea_terms(enrich_dir, n=10)
    top = df.set_index("Term").reindex(terms).dropna(subset=["NES"]).reset_index().iloc[::-1].copy()

    fig, ax = plt.subplots(figsize=(6.35, 4.35))
    scores = ranked["score"].to_numpy(dtype=float)
    violins = ax.violinplot(
        [scores for _ in range(len(top))],
        positions=np.arange(len(top)),
        vert=False,
        widths=0.7,
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
            ax.scatter(
                member_scores,
                i + rng.uniform(-0.12, 0.12, size=len(member_scores)),
                s=q_to_size_panel_b(float(row.FDR)),
                color=color,
                alpha=0.58,
                edgecolors="none",
                zorder=3,
                rasterized=True,
            )
        ax.text(
            1.01,
            i,
            f"{row.NES:+.2f}",
            va="center",
            ha="left",
            fontsize=9.4,
            transform=ax.get_yaxis_transform(),
            clip_on=False,
        )
        ax.text(
            1.13,
            i,
            q_text(float(row.FDR)),
            va="center",
            ha="left",
            fontsize=8.2,
            transform=ax.get_yaxis_transform(),
            clip_on=False,
        )

    ax.axvline(0, color="0.25", lw=0.9)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([TERM_LABELS.get(t, t) for t in top["TermName"]], fontsize=8.4)
    ax.set_xlabel("r")
    lim = float(np.nanquantile(np.abs(scores), 0.995))
    ax.set_xlim(-lim * 1.08, lim * 1.08)
    ax.text(1.01, len(top) - 0.2, "NES", ha="left", va="bottom", fontsize=10.2, transform=ax.get_yaxis_transform())
    ax.text(1.13, len(top) - 0.2, "FDR", ha="left", va="bottom", fontsize=10.2, transform=ax.get_yaxis_transform())

    fdr_y = -0.155
    ax.text(0.0, fdr_y, "FDR", transform=ax.transAxes, ha="left", va="center", fontsize=9.5, clip_on=False)
    for x_pos, q in zip([0.14, 0.34, 0.58], [0.05, 0.01, 0.001]):
        ax.scatter(
            [x_pos],
            [fdr_y],
            s=q_to_size_panel_b(q),
            color="0.35",
            alpha=0.75,
            edgecolors="none",
            transform=ax.transAxes,
            clip_on=False,
            zorder=10,
        )
        ax.text(x_pos + 0.035, fdr_y, f"q={q:g}", transform=ax.transAxes, ha="left", va="center", fontsize=7.8, clip_on=False)

    ax.legend(
        handles=[Patch(color=CATEGORY_COLORS[c], label=c) for c in CATEGORY_ORDER],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.275),
        ncol=5,
        handlelength=1.05,
        columnspacing=0.85,
        handletextpad=0.35,
        fontsize=7.8,
    )
    save_pdf(fig, out_file)


def plot_gene_coefficient_cloud(trans_dir: Path, out_file: Path) -> None:
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
    sc = ax.scatter(
        x,
        y,
        c=df["r"],
        cmap="coolwarm",
        vmin=-max_abs,
        vmax=max_abs,
        s=25,
        alpha=0.58,
        edgecolors="#454545",
        linewidths=0.22,
        rasterized=True,
    )
    p2a, p2b = segment(center, pc2, x, y, 0.23)
    p1a, p1b = segment(center, pc1, x, y, 0.36)
    ax.plot([p2a[0], p2b[0]], [p2a[1], p2b[1]], color="#8c8c8c", lw=5, solid_capstyle="round")
    ax.plot([p1a[0], p1b[0]], [p1a[1], p1b[1]], color="black", lw=5, solid_capstyle="round")
    ax.set_xlabel(r"gene coefficient for $z(\beta)$", fontsize=12.8)
    ax.set_ylabel(r"gene coefficient for $z(\gamma)$", fontsize=12.8)
    ax.set_title(rf"PC1 defines $\eta = {pc1[0]:.2f}z(\beta) + {pc1[1]:.2f}z(\gamma)$", fontsize=12.5)
    ax.legend(
        handles=[
            Line2D([0], [0], color="#8c8c8c", lw=4.8, solid_capstyle="round", label="PC2"),
            Line2D([0], [0], color="black", lw=4.8, solid_capstyle="round", label="PC1"),
        ],
        frameon=False,
        loc="lower left",
        handlelength=1.9,
        handletextpad=0.32,
        fontsize=11.2,
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.025)
    cbar.set_label(r"corr(gene, $\eta$)", fontsize=12.2)
    cbar.ax.tick_params(labelsize=10.5)
    ax.text(
        0.98,
        0.97,
        f"PC1 {100 * pc1_row['explained_variance_ratio']:.1f}%",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11.4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 3.5},
    )
    save_pdf(fig, out_file)


def plot_category_odds(enrich_dir: Path, out_file: Path) -> None:
    df = load_gsea(enrich_dir / "gsea_results_all.tsv")
    summary = category_summary(df)
    summary.to_csv(enrich_dir / "category_enrichment_fisher.csv", index=False)
    plot_df = summary.iloc[::-1].copy()

    fig, ax = plt.subplots(figsize=(3.55, 2.85))
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["log2_odds_ratio"], color=[CATEGORY_COLORS[c] for c in plot_df["Category"]])
    ax.axvline(0, color="0.25", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([CATEGORY_LABELS[c] for c in plot_df["Category"]], fontsize=8.5)
    ax.set_xlabel(r"$\log_2$(odds ratio)", fontsize=10.5)
    xmin = min(-2.0, float(plot_df["log2_odds_ratio"].min()) - 0.3)
    xmax = max(2.0, float(plot_df["log2_odds_ratio"].max()) + 0.55)
    ax.set_xlim(xmin, xmax)
    span = xmax - xmin
    for yi, row in zip(y, plot_df.itertuples(index=False)):
        label = f"q={row.q:.3g}"
        stars = q_stars(float(row.q))
        if row.log2_odds_ratio >= 0:
            xpos = -0.045 * span
            ha = "right"
            star_x = row.log2_odds_ratio + 0.035 * span
            star_ha = "left"
        else:
            xpos = 0.05 * span
            ha = "left"
            star_x = row.log2_odds_ratio - 0.035 * span
            star_ha = "right"
        ax.text(xpos, yi, label, va="center", ha=ha, fontsize=7.7)
        if stars:
            ax.text(star_x, yi, stars, va="center", ha=star_ha, fontsize=14.5, fontweight="bold")
    save_pdf(fig, out_file)


def plot_category_nes_violin(enrich_dir: Path, out_file: Path) -> None:
    df = load_gsea(enrich_dir / "gsea_results_all.tsv")
    cats = [c for c in CATEGORY_ORDER if c in set(df["Category"])]
    data = [df.loc[df["Category"] == c, "NES"].dropna().to_numpy(dtype=float) for c in cats]

    fig, ax = plt.subplots(figsize=(4.25, 2.9))
    vp = ax.violinplot(data, showmeans=False, showmedians=True, showextrema=False)
    for body, cat in zip(vp["bodies"], cats):
        body.set_facecolor(CATEGORY_COLORS[cat])
        body.set_edgecolor("none")
        body.set_alpha(0.82)
    if "cmedians" in vp:
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.0)
    ax.axhline(0, color="0.75", lw=0.9)
    ax.set_ylabel("NES", fontsize=10.5)
    ax.set_xticks(np.arange(1, len(cats) + 1))
    ax.set_xticklabels(cats, rotation=22, ha="right", fontsize=8.4)
    save_pdf(fig, out_file)


def plot_cell_type_correlations(cell_dir: Path, out_file: Path) -> None:
    df = pd.read_csv(cell_dir / "eta_celltype_correlations.csv").set_index("cell_type")
    df = df.reindex([c for c in CELL_ORDER if c in df.index]).reset_index()
    labels = df["cell_type"].str.replace("frac_", "", regex=False)
    y = np.arange(len(df))[::-1]

    fig, ax = plt.subplots(figsize=(4.35, 3.65))
    ax.barh(y, df["spearman_rho"], color="#8c61aa")
    ax.axvline(0, color="0.25", lw=0.9)
    xmin = min(-0.24, float(df["spearman_rho"].min()) - 0.06)
    xmax = max(0.42, float(df["spearman_rho"].max()) + 0.08)
    ax.set_xlim(xmin, xmax)
    span = xmax - xmin
    for yi, row in zip(y, df.itertuples(index=False)):
        label = f"q={row.p_fdr:.3g}"
        stars = q_stars(float(row.p_fdr))
        if row.spearman_rho >= 0:
            xpos = -0.04 * span
            ha = "right"
            star_x = row.spearman_rho + 0.028 * span
            star_ha = "left"
        else:
            xpos = 0.045 * span
            ha = "left"
            star_x = row.spearman_rho - 0.028 * span
            star_ha = "right"
        ax.text(xpos, yi, label, va="center", ha=ha, fontsize=7.7)
        if stars:
            ax.text(star_x, yi, stars, va="center", ha=star_ha, fontsize=14.5, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r"Spearman $\rho(\eta,\mathrm{CLR(cell\ type)})$", fontsize=10.5)
    save_pdf(fig, out_file)


def plot_monoaminergic_eta(cell_dir: Path, out_file: Path) -> None:
    joint = pd.read_csv(cell_dir / "eta_celltype_joint_table.csv")
    stat = pd.read_csv(cell_dir / "eta_monoaminergic_stats.csv").iloc[0]
    x = joint["monoaminergic_score"].to_numpy(dtype=float)
    y = joint["eta"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(4.15, 3.6))
    ax.scatter(x, y, s=48, color="#8c61aa", alpha=0.78, edgecolor="white", linewidth=0.4)
    if len(joint) >= 3:
        slope, intercept, *_ = stats.linregress(x, y)
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 120)
        ax.plot(xx, intercept + slope * xx, color="black", lw=2.0)
    ax.set_xlabel("Monoaminergic score", fontsize=10.5)
    ax.set_ylabel(r"vulnerability axis $\eta$", fontsize=10.5)
    p_perm = float(stat["p_perm"])
    ax.set_title(
        rf"Spearman $\rho$={stat['spearman_rho']:.3f}, {p_text(p_perm)}",
        fontsize=10.5,
        pad=7,
    )
    add_bold_title_stars(ax, p_stars(p_perm))
    save_pdf(fig, out_file)


def plot_gene_eta_comparison(out_file: Path) -> None:
    s_eta = pd.read_csv(STRIATUM_TRANS_DIR / "gene_eta_correlations.csv")[["gene", "r"]]
    h_eta = pd.read_csv(HIPPO_TRANS_DIR / "gene_eta_correlations.csv")[["gene", "r"]]
    df = s_eta.merge(h_eta, on="gene", suffixes=("_striatum", "_hippocampus"))
    x = df["r_striatum"].to_numpy(dtype=float)
    y = df["r_hippocampus"].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    rho, p = stats.spearmanr(x, y)
    lim = float(np.nanquantile(np.abs(np.concatenate([x, y])), 0.995))
    lim = max(lim, 0.05)

    fig, ax = plt.subplots(figsize=(4.15, 3.9))
    ax.scatter(x, y, s=18, color="#1f8a8a", alpha=0.24, edgecolor="#123f3f", linewidth=0.18, rasterized=True)
    slope, intercept, *_ = stats.linregress(x, y)
    xx = np.linspace(-lim, lim, 120)
    ax.plot(xx, intercept + slope * xx, color="black", lw=2.2)
    ax.axhline(0, color="0.55", lw=1.1, ls=(0, (6, 5)))
    ax.axvline(0, color="0.55", lw=1.1, ls=(0, (6, 5)))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"Striatum gene correlation with $\eta$", fontsize=10.5)
    ax.set_ylabel(r"Hippocampus gene correlation with $\eta$", fontsize=10.5)
    ax.set_title(rf"Spearman $\rho$={rho:.3f}, {p_text(float(p))}", fontsize=10.5, pad=7)
    add_bold_title_stars(ax, p_stars(float(p)))
    save_pdf(fig, out_file)


def plot_pc1_direction(out_file: Path) -> None:
    df = pd.read_csv(PC1_DIR_CSV)
    fig, ax = plt.subplots(figsize=(3.35, 3.25))
    theta = np.linspace(0, np.pi, 240)
    ax.plot(np.cos(theta), np.sin(theta), color="0.88", lw=1)
    ax.axhline(0, color="0.86", lw=0.9)
    ax.axvline(0, color="0.86", lw=0.9)
    for row in df.itertuples(index=False):
        ax.annotate(
            "",
            xy=(row.loading_beta, row.loading_gamma),
            xytext=(0, 0),
            arrowprops={"arrowstyle": "-|>", "lw": 2.2, "color": row.color, "mutation_scale": 12},
        )
        ax.plot([], [], color=row.color, lw=2.2, label=f"{row.dataset} PC1")
    ax.text(0.06, 0.90, f"angle = {df.angle_degrees.iloc[0]:.1f} deg", transform=ax.transAxes, ha="left", va="top", fontsize=10.2)
    ax.set_xlim(-1.04, 0.10)
    ax.set_ylim(-0.03, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"PC1 loading on $z(\beta)$", fontsize=10.5)
    ax.set_ylabel(r"PC1 loading on $z(\gamma)$", fontsize=10.5)
    ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.28, 0.84),
        fontsize=10.2,
        handlelength=1.55,
        handletextpad=0.5,
        labelspacing=0.3,
        borderaxespad=0.0,
    )
    save_pdf(fig, out_file)


def write_figure7_panels(out_dir: Path, replicated_terms: list[str]) -> None:
    plot_gene_coefficient_cloud(HIPPO_TRANS_DIR, out_dir / "A_gene_coefficient_cloud.pdf")
    plot_gsea_top(HIPPO_ENRICH_DIR, out_dir / "D_gsea_top_pathways.pdf", terms=replicated_terms)
    plot_gsea_top(HIPPO_ENRICH_DIR, out_dir / "G_gsea_hippocampus_top_pathways.pdf")
    plot_pc1_direction(out_dir / "B_pc1_direction_comparison.pdf")
    plot_gene_eta_comparison(out_dir / "C_gene_eta_correlation_comparison.pdf")
    plot_cell_type_correlations(HIPPO_CELL_DIR, out_dir / "E_cell_type_correlations.pdf")
    plot_monoaminergic_eta(HIPPO_CELL_DIR, out_dir / "F_monoaminergic_eta.pdf")


def main() -> None:
    setup_style()
    replicated_terms = replicated_gsea_terms(ENRICH_DIR, HIPPO_ENRICH_DIR, n=10)

    plot_gene_coefficient_cloud(STRIATUM_TRANS_DIR, OUT_DIR / "A_gene_coefficient_cloud.pdf")
    plot_gsea_top(ENRICH_DIR, OUT_DIR / "B_gsea_top_pathways.pdf")
    plot_category_odds(ENRICH_DIR, OUT_DIR / "C_gsea_category_odds.pdf")
    plot_category_nes_violin(ENRICH_DIR, OUT_DIR / "D_gsea_category_nes_violin.pdf")
    plot_cell_type_correlations(CELL_DIR, OUT_DIR / "E_cell_type_correlations.pdf")
    plot_monoaminergic_eta(CELL_DIR, OUT_DIR / "F_monoaminergic_eta.pdf")

    for fig7_dir in [FIG7_OUT_DIR, FIG7_EXTRA_OUT_DIR]:
        write_figure7_panels(fig7_dir, replicated_terms)
    print(OUT_DIR)
    print(FIG7_OUT_DIR)
    print(FIG7_EXTRA_OUT_DIR)


if __name__ == "__main__":
    main()
