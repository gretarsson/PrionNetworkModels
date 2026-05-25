#!/usr/bin/env python3
"""Generate manuscript-style biological figure panels from paper/results tables."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import fisher_exact


CATEGORY_ORDER = [
    "Metabolism",
    "Protein homeostasis",
    "Synaptic function",
    "Neurodegenerative disease",
    "Other",
]

CATEGORY_COLORS = {
    "Metabolism": "#2CA02C",
    "Protein homeostasis": "#FF7F0E",
    "Synaptic function": "#9467BD",
    "Neurodegenerative disease": "#C44E52",
    "Other": "#7F7F7F",
}


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 300,
        "savefig.dpi": 300,
    })


def save(fig, out_base: Path) -> None:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def categorize_term(term_name: str) -> str:
    s = str(term_name).lower()
    if any(k in s for k in ["parkinson", "alzheimer", "huntington", "prion", "amyotrophic", "neurodegenerative", "als"]):
        return "Neurodegenerative disease"
    if any(k in s for k in ["synapse", "synaptic", "neurotransmitter", "vesicle", "long-term potentiation", "long-term depression"]):
        return "Synaptic function"
    if any(k in s for k in ["proteasome", "ubiquitin", "autophagy", "lysosome", "endoplasmic reticulum", "chaperone", "protein folding", "protein processing", "ribosome", "translation", "mrna", "degradation", "quality control", "mitophagy"]):
        return "Protein homeostasis"
    if any(k in s for k in ["metabolism", "metabolic", "oxidative phosphorylation", "mitochond", "glycolysis", "tca", "citrate cycle", "respiratory chain", "fatty acid", "lipid", "cholesterol", "biosynthesis", "amino acid", "nucleotide", "energy"]):
        return "Metabolism"
    return "Other"


def parse_term(term: str) -> tuple[str, str]:
    if "__" in str(term):
        library, name = str(term).split("__", 1)
    else:
        library, name = "UNKNOWN", str(term)
    return library, name.replace("_", " ")


def find_fdr_col(df: pd.DataFrame) -> str:
    for col in ["FDR q-val", "FDR", "fdr", "qval", "Adjusted P-value"]:
        if col in df.columns:
            return col
    raise ValueError(f"Could not find FDR column in {list(df.columns)}")


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * len(pvalues) / np.arange(1, len(pvalues) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def load_gsea(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    fdr_col = find_fdr_col(df)
    libraries, names = zip(*[parse_term(t) for t in df["Term"].astype(str)])
    df = df.copy()
    df["Library"] = libraries
    df["TermName"] = names
    df["FDRcol"] = pd.to_numeric(df[fdr_col], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df["Category"] = df["TermName"].map(categorize_term)
    return df.dropna(subset=["NES", "FDRcol"])


def plot_pca_coefficients(result_dir: Path, out_dir: Path, label: str) -> None:
    coefs = pd.read_csv(result_dir / "gene_parameter_coefficients.csv")
    corr = pd.read_csv(result_dir / "gene_eta_correlations.csv")[["gene", "r"]]
    pca = pd.read_csv(result_dir / "pca_summary.csv").iloc[0]
    df = coefs.merge(corr, on="gene", how="left")
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    sc = ax.scatter(df["coef_beta"], df["coef_gamma"], c=df["r"], cmap="coolwarm", s=5, alpha=0.65, rasterized=True)
    scale = max(df["coef_beta"].abs().quantile(0.995), df["coef_gamma"].abs().quantile(0.995))
    ax.plot([0, pca["loading_beta"] * scale], [0, pca["loading_gamma"] * scale], color="black", lw=2.0)
    ax.axhline(0, color="0.75", lw=0.7)
    ax.axvline(0, color="0.75", lw=0.7)
    ax.set_xlabel(r"coefficient for $z(\beta)$")
    ax.set_ylabel(r"coefficient for $z(\gamma)$")
    ax.set_title(f"{label}: gene coefficient PCA")
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r"corr(gene, $\eta$)")
    style_axis(ax)
    save(fig, out_dir / "pca_gene_coefficients")


def plot_beta_gamma_eta(result_dir: Path, out_dir: Path, label: str) -> None:
    df = pd.read_csv(result_dir / "region_axis.csv")
    pca = pd.read_csv(result_dir / "pca_summary.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    sc = ax.scatter(df["z_beta"], df["z_gamma"], c=df["eta"], cmap="viridis", s=45, edgecolor="white", linewidth=0.4)
    ax.axhline(0, color="0.75", lw=0.7)
    ax.axvline(0, color="0.75", lw=0.7)
    ax.set_xlabel(r"$z(\beta)$")
    ax.set_ylabel(r"$z(\gamma)$")
    ax.set_title(f"{label}: vulnerability axis")
    ax.text(
        0.02,
        0.98,
        rf"$\eta={pca['loading_beta']:.2f}z(\beta)+{pca['loading_gamma']:.2f}z(\gamma)$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r"$\eta$")
    style_axis(ax)
    save(fig, out_dir / "beta_gamma_colored_by_eta")


def plot_celltype_bar(cell_dir: Path, out_dir: Path) -> None:
    df = pd.read_csv(cell_dir / "eta_celltype_correlations.csv")
    df = df.sort_values("spearman_rho")
    labels = [c.replace("frac_", "") for c in df["cell_type"]]
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    colors = ["#C44E52" if x < 0 else "#4C72B0" for x in df["spearman_rho"]]
    ax.barh(labels, df["spearman_rho"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(r"Spearman $\rho$ with $\eta$")
    style_axis(ax)
    save(fig, out_dir / "celltype_eta_correlations")


def plot_monoaminergic(cell_dir: Path, out_dir: Path, label: str) -> None:
    joint = pd.read_csv(cell_dir / "eta_celltype_joint_table.csv")
    stats_df = pd.read_csv(cell_dir / "eta_monoaminergic_stats.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(4.0, 3.6))
    ax.scatter(joint["eta"], joint["monoaminergic_score"], s=50, color="#4C72B0", edgecolor="white", linewidth=0.4)
    slope, intercept, *_ = stats.linregress(joint["eta"], joint["monoaminergic_score"])
    x = np.linspace(joint["eta"].min(), joint["eta"].max(), 100)
    ax.plot(x, intercept + slope * x, color="black", lw=1.5)
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("monoaminergic score")
    ax.set_title(label)
    ax.text(
        0.04,
        0.96,
        rf"$\rho={stats_df['spearman_rho']:.2f}$" + f"\nperm p={stats_df['p_perm']:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
    )
    style_axis(ax)
    save(fig, out_dir / "monoaminergic_score_vs_eta")


def plot_gsea_dot(gsea_path: Path, out_dir: Path) -> None:
    if not gsea_path.exists():
        return
    df = load_gsea(gsea_path)
    df = df[df["FDRcol"] <= 0.05].copy()
    if df.empty:
        return
    df["absNES"] = df["NES"].abs()
    df = df.sort_values("absNES", ascending=False).head(14)
    labels = ["\n".join(textwrap.wrap(x, width=30)) for x in df["TermName"]]
    weights = -np.log10(df["FDRcol"].clip(lower=1e-300))
    sizes = 30 + 120 * weights / max(weights.max(), 1)
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(5.8, max(3.0, 0.32 * len(df) + 1.2)))
    ax.scatter(df["NES"], y, s=sizes, color="black")
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("NES")
    style_axis(ax)
    save(fig, out_dir / "gsea_dotplot_top_absNES")


def category_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sig = df["FDRcol"] <= 0.05
    for cat in CATEGORY_ORDER:
        in_cat = df["Category"] == cat
        a = int((sig & in_cat).sum())
        b = int((sig & ~in_cat).sum())
        c = int((~sig & in_cat).sum())
        d = int((~sig & ~in_cat).sum())
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds_ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append((cat, a, int(in_cat.sum()), odds, np.log2(odds_ha), p))
    out = pd.DataFrame(rows, columns=["Category", "n_sig", "n_total", "odds_ratio", "log2_odds_ratio_ha", "p"])
    out["p_fdr"] = bh_fdr(out["p"].to_numpy())
    return out


def plot_gsea_categories(gsea_path: Path, out_dir: Path) -> None:
    if not gsea_path.exists():
        return
    df = load_gsea(gsea_path)
    summary = category_summary(df)
    summary.to_csv(out_dir / "category_enrichment_fisher.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    plot_df = summary.sort_values("log2_odds_ratio_ha")
    ax.barh(plot_df["Category"], plot_df["log2_odds_ratio_ha"], color=[CATEGORY_COLORS[c] for c in plot_df["Category"]])
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(r"$\log_2$(odds ratio)")
    style_axis(ax)
    save(fig, out_dir / "category_fisher_log2odds")

    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    cats = [c for c in CATEGORY_ORDER if c in set(df["Category"])]
    data = [df.loc[df["Category"] == c, "NES"].to_numpy() for c in cats]
    vp = ax.violinplot(data, showmedians=True, showextrema=False)
    for body, cat in zip(vp["bodies"], cats):
        body.set_facecolor(CATEGORY_COLORS[cat])
        body.set_edgecolor("none")
        body.set_alpha(0.85)
    if "cmedians" in vp:
        vp["cmedians"].set_color("black")
    ax.set_xticks(np.arange(1, len(cats) + 1))
    ax.set_xticklabels(cats, rotation=25, ha="right")
    ax.set_ylabel("NES")
    style_axis(ax)
    save(fig, out_dir / "category_nes_violin")


def plot_axis_comparison(compare_dir: Path, out_dir: Path) -> None:
    comp = compare_dir / "gene_eta_comparison.csv"
    if not comp.exists():
        return
    df = pd.read_csv(comp)
    summary = pd.read_csv(compare_dir / "axis_comparison_summary.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    ax.scatter(df["r_striatum"], df["r_hippocampus"], s=5, alpha=0.4, color="0.25", rasterized=True)
    lim = max(abs(df["r_striatum"]).quantile(0.995), abs(df["r_hippocampus"]).quantile(0.995))
    ax.plot([-lim, lim], [-lim, lim], color="0.6", lw=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"striatal corr(gene, $\eta$)")
    ax.set_ylabel(r"hippocampal corr(gene, $\eta$)")
    ax.text(0.04, 0.96, rf"$r={summary['gene_eta_pearson_r']:.2f}$", transform=ax.transAxes, ha="left", va="top")
    style_axis(ax)
    save(fig, out_dir / "gene_eta_striatum_vs_hippocampus")


def comparison_text(x: np.ndarray, y: np.ndarray) -> str:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return f"n={mask.sum()}"
    r, p = stats.pearsonr(x[mask], y[mask])
    return f"r={r:.2f}\np={p:.2g}\nn={mask.sum()}"


def shared_lims(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    vals = pd.concat([x, y], ignore_index=True)
    lo = vals.quantile(0.01)
    hi = vals.quantile(0.99)
    pad = 0.06 * max(hi - lo, 1e-12)
    return float(lo - pad), float(hi + pad)


def plot_gene_parameter_comparison(root: Path, out_dir: Path) -> None:
    str_path = root / "transcriptomics/striatum/gene_parameter_coefficients.csv"
    hip_path = root / "transcriptomics/hippocampus/gene_parameter_coefficients.csv"
    if not str_path.exists() or not hip_path.exists():
        return

    striatum = pd.read_csv(str_path)[["gene", "coef_beta", "coef_gamma"]]
    hippocampus = pd.read_csv(hip_path)[["gene", "coef_beta", "coef_gamma"]]
    df = striatum.merge(hippocampus, on="gene", suffixes=("_striatum", "_hippocampus"))

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.6))
    for ax, param, title in zip(axes, ["beta", "gamma"], [r"$\beta$ coefficients", r"$\gamma$ coefficients"]):
        x = df[f"coef_{param}_striatum"]
        y = df[f"coef_{param}_hippocampus"]
        ax.scatter(x, y, s=5, alpha=0.35, color="0.25", rasterized=True)
        lo, hi = shared_lims(x, y)
        ax.plot([lo, hi], [lo, hi], color="0.6", lw=1.0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel(f"striatal {title}")
        ax.set_ylabel(f"hippocampal {title}")
        ax.text(
            0.04,
            0.96,
            comparison_text(x.to_numpy(), y.to_numpy()),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
        style_axis(ax)
    fig.tight_layout()
    save(fig, out_dir / "gene_beta_gamma_coefficients_striatum_vs_hippocampus")


def load_parameter_family(root: Path, dataset: str, family: str) -> pd.DataFrame:
    path = root / f"parameters/{dataset}_diff_rf/{family}.csv"
    df = pd.read_csv(path)
    return df.rename(
        columns={
            "mean_post": family,
            "ks_pvalue": f"{family}_ks_pvalue",
            "updated": f"{family}_updated",
        }
    )[["region", family, f"{family}_ks_pvalue", f"{family}_updated"]]


def load_regional_parameters(root: Path, dataset: str) -> pd.DataFrame:
    beta = load_parameter_family(root, dataset, "beta")
    gamma = load_parameter_family(root, dataset, "gamma")
    df = beta.merge(gamma, on="region", how="inner")
    df["region_base"] = df["region"].astype(str).map(lambda x: x[1:] if len(x) > 1 and x[0] in {"i", "c"} else x)
    df["hemi"] = df["region"].astype(str).str[0]
    return df


def plot_regional_parameter_comparison(root: Path, out_dir: Path) -> None:
    try:
        striatum = load_regional_parameters(root, "striatum")
        hippocampus = load_regional_parameters(root, "hippocampus")
    except FileNotFoundError:
        return

    df = striatum.merge(hippocampus, on="region", suffixes=("_striatum", "_hippocampus"))
    df["used_in_both"] = (
        (df["beta_striatum"] > 0)
        & (df["beta_hippocampus"] > 0)
        & (df["beta_ks_pvalue_striatum"] < 0.001)
        & (df["gamma_ks_pvalue_striatum"] < 0.001)
        & (df["beta_ks_pvalue_hippocampus"] < 0.001)
        & (df["gamma_ks_pvalue_hippocampus"] < 0.001)
    )

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.6))
    for ax, param, title in zip(axes, ["beta", "gamma"], [r"$\beta$", r"$\gamma$"]):
        x = df[f"{param}_striatum"]
        y = df[f"{param}_hippocampus"]
        ax.scatter(x, y, s=16, alpha=0.28, color="0.65", linewidth=0)
        highlight = df["used_in_both"]
        if highlight.any():
            ax.scatter(
                x[highlight],
                y[highlight],
                s=34,
                color="#4C72B0",
                edgecolor="white",
                linewidth=0.4,
                label="used in both",
            )
        lo, hi = shared_lims(x, y)
        ax.plot([lo, hi], [lo, hi], color="0.5", lw=1.0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel(f"striatal {title}")
        ax.set_ylabel(f"hippocampal {title}")
        ax.text(
            0.04,
            0.96,
            comparison_text(x.to_numpy(), y.to_numpy()),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
        style_axis(ax)
    if df["used_in_both"].any():
        axes[0].legend(frameon=False, loc="lower right", fontsize=8)
    fig.tight_layout()
    save(fig, out_dir / "regional_beta_gamma_striatum_vs_hippocampus")


def plot_filtering_robustness(root: Path, dataset: str, out_dir: Path, label: str) -> None:
    variants = [
        ("all", "All regions", root / f"transcriptomics/{dataset}_filter_all"),
        ("beta_positive", r"$\beta>0$", root / f"transcriptomics/{dataset}_filter_beta_positive"),
        ("updated", r"$\beta>0$ + updated", root / f"transcriptomics/{dataset}"),
    ]
    variants = [(key, name, path) for key, name, path in variants if (path / "pca_summary.csv").exists()]
    if len(variants) < 2:
        return

    colors = {
        "all": "#7F7F7F",
        "beta_positive": "#4C72B0",
        "updated": "#C44E52",
    }

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    load_rows = []
    for key, name, path in variants:
        pca = pd.read_csv(path / "pca_summary.csv").iloc[0]
        load_rows.append(
            {
                "key": key,
                "name": name,
                "loading_beta": pca["loading_beta"],
                "loading_gamma": pca["loading_gamma"],
                "explained": pca["explained_variance_ratio"],
                "n_regions": int(pca["n_regions"]),
            }
        )
        axes[0].arrow(
            0,
            0,
            pca["loading_beta"],
            pca["loading_gamma"],
            width=0.006,
            head_width=0.05,
            length_includes_head=True,
            color=colors.get(key, "0.3"),
            alpha=0.95,
            label=f"{name} (n={int(pca['n_regions'])})",
        )
    axes[0].axhline(0, color="0.8", lw=0.8)
    axes[0].axvline(0, color="0.8", lw=0.8)
    axes[0].set_xlim(-1.05, 1.05)
    axes[0].set_ylim(-1.05, 1.05)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel(r"PC1 loading on $z(\beta)$")
    axes[0].set_ylabel(r"PC1 loading on $z(\gamma)$")
    axes[0].set_title(f"{label}: PC1 directions")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    style_axis(axes[0])

    for key, name, path in variants:
        region = pd.read_csv(path / "region_axis.csv")
        axes[1].scatter(
            region["z_beta"],
            region["z_gamma"],
            s=24,
            alpha=0.5 if key != "updated" else 0.9,
            color=colors.get(key, "0.3"),
            edgecolor="white" if key == "updated" else "none",
            linewidth=0.3,
            label=name,
        )
    axes[1].axhline(0, color="0.8", lw=0.8)
    axes[1].axvline(0, color="0.8", lw=0.8)
    axes[1].set_xlabel(r"$z(\beta)$")
    axes[1].set_ylabel(r"$z(\gamma)$")
    axes[1].set_title(f"{label}: retained regions")
    style_axis(axes[1])
    fig.tight_layout()
    save(fig, out_dir / f"{dataset}_pca_filtering_robustness")

    pd.DataFrame(load_rows).to_csv(out_dir / f"{dataset}_pca_filtering_summary.csv", index=False)


def plot_dataset(result_dir: Path, cell_dir: Path, enrichment_dir: Path | None, out_dir: Path, label: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_pca_coefficients(result_dir, out_dir, label)
    plot_beta_gamma_eta(result_dir, out_dir, label)
    plot_celltype_bar(cell_dir, out_dir)
    plot_monoaminergic(cell_dir, out_dir, label)
    if enrichment_dir is not None:
        gsea_path = enrichment_dir / "gsea_results_all.tsv"
        plot_gsea_dot(gsea_path, out_dir)
        plot_gsea_categories(gsea_path, out_dir)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", default="paper/results")
    parser.add_argument("--out-dir", default="paper/figures/biological")
    parser.add_argument("--regional-only", action="store_true")
    args = parser.parse_args()

    setup_style()
    root = Path(args.results_root)
    out = Path(args.out_dir)

    if args.regional_only:
        plot_regional_parameter_comparison(root, out)
    else:
        plot_dataset(
            root / "transcriptomics/striatum",
            root / "cell_types/striatum",
            root / "enrichment/striatum",
            out / "striatum",
            "Striatum",
        )
        plot_dataset(
            root / "transcriptomics/hippocampus",
            root / "cell_types/hippocampus",
            root / "enrichment/hippocampus",
            out / "hippocampus",
            "Hippocampus",
        )
        plot_axis_comparison(root / "transcriptomics/striatum_vs_hippocampus", out / "comparison")
        plot_gene_parameter_comparison(root, out / "comparison")
        plot_filtering_robustness(root, "striatum", out / "filtering", "Striatum")
        plot_filtering_robustness(root, "hippocampus", out / "filtering", "Hippocampus")
    print(out)


if __name__ == "__main__":
    main()
