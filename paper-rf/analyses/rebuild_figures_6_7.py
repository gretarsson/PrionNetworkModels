#!/usr/bin/env python3
"""Rebuild clean Figure 6/7 biological outputs.

This script uses the seed-included striatum DIFF-RF run and the converged
retained hippocampus global-prior DIFF-RF mode as the source inference.
It writes a compact, navigable set of results and independent figure panels
for the all-region and beta-positive editions.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import textwrap
from dataclasses import dataclass
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr
from sklearn.decomposition import PCA


PARAM_RE = re.compile(r"^(beta|gamma)\[(\d+)\]$")
MA_CELLTYPES = ["frac_Dopa", "frac_Nora", "frac_Sero", "frac_Hist"]


@dataclass(frozen=True)
class DatasetSpec:
    key: str
    label: str
    run_dir: Path
    observations: Path
    color: str


DATASETS = [
    DatasetSpec(
        "striatum",
        "Striatum",
        Path("runs/striatum_DIFF-RF_RETRO"),
        Path("paper-rf/data/striatum/observations.csv"),
        "#2b6cb0",
    ),
    DatasetSpec(
        "hippocampus_C3_C4",
        "Hippocampus",
        Path("runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C3_C4"),
        Path("paper-rf/data/hippocampus/observations.csv"),
        "#c2410c",
    ),
]


def strip_hemi(region: str) -> str:
    region = str(region).strip()
    return region[1:] if len(region) > 1 and region[0] in {"i", "c"} else region


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def load_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path)
    expr = expr.rename(columns={expr.columns[0]: "region_base"})
    expr["region_base"] = expr["region_base"].astype(str)
    return expr


def export_parameter_tables(spec: DatasetSpec, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "julia",
            "--project=.",
            "paper-rf/analyses/model_parameters/export_parameter_tables.jl",
            "--run",
            str(spec.run_dir),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
    )


def load_parameters(param_dir: Path) -> pd.DataFrame:
    beta = pd.read_csv(param_dir / "beta.csv")
    gamma = pd.read_csv(param_dir / "gamma.csv")
    beta = beta.rename(columns={"mean_post": "beta", "ks_pvalue": "beta_ks_pvalue", "updated": "beta_updated"})
    gamma = gamma.rename(columns={"mean_post": "gamma", "ks_pvalue": "gamma_ks_pvalue", "updated": "gamma_updated"})
    params = beta[["region", "beta", "beta_ks_pvalue", "beta_updated"]].merge(
        gamma[["region", "gamma", "gamma_ks_pvalue", "gamma_updated"]],
        on="region",
        how="inner",
    )
    params["region_base"] = params["region"].map(strip_hemi)
    params["hemi"] = params["region"].astype(str).str[0]
    return params


def regress_gene_coefficients(df: pd.DataFrame, gene_cols: list[str]) -> pd.DataFrame:
    x = np.column_stack([np.ones(len(df)), df["z_beta"].to_numpy(dtype=float), df["z_gamma"].to_numpy(dtype=float)])
    rows = []
    for gene in gene_cols:
        y = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if mask.sum() < 4:
            continue
        coef, *_ = np.linalg.lstsq(x[mask], y[mask], rcond=None)
        r_eta_placeholder = np.nan
        rows.append((gene, coef[1], coef[2], r_eta_placeholder, int(mask.sum())))
    return pd.DataFrame(rows, columns=["gene", "coef_beta", "coef_gamma", "pc1_score", "n_used"])


def gene_correlations(df: pd.DataFrame, gene_cols: list[str], value_col: str, out_name: str) -> pd.DataFrame:
    x = df[value_col].to_numpy(dtype=float)
    rows = []
    for gene in gene_cols:
        y = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 5:
            continue
        r, p = stats.pearsonr(y[mask], x[mask])
        rows.append((gene, r, p, int(mask.sum())))
    out = pd.DataFrame(rows, columns=["gene", "r", "p_un", "n_used"])
    out["p_fdr"] = bh_fdr(out["p_un"].to_numpy())
    out = out.rename(columns={"r": out_name})
    return out.sort_values(out_name, ascending=False)


def run_transcriptomics(
    spec: DatasetSpec,
    filter_key: str,
    params: pd.DataFrame,
    expression: pd.DataFrame,
    results_dir: Path,
    figures_dir: Path,
) -> dict[str, float | str | int]:
    keep = np.isfinite(params["beta"].to_numpy(dtype=float)) & np.isfinite(params["gamma"].to_numpy(dtype=float))
    if filter_key == "beta_positive":
        keep &= params["beta"].to_numpy(dtype=float) > 0
    params = params.loc[keep].copy()
    params["z_beta"] = zscore(params["beta"].to_numpy(dtype=float))
    params["z_gamma"] = zscore(params["gamma"].to_numpy(dtype=float))

    merged = params.merge(expression, on="region_base", how="inner")
    meta_cols = set(params.columns)
    gene_cols = [c for c in merged.columns if c not in meta_cols]
    coefs = regress_gene_coefficients(merged, gene_cols)

    pca = PCA(n_components=2).fit(coefs[["coef_beta", "coef_gamma"]].to_numpy(dtype=float))
    pc1 = pca.components_[0].copy()
    if pc1[1] < 0:
        pc1 *= -1
    pc2 = np.array([-pc1[1], pc1[0]])

    zmat = merged[["z_beta", "z_gamma"]].to_numpy(dtype=float)
    merged["eta"] = zmat @ pc1
    coefs["pc1_score"] = coefs[["coef_beta", "coef_gamma"]].to_numpy(dtype=float) @ pc1

    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    params.to_csv(results_dir / "filtered_parameters.csv", index=False)
    merged[["region", "region_base", "hemi", "beta", "gamma", "z_beta", "z_gamma", "eta"]].to_csv(
        results_dir / "region_axis.csv", index=False
    )
    coefs.to_csv(results_dir / "gene_parameter_coefficients.csv", index=False)
    gene_correlations(merged, gene_cols, "eta", "r").to_csv(results_dir / "gene_eta_correlations.csv", index=False)
    gene_correlations(merged, gene_cols, "z_beta", "r_beta").to_csv(
        results_dir / "gene_zbeta_correlations.csv", index=False
    )
    gene_correlations(merged, gene_cols, "z_gamma", "r_gamma").to_csv(
        results_dir / "gene_zgamma_correlations.csv", index=False
    )

    pd.DataFrame(
        [
            {
                "component": "PC1",
                "loading_beta": pc1[0],
                "loading_gamma": pc1[1],
                "explained_variance_ratio": pca.explained_variance_ratio_[0],
                "n_regions": len(merged),
                "n_genes": len(coefs),
            },
            {
                "component": "PC2",
                "loading_beta": pc2[0],
                "loading_gamma": pc2[1],
                "explained_variance_ratio": pca.explained_variance_ratio_[1],
                "n_regions": len(merged),
                "n_genes": len(coefs),
            },
        ]
    ).to_csv(results_dir / "pca_summary.csv", index=False)

    plot_gene_cloud(results_dir, figures_dir, spec.label)
    plot_region_axis(results_dir, figures_dir, spec.label)
    return {
        "dataset": spec.key,
        "filter": filter_key,
        "pc1_beta": pc1[0],
        "pc1_gamma": pc1[1],
        "pc1_explained": float(pca.explained_variance_ratio_[0]),
        "n_regions": int(len(merged)),
        "n_genes": int(len(coefs)),
    }


def clr_transform(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    x = np.asarray(x, dtype=float) + eps
    logx = np.log(x)
    return logx - logx.mean(axis=1, keepdims=True)


def cluster_permutation_p(x: np.ndarray, y: np.ndarray, clusters: np.ndarray, rho_obs: float, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    unique = pd.Index(pd.unique(clusters))
    first = np.array([np.where(clusters == c)[0][0] for c in unique])
    x_unique = x[first]
    exceed = 0
    for _ in range(n_perm):
        perm_values = rng.permutation(x_unique)
        lookup = dict(zip(unique, perm_values))
        xp = np.array([lookup[c] for c in clusters])
        rho, _ = spearmanr(xp, y)
        if abs(rho) >= abs(rho_obs):
            exceed += 1
    return (exceed + 1) / (n_perm + 1)


def run_cell_types(axis_path: Path, cell_path: Path, out_dir: Path, fig_dir: Path, label: str, seed: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    axis = pd.read_csv(axis_path)
    cells = pd.read_csv(cell_path)
    frac_cols = [c for c in cells.columns if c.startswith("frac_") and c != "frac_row_sum"]
    cells = cells.drop_duplicates("region_base")[["region_base"] + frac_cols]
    df = axis.merge(cells, on="region_base", how="inner")

    clr = clr_transform(df[frac_cols].to_numpy())
    clusters = df["region_base"].astype(str).to_numpy()
    y = df["eta"].to_numpy(dtype=float)
    rows = []
    for i, col in enumerate(frac_cols):
        x = clr[:, i]
        rho, p_scipy = spearmanr(x, y)
        p_perm = cluster_permutation_p(x, y, clusters, float(rho), 10_000, seed + i)
        rows.append((col, rho, p_scipy, p_perm))
    table = pd.DataFrame(rows, columns=["cell_type", "spearman_rho", "p_scipy", "p_perm"])
    table["p_fdr"] = bh_fdr(table["p_perm"].to_numpy())
    table.to_csv(out_dir / "eta_celltype_correlations.csv", index=False)

    ma_cols = [c for c in MA_CELLTYPES if c in frac_cols]
    ma_idx = [frac_cols.index(c) for c in ma_cols]
    ma_score = clr[:, ma_idx].mean(axis=1)
    rho, p_scipy = spearmanr(ma_score, y)
    p_perm = cluster_permutation_p(ma_score, y, clusters, float(rho), 10_000, seed + 10_000)
    mono = pd.DataFrame(
        [
            {
                "score": "monoaminergic",
                "components": ";".join(ma_cols),
                "spearman_rho": rho,
                "p_scipy": p_scipy,
                "p_perm": p_perm,
                "n": len(df),
            }
        ]
    )
    mono.to_csv(out_dir / "eta_monoaminergic_stats.csv", index=False)
    joint = df[["region", "region_base", "hemi", "eta"]].copy()
    joint["monoaminergic_score"] = ma_score
    joint.to_csv(out_dir / "eta_celltype_joint_table.csv", index=False)
    plot_celltype_bar(table, fig_dir)
    plot_monoaminergic(joint, mono.iloc[0], fig_dir, label)


def run_gsea(input_path: Path, out_dir: Path, fig_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_path).dropna(subset=["gene", "r"])
    rank = df[["gene", "r"]].sort_values("r", ascending=False)
    rank_path = out_dir / "ranked_genes.rnk"
    rank.to_csv(rank_path, sep="\t", index=False, header=False)
    pre = gp.prerank(
        rnk=str(rank_path),
        gene_sets="KEGG_2019_Mouse",
        organism="Mouse",
        permutation_num=1000,
        min_size=10,
        max_size=500,
        seed=42,
        outdir=str(out_dir / "gseapy"),
        verbose=False,
    )
    res = pre.res2d.copy()
    res.to_csv(out_dir / "gsea_results_all.tsv", sep="\t", float_format="%.16e")
    sig = res[pd.to_numeric(res["FDR q-val"], errors="coerce") <= 0.05]
    sig.to_csv(out_dir / "gsea_results_significant.tsv", sep="\t", float_format="%.16e")
    plot_gsea_dot(out_dir / "gsea_results_all.tsv", fig_dir)


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 150,
        }
    )


def plot_gene_cloud(result_dir: Path, fig_dir: Path, label: str) -> None:
    coefs = pd.read_csv(result_dir / "gene_parameter_coefficients.csv")
    corr = pd.read_csv(result_dir / "gene_eta_correlations.csv")[["gene", "r"]]
    pca = pd.read_csv(result_dir / "pca_summary.csv").iloc[0]
    df = coefs.merge(corr, on="gene", how="left")
    x = df["coef_beta"].to_numpy(dtype=float)
    y = df["coef_gamma"].to_numpy(dtype=float)
    lim = float(np.nanquantile(np.abs(np.concatenate([x, y])), 0.995))
    lim = max(lim, 0.04)
    fig, ax = plt.subplots(figsize=(4.8, 4.3))
    sc = ax.scatter(x, y, c=df["r"], cmap="coolwarm", s=5, alpha=0.65, rasterized=True)
    ax.axhline(0, color="0.75", lw=0.8)
    ax.axvline(0, color="0.75", lw=0.8)
    vec = np.array([pca["loading_beta"], pca["loading_gamma"]], dtype=float)
    ax.arrow(0, 0, vec[0] * lim * 0.82, vec[1] * lim * 0.82, color="black", lw=2.2, head_width=lim * 0.045, length_includes_head=True)
    ax.text(0.03, 0.97, f"PC1 {100*pca['explained_variance_ratio']:.1f}%\nn={int(pca['n_regions'])}", transform=ax.transAxes, ha="left", va="top", fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 3})
    ax.set_xlabel("coefficient for z(beta)")
    ax.set_ylabel("coefficient for z(gamma)")
    ax.set_title(f"{label}: gene coefficient PCA")
    fig.colorbar(sc, ax=ax, pad=0.02).set_label("corr(gene, eta)")
    save(fig, fig_dir / "01_gene_coefficient_cloud")


def plot_region_axis(result_dir: Path, fig_dir: Path, label: str) -> None:
    df = pd.read_csv(result_dir / "region_axis.csv")
    pca = pd.read_csv(result_dir / "pca_summary.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(4.4, 4.1))
    sc = ax.scatter(df["z_beta"], df["z_gamma"], c=df["eta"], cmap="viridis", s=32, alpha=0.85, edgecolor="white", linewidth=0.3)
    ax.axhline(0, color="0.75", lw=0.8)
    ax.axvline(0, color="0.75", lw=0.8)
    ax.set_xlabel("z(beta)")
    ax.set_ylabel("z(gamma)")
    ax.set_title(f"{label}: regional vulnerability axis")
    ax.text(0.03, 0.97, f"eta={pca['loading_beta']:.2f} z(beta)+{pca['loading_gamma']:.2f} z(gamma)", transform=ax.transAxes, ha="left", va="top", fontsize=8.5, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 3})
    fig.colorbar(sc, ax=ax, pad=0.02).set_label("eta")
    save(fig, fig_dir / "02_regional_beta_gamma_eta")


def plot_celltype_bar(table: pd.DataFrame, fig_dir: Path) -> None:
    order = ["frac_Dopa", "frac_Nora", "frac_Sero", "frac_Hist", "frac_Chol", "frac_GABA-Glyc", "frac_Glut", "frac_GABA", "frac_Glut-GABA", "frac_Unknown"]
    plot_df = table.set_index("cell_type").reindex([x for x in order if x in set(table["cell_type"])]).reset_index()
    labels = plot_df["cell_type"].str.replace("frac_", "", regex=False)
    colors = ["#c2410c" if x in MA_CELLTYPES else "#4b5563" for x in plot_df["cell_type"]]
    fig, ax = plt.subplots(figsize=(5.3, 3.4))
    bars = ax.bar(labels, plot_df["spearman_rho"], color=colors, alpha=0.88)
    for bar, p in zip(bars, plot_df["p_fdr"]):
        if p < 0.05:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + (0.025 if bar.get_height() >= 0 else -0.045), "*", ha="center", va="center", fontsize=12)
    ax.axhline(0, color="0.35", lw=0.8)
    ax.set_ylabel("Spearman rho")
    ax.set_title("Cell-type association with eta", pad=18)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    save(fig, fig_dir / "03_cell_type_correlations")


def plot_monoaminergic(joint: pd.DataFrame, stats_row: pd.Series, fig_dir: Path, label: str) -> None:
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    ax.scatter(joint["monoaminergic_score"], joint["eta"], s=30, color="#c2410c", alpha=0.75, edgecolor="white", linewidth=0.35)
    x = joint["monoaminergic_score"].to_numpy(dtype=float)
    y = joint["eta"].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    ax.text(0.03, 0.97, f"rho={stats_row['spearman_rho']:.2f}\nperm p={stats_row['p_perm']:.3g}\nn={int(stats_row['n'])}", transform=ax.transAxes, ha="left", va="top", fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 3})
    ax.set_xlabel("monoaminergic score")
    ax.set_ylabel("eta")
    save(fig, fig_dir / "04_monoaminergic_eta")


def plot_gsea_dot(gsea_path: Path, fig_dir: Path) -> None:
    df = pd.read_csv(gsea_path, sep="\t")
    df = df.rename(columns={"FDR q-val": "FDR"})
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df[df["FDR"] <= 0.05].copy()
    if df.empty:
        return
    df["absNES"] = df["NES"].abs()
    df = df.sort_values("absNES", ascending=False).head(14)
    labels = ["\n".join(textwrap.wrap(x, width=30)) for x in df["Term"]]
    weights = -np.log10(df["FDR"].clip(lower=1e-300))
    sizes = 30 + 110 * weights / max(weights.max(), 1)
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(5.8, max(3.2, 0.32 * len(df) + 1.2)))
    ax.scatter(df["NES"], y, s=sizes, color="black")
    ax.axvline(0, color="0.55", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("NES")
    ax.set_title("Top significant KEGG GSEA terms")
    save(fig, fig_dir / "05_gsea_top_terms")


def compare_outputs(filter_key: str, results_root: Path, figures_root: Path) -> dict[str, float | str]:
    sdir = results_root / filter_key / "striatum" / "transcriptomics"
    hdir = results_root / filter_key / "hippocampus_C3_C4" / "transcriptomics"
    out = figures_root / filter_key / "comparisons"
    out.mkdir(parents=True, exist_ok=True)

    sp = pd.read_csv(sdir / "pca_summary.csv").iloc[0]
    hp = pd.read_csv(hdir / "pca_summary.csv").iloc[0]
    sv = np.array([sp["loading_beta"], sp["loading_gamma"]], dtype=float)
    hv = np.array([hp["loading_beta"], hp["loading_gamma"]], dtype=float)
    cosine = float(np.dot(sv, hv) / (np.linalg.norm(sv) * np.linalg.norm(hv)))
    angle = float(np.degrees(np.arccos(np.clip(abs(cosine), -1, 1))))

    fig, ax = plt.subplots(figsize=(4.5, 4.3))
    ax.axhline(0, color="0.82", lw=0.8)
    ax.axvline(0, color="0.82", lw=0.8)
    for label, vec, color in [("Striatum", sv, "#2b6cb0"), ("Hippocampus", hv, "#c2410c")]:
        ax.arrow(0, 0, vec[0], vec[1], color=color, lw=2.5, head_width=0.04, length_includes_head=True)
        ax.text(vec[0] * 1.12, vec[1] * 1.12, label, color=color, ha="center", va="center", fontsize=9)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel("PC1 loading beta")
    ax.set_ylabel("PC1 loading gamma")
    ax.set_title(f"PC1 direction comparison\nangle={angle:.1f} deg")
    save(fig, out / "01_pc1_direction_comparison")

    eta_r = scatter_compare(
        sdir / "region_axis.csv",
        hdir / "region_axis.csv",
        "eta",
        "eta",
        "region_base",
        "Striatum eta",
        "Hippocampus eta",
        out / "02_region_eta_correlation",
    )
    beta_r = scatter_compare(
        sdir / "gene_zbeta_correlations.csv",
        hdir / "gene_zbeta_correlations.csv",
        "r_beta",
        "r_beta",
        "gene",
        "Striatum corr(gene, z(beta))",
        "Hippocampus corr(gene, z(beta))",
        out / "03_gene_zbeta_correlation_comparison",
    )
    eta_gene_r = scatter_compare(
        sdir / "gene_eta_correlations.csv",
        hdir / "gene_eta_correlations.csv",
        "r",
        "r",
        "gene",
        "Striatum corr(gene, eta)",
        "Hippocampus corr(gene, eta)",
        out / "04_gene_eta_correlation_comparison",
    )
    return {"filter": filter_key, "pc1_angle_deg": angle, "region_eta_r": eta_r, "gene_zbeta_r": beta_r, "gene_eta_r": eta_gene_r}


def scatter_compare(
    left_path: Path,
    right_path: Path,
    left_col: str,
    right_col: str,
    key: str,
    xlabel: str,
    ylabel: str,
    out_path: Path,
) -> float:
    left = pd.read_csv(left_path)[[key, left_col]].rename(columns={left_col: "x"})
    right = pd.read_csv(right_path)[[key, right_col]].rename(columns={right_col: "y"})
    if key == "region_base":
        left = left.groupby(key, as_index=False)["x"].mean()
        right = right.groupby(key, as_index=False)["y"].mean()
    df = left.merge(right, on=key, how="inner").dropna()
    r, p = stats.pearsonr(df["x"], df["y"]) if len(df) >= 3 else (np.nan, np.nan)
    if len(df) == 0:
        lim = 1.0
    else:
        lim = float(np.nanquantile(np.abs(df[["x", "y"]].to_numpy(dtype=float)), 0.995))
        lim = max(lim, 0.05)
    fig, ax = plt.subplots(figsize=(4.3, 4.1))
    ax.scatter(df["x"], df["y"], s=8, color="0.25", alpha=0.45, rasterized=True)
    ax.plot([-lim, lim], [-lim, lim], color="0.65", lw=1.0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.text(0.04, 0.96, f"r={r:.2f}\np={p:.2g}\nn={len(df)}", transform=ax.transAxes, ha="left", va="top", fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 3})
    save(fig, out_path)
    df.to_csv(out_path.with_suffix(".csv"), index=False)
    return float(r)


def combined_gsea(filter_key: str, results_root: Path, figures_root: Path) -> None:
    s = pd.read_csv(results_root / filter_key / "striatum" / "enrichment" / "gsea_results_all.tsv", sep="\t").assign(dataset="Striatum")
    h = pd.read_csv(results_root / filter_key / "hippocampus_C3_C4" / "enrichment" / "gsea_results_all.tsv", sep="\t").assign(dataset="Hippocampus")
    df = pd.concat([s, h], ignore_index=True).rename(columns={"FDR q-val": "FDR"})
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df["absNES"] = df["NES"].abs()
    term_stats = df.groupby("Term").agg(sig_count=("FDR", lambda x: int((x <= 0.05).sum())), mean_abs=("absNES", "mean"), max_abs=("absNES", "max")).sort_values(["sig_count", "mean_abs", "max_abs"], ascending=False)
    terms = term_stats.head(18).index.tolist()
    plot = df[df["Term"].isin(terms)].copy()
    order = plot.groupby("Term")["absNES"].mean().sort_values(ascending=True).index.tolist()
    ymap = {t: i for i, t in enumerate(order)}
    fig, ax = plt.subplots(figsize=(8.4, max(4.4, 0.34 * len(order) + 1.3)))
    colors = {"Striatum": "#2b6cb0", "Hippocampus": "#c2410c"}
    offsets = {"Striatum": -0.16, "Hippocampus": 0.16}
    for dataset in ["Striatum", "Hippocampus"]:
        sub = plot[plot["dataset"] == dataset]
        y = np.array([ymap[t] + offsets[dataset] for t in sub["Term"]])
        weights = -np.log10(sub["FDR"].clip(lower=1e-300))
        sizes = 28 + 95 * weights / max(weights.max(), 1)
        ax.scatter(sub["NES"], y, s=sizes, color=colors[dataset], alpha=0.88, label=dataset, edgecolor="white", linewidth=0.35)
    ax.axvline(0, color="0.55", lw=0.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlabel("NES")
    ax.set_title(f"{filter_key}: striatum vs hippocampus KEGG GSEA")
    ax.legend(frameon=False, loc="lower right")
    out = figures_root / filter_key / "comparisons"
    save(fig, out / "05_gsea_striatum_vs_hippocampus")


def combined_panel_contact_sheets(filter_key: str, figures_root: Path) -> None:
    out = figures_root / filter_key / "comparisons"
    out.mkdir(parents=True, exist_ok=True)
    panels = [
        ("Gene coefficient cloud", "transcriptomics/01_gene_coefficient_cloud.png"),
        ("Regional eta", "transcriptomics/02_regional_beta_gamma_eta.png"),
        ("Cell types", "cell_types/03_cell_type_correlations.png"),
        ("Monoaminergic eta", "cell_types/04_monoaminergic_eta.png"),
        ("GSEA", "enrichment/05_gsea_top_terms.png"),
    ]
    datasets = [("Striatum", "striatum"), ("Hippocampus", "hippocampus_C3_C4")]
    for title, rel in panels:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
        for ax, (dataset_label, dataset_key) in zip(axes, datasets):
            img_path = figures_root / filter_key / dataset_key / rel
            ax.imshow(mpimg.imread(img_path))
            ax.set_title(dataset_label, fontsize=11)
            ax.axis("off")
        fig.suptitle(f"{filter_key}: {title}", y=0.98, fontsize=12)
        fig.tight_layout()
        slug = title.lower().replace(" ", "_")
        save(fig, out / f"06_side_by_side_{slug}")


def diagnostics_summary(results_root: Path) -> None:
    rows = []
    for spec in DATASETS:
        chain_fit = spec.run_dir / "plots" / "diagnostics" / "chain_fit_metrics.csv"
        rhat = spec.run_dir / "plots" / "diagnostics" / "rhat_summary.csv"
        if rhat.exists():
            rh = pd.read_csv(rhat)
            max_rhat = float(rh["rhat"].max())
            median_rhat = float(rh["rhat"].median())
            n_gt = int((rh["rhat"] > 1.05).sum())
            n_params = int(len(rh))
        else:
            max_rhat = median_rhat = np.nan
            n_gt = n_params = 0
        if chain_fit.exists():
            cf = pd.read_csv(chain_fit)
            loglik = "; ".join(f"C{int(c)}={v:.1f}" for c, v in zip(cf["chain"], cf["loglik_all"]))
        else:
            loglik = ""
        rows.append({"dataset": spec.key, "run_dir": str(spec.run_dir), "max_rhat": max_rhat, "median_rhat": median_rhat, "n_rhat_gt_1p05": n_gt, "n_parameters": n_params, "chain_logliks": loglik})
    pd.DataFrame(rows).to_csv(results_root / "diagnostics_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", default="paper-rf/results/figure_6_7")
    parser.add_argument("--figures-root", default="paper-rf/figures/figure_6_7")
    parser.add_argument("--expression", default="paper-rf/data/transcriptomics/avg_Pangea_exp.csv")
    parser.add_argument("--cell-types", default="paper-rf/data/cell_types/connectome_celltype.csv")
    args = parser.parse_args()

    setup_style()
    results_root = Path(args.results_root)
    figures_root = Path(args.figures_root)
    expression = load_expression(Path(args.expression))
    cell_types = Path(args.cell_types)

    param_dirs = {}
    for spec in DATASETS:
        param_dir = results_root / "parameters" / spec.key
        export_parameter_tables(spec, param_dir)
        param_dirs[spec.key] = param_dir
    params = {key: load_parameters(path) for key, path in param_dirs.items()}
    rows = []
    comparison_rows = []
    for filter_key in ["all", "beta_positive"]:
        for spec in DATASETS:
            base_results = results_root / filter_key / spec.key
            base_figures = figures_root / filter_key / spec.key
            row = run_transcriptomics(
                spec,
                filter_key,
                params[spec.key],
                expression,
                base_results / "transcriptomics",
                base_figures / "transcriptomics",
            )
            rows.append(row)
            run_cell_types(
                base_results / "transcriptomics" / "region_axis.csv",
                cell_types,
                base_results / "cell_types",
                base_figures / "cell_types",
                spec.label,
                seed=0 if spec.key == "striatum" else 100,
            )
            run_gsea(
                base_results / "transcriptomics" / "gene_eta_correlations.csv",
                base_results / "enrichment",
                base_figures / "enrichment",
            )
        comparison_rows.append(compare_outputs(filter_key, results_root, figures_root))
        combined_gsea(filter_key, results_root, figures_root)
        combined_panel_contact_sheets(filter_key, figures_root)

    results_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(results_root / "pca_summary_by_dataset_filter.csv", index=False)
    pd.DataFrame(comparison_rows).to_csv(results_root / "comparison_summary.csv", index=False)
    diagnostics_summary(results_root)
    print(figures_root)


if __name__ == "__main__":
    main()
