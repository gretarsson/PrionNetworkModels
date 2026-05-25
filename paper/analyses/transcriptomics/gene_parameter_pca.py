#!/usr/bin/env python3
"""Gene-parameter PCA vulnerability-axis analysis.

This script ports the manuscript-critical gene analysis into a reproducible,
argument-driven workflow:

1. align regional gene expression with hemisphere-specific beta/gamma maps,
2. filter to beta-positive and posterior-updated regions,
3. fit expression_g ~ z(beta) + z(gamma) for every gene,
4. run PCA on the gene-level coefficient pairs,
5. define eta = PC1_beta * z(beta) + PC1_gamma * z(gamma),
6. rank genes by Pearson correlation with eta.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA


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


def load_parameter_pair(
    beta_path: Path,
    gamma_path: Path,
    update_alpha: float,
    beta_min: float,
    filter_updated: bool,
) -> pd.DataFrame:
    beta = pd.read_csv(beta_path)
    gamma = pd.read_csv(gamma_path)
    beta = beta.rename(columns={"mean_post": "beta", "ks_pvalue": "beta_ks_pvalue", "updated": "beta_updated"})
    gamma = gamma.rename(columns={"mean_post": "gamma", "ks_pvalue": "gamma_ks_pvalue", "updated": "gamma_updated"})
    keep_beta_cols = ["region", "beta", "beta_ks_pvalue", "beta_updated"]
    keep_gamma_cols = ["region", "gamma", "gamma_ks_pvalue", "gamma_updated"]
    params = beta[keep_beta_cols].merge(gamma[keep_gamma_cols], on="region", how="inner")
    params["region_base"] = params["region"].map(strip_hemi)
    params["hemi"] = params["region"].astype(str).str[0]
    params = params[(params["beta"] > beta_min)].copy()
    if filter_updated:
        params = params[(params["beta_ks_pvalue"] < update_alpha) & (params["gamma_ks_pvalue"] < update_alpha)].copy()
    params["z_beta"] = zscore(params["beta"].to_numpy())
    params["z_gamma"] = zscore(params["gamma"].to_numpy())
    return params


def load_expression(expression_path: Path) -> pd.DataFrame:
    expr = pd.read_csv(expression_path)
    region_col = expr.columns[0]
    expr = expr.rename(columns={region_col: "region_base"})
    expr["region_base"] = expr["region_base"].astype(str)
    return expr


def align_expression(params: pd.DataFrame, expression: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    merged = params.merge(expression, on="region_base", how="inner")
    meta_cols = {
        "region", "region_base", "hemi", "beta", "gamma", "beta_ks_pvalue",
        "gamma_ks_pvalue", "beta_updated", "gamma_updated", "z_beta", "z_gamma",
    }
    gene_cols = [c for c in merged.columns if c not in meta_cols]
    return merged, gene_cols


def regress_gene_coefficients(df: pd.DataFrame, gene_cols: list[str]) -> pd.DataFrame:
    x = np.column_stack([df["z_beta"].to_numpy(), df["z_gamma"].to_numpy()])
    # No intercept after centering expression; coefficients are comparable to standardized parameter axes.
    rows = []
    for gene in gene_cols:
        y = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if mask.sum() < 5:
            continue
        yz = zscore(y[mask])
        coef, *_ = np.linalg.lstsq(x[mask], yz, rcond=None)
        rows.append((gene, coef[0], coef[1], int(mask.sum())))
    return pd.DataFrame(rows, columns=["gene", "coef_beta", "coef_gamma", "n_used"])


def orient_pc1(loadings: np.ndarray) -> np.ndarray:
    # Manuscript convention: increasing eta tracks stronger fall dynamics.
    return -loadings if loadings[1] < 0 else loadings


def gene_correlations(df: pd.DataFrame, gene_cols: list[str], eta: np.ndarray) -> pd.DataFrame:
    rows = []
    for gene in gene_cols:
        y = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y) & np.isfinite(eta)
        if mask.sum() < 5:
            continue
        r, p = stats.pearsonr(y[mask], eta[mask])
        rows.append((gene, r, p, int(mask.sum())))
    out = pd.DataFrame(rows, columns=["gene", "r", "p_un", "n_used"])
    out["p_bonf"] = np.minimum(out["p_un"] * len(out), 1.0)
    out["p_fdr"] = bh_fdr(out["p_un"].to_numpy())
    return out.sort_values("r", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True)
    parser.add_argument("--beta", required=True)
    parser.add_argument("--gamma", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--update-alpha", type=float, default=0.001)
    parser.add_argument("--beta-min", type=float, default=0.0)
    parser.add_argument("--no-update-filter", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = load_parameter_pair(
        Path(args.beta),
        Path(args.gamma),
        args.update_alpha,
        args.beta_min,
        filter_updated=not args.no_update_filter,
    )
    expr = load_expression(Path(args.expression))
    merged, gene_cols = align_expression(params, expr)
    coefs = regress_gene_coefficients(merged, gene_cols)

    pca = PCA(n_components=2)
    pca.fit(coefs[["coef_beta", "coef_gamma"]].to_numpy())
    pc1 = orient_pc1(pca.components_[0].copy())
    eta = pc1[0] * merged["z_beta"].to_numpy() + pc1[1] * merged["z_gamma"].to_numpy()
    merged["eta"] = eta

    corrs = gene_correlations(merged, gene_cols, eta)
    coefs["pc1_score"] = coefs[["coef_beta", "coef_gamma"]].to_numpy() @ pc1

    params.to_csv(out_dir / "filtered_parameters.csv", index=False)
    merged[["region", "region_base", "hemi", "beta", "gamma", "z_beta", "z_gamma", "eta"]].to_csv(
        out_dir / "region_axis.csv", index=False
    )
    coefs.to_csv(out_dir / "gene_parameter_coefficients.csv", index=False)
    corrs.to_csv(out_dir / "gene_eta_correlations.csv", index=False)
    pd.DataFrame({
        "component": ["PC1"],
        "loading_beta": [pc1[0]],
        "loading_gamma": [pc1[1]],
        "explained_variance_ratio": [pca.explained_variance_ratio_[0]],
        "n_regions": [len(merged)],
        "n_genes": [len(coefs)],
    }).to_csv(out_dir / "pca_summary.csv", index=False)

    print(out_dir)


if __name__ == "__main__":
    main()
