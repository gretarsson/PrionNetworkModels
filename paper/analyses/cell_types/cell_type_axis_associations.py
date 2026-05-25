#!/usr/bin/env python3
"""Cell-type associations with the PCA vulnerability axis eta."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


MA_CELLTYPES = ["frac_Dopa", "frac_Nora", "frac_Sero", "frac_Hist"]


def clr_transform(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    x = np.asarray(x, dtype=float) + eps
    logx = np.log(x)
    return logx - logx.mean(axis=1, keepdims=True)


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--axis", required=True, help="region_axis.csv from gene_parameter_pca.py")
    parser.add_argument("--cell-types", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-perm", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    axis = pd.read_csv(args.axis)
    cells = pd.read_csv(args.cell_types)
    if "region_base" not in cells.columns:
        raise SystemExit("Cell-type table must contain region_base")
    frac_cols = [c for c in cells.columns if c.startswith("frac_") and c != "frac_row_sum"]
    cells = cells.drop_duplicates("region_base")[["region_base"] + frac_cols]
    df = axis.merge(cells, on="region_base", how="inner")

    clr = clr_transform(df[frac_cols].to_numpy())
    rows = []
    clusters = df["region_base"].astype(str).to_numpy()
    y = df["eta"].to_numpy(dtype=float)
    for i, col in enumerate(frac_cols):
        x = clr[:, i]
        rho, p_scipy = spearmanr(x, y)
        p_perm = cluster_permutation_p(x, y, clusters, float(rho), args.n_perm, args.seed + i)
        rows.append((col, rho, p_scipy, p_perm))

    table = pd.DataFrame(rows, columns=["cell_type", "spearman_rho", "p_scipy", "p_perm"])
    table["p_fdr"] = bh_fdr(table["p_perm"].to_numpy())
    table.to_csv(out_dir / "eta_celltype_correlations.csv", index=False)

    ma_cols = [c for c in MA_CELLTYPES if c in frac_cols]
    ma_idx = [frac_cols.index(c) for c in ma_cols]
    ma_score = clr[:, ma_idx].mean(axis=1)
    rho, p_scipy = spearmanr(ma_score, y)
    p_perm = cluster_permutation_p(ma_score, y, clusters, float(rho), args.n_perm, args.seed + 10_000)
    pd.DataFrame([{
        "score": "monoaminergic",
        "components": ";".join(ma_cols),
        "spearman_rho": rho,
        "p_scipy": p_scipy,
        "p_perm": p_perm,
        "n": len(df),
    }]).to_csv(out_dir / "eta_monoaminergic_stats.csv", index=False)

    joint = df[["region", "region_base", "hemi", "eta"]].copy()
    joint["monoaminergic_score"] = ma_score
    joint.to_csv(out_dir / "eta_celltype_joint_table.csv", index=False)
    print(out_dir)


if __name__ == "__main__":
    main()
