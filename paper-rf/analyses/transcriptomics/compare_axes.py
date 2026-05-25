#!/usr/bin/env python3
"""Compare striatum and hippocampus gene-axis outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def read_pc(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    return df.loc[0, ["loading_beta", "loading_gamma"]].to_numpy(dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--striatum-dir", required=True)
    parser.add_argument("--hippocampus-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    sdir = Path(args.striatum_dir)
    hdir = Path(args.hippocampus_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    spc = read_pc(sdir / "pca_summary.csv")
    hpc = read_pc(hdir / "pca_summary.csv")
    cosine = float(np.dot(spc, hpc) / (np.linalg.norm(spc) * np.linalg.norm(hpc)))
    angle = float(np.degrees(np.arccos(np.clip(abs(cosine), -1, 1))))

    sg = pd.read_csv(sdir / "gene_eta_correlations.csv")[["gene", "r"]].rename(columns={"r": "r_striatum"})
    hg = pd.read_csv(hdir / "gene_eta_correlations.csv")[["gene", "r"]].rename(columns={"r": "r_hippocampus"})
    merged = sg.merge(hg, on="gene", how="inner").dropna()
    r, p = stats.pearsonr(merged["r_striatum"], merged["r_hippocampus"])
    slope, intercept, *_ = stats.linregress(merged["r_striatum"], merged["r_hippocampus"])
    sign_agreement = float((np.sign(merged["r_striatum"]) == np.sign(merged["r_hippocampus"])).mean())

    merged.to_csv(out / "gene_eta_comparison.csv", index=False)
    pd.DataFrame([{
        "pc_cosine_similarity": cosine,
        "pc_absolute_angle_degrees": angle,
        "gene_eta_pearson_r": r,
        "gene_eta_pearson_p": p,
        "gene_eta_ols_slope": slope,
        "gene_eta_ols_intercept": intercept,
        "gene_eta_sign_agreement": sign_agreement,
        "n_genes": len(merged),
    }]).to_csv(out / "axis_comparison_summary.csv", index=False)
    print(out)


if __name__ == "__main__":
    main()
