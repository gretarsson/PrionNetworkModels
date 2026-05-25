#!/usr/bin/env python3
"""Pre-ranked GSEA for gene-axis correlations."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="gene_eta_correlations.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--gene-sets", default="KEGG_2019_Mouse")
    parser.add_argument("--min-size", type=int, default=10)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input).dropna(subset=["gene", "r"])
    rank = df[["gene", "r"]].sort_values("r", ascending=False)
    rank_path = out_dir / "ranked_genes.rnk"
    rank.to_csv(rank_path, sep="\t", index=False, header=False)

    pre = gp.prerank(
        rnk=str(rank_path),
        gene_sets=args.gene_sets,
        organism="Mouse",
        permutation_num=args.permutations,
        min_size=args.min_size,
        max_size=args.max_size,
        seed=args.seed,
        outdir=str(out_dir / "gseapy"),
        verbose=True,
    )
    res = pre.res2d.copy()
    res.to_csv(out_dir / "gsea_results_all.tsv", sep="\t", float_format="%.16e")
    if "FDR q-val" in res.columns:
        sig = res[pd.to_numeric(res["FDR q-val"], errors="coerce") <= 0.05]
        sig.to_csv(out_dir / "gsea_results_significant.tsv", sep="\t", float_format="%.16e")
    print(out_dir)


if __name__ == "__main__":
    main()
