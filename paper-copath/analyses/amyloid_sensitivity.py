#!/usr/bin/env python3
"""Amyloid sensitivity analyses for co-pathology REGION-RF fits."""

from __future__ import annotations

import argparse
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA


PARAMETERS = ["alpha", "beta", "gamma"]
AMYLOIDS = ["ab40", "ab42"]
MA_CELLTYPES = ["frac_Dopa", "frac_Nora", "frac_Sero", "frac_Hist"]
PROTEIN_LABEL = {"syn": "Synuclein", "tau": "Tau"}
PROTEIN_COLOR = {"syn": "#0047AB", "tau": "#C43616"}


def strip_hemi(region: str) -> str:
    region = str(region).strip()
    return region[1:] if len(region) > 1 and region[0] in {"i", "c"} else region


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values, dtype=float)
    return (values - np.nanmean(values)) / sd


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    out = np.full(n, np.nan)
    finite = np.isfinite(pvalues)
    if finite.sum() == 0:
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(pvalues[finite])]
    ranked = pvalues[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig: plt.Figure, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def truncate_label(label: str, max_len: int = 40) -> str:
    label = str(label)
    return label if len(label) <= max_len else f"{label[: max_len - 1]}..."


def clr_transform(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    x = np.asarray(x, dtype=float) + eps
    logx = np.log(x)
    return logx - logx.mean(axis=1, keepdims=True)


def corr_stats(x, y) -> dict[str, float | int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return {"n": int(mask.sum()), "pearson_r": np.nan, "pearson_p": np.nan, "spearman_r": np.nan, "spearman_p": np.nan}
    pr = stats.pearsonr(x[mask], y[mask])
    sr = stats.spearmanr(x[mask], y[mask])
    return {
        "n": int(mask.sum()),
        "pearson_r": float(pr.statistic),
        "pearson_p": float(pr.pvalue),
        "spearman_r": float(sr.statistic),
        "spearman_p": float(sr.pvalue),
    }


def cluster_permutation_p(x: np.ndarray, y: np.ndarray, clusters: np.ndarray, rho_obs: float, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    clusters = clusters[mask]
    if len(x) < 3:
        return np.nan
    unique = pd.Index(pd.unique(clusters))
    first = np.array([np.where(clusters == c)[0][0] for c in unique])
    x_unique = x[first]
    cluster_to_idx = {cluster: i for i, cluster in enumerate(unique)}
    inverse = np.array([cluster_to_idx[c] for c in clusters])

    perm_unique = rng.permuted(np.tile(x_unique, (n_perm, 1)), axis=1)
    permuted_x = perm_unique[:, inverse]
    rx = stats.rankdata(permuted_x, axis=1)
    ry = stats.rankdata(y)
    rx = rx - rx.mean(axis=1, keepdims=True)
    ry = ry - ry.mean()
    denom = np.sqrt(np.sum(rx**2, axis=1) * np.sum(ry**2))
    rho_perm = np.divide(rx @ ry, denom, out=np.zeros(n_perm), where=denom > 0)
    exceed = int(np.sum(np.abs(rho_perm) >= abs(rho_obs)))
    return (exceed + 1) / (n_perm + 1)


def regression_table(df: pd.DataFrame, protein: str, amyloid: str) -> pd.DataFrame:
    rows = []
    amy_col = f"{amyloid}_treatment_mean_prelimval"
    for parameter in PARAMETERS + ["delta_pc1"]:
        y_col = f"{parameter}_diff_app_minus_mapt" if parameter in PARAMETERS else parameter
        sub = df[[amy_col, y_col, "peak_mapt", "active_rhat"]].copy()
        sub = sub[sub["active_rhat"]]
        sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
        if len(sub) < 10:
            continue
        x = np.column_stack(
            [
                np.ones(len(sub)),
                zscore(sub[amy_col].to_numpy()),
                zscore(sub["peak_mapt"].to_numpy()),
            ]
        )
        y = zscore(sub[y_col].to_numpy())
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        resid = y - x @ beta
        dof = len(y) - x.shape[1]
        sigma2 = np.sum(resid**2) / dof
        cov = sigma2 * np.linalg.inv(x.T @ x)
        se = np.sqrt(np.diag(cov))
        t = beta / se
        p = 2 * stats.t.sf(np.abs(t), dof)
        rows.append(
            {
                "protein": protein,
                "amyloid": amyloid,
                "outcome": parameter,
                "n_regions": len(sub),
                "amyloid_standardized_beta": beta[1],
                "amyloid_t": t[1],
                "amyloid_p": p[1],
                "baseline_peak_standardized_beta": beta[2],
                "baseline_peak_p": p[2],
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["amyloid_p_fdr"] = bh_fdr(out["amyloid_p"].to_numpy())
    return out


def load_comparison_tables(result_dir: Path) -> dict[str, pd.DataFrame]:
    tables = {}
    for protein in ["syn", "tau"]:
        df = pd.read_csv(result_dir / f"{protein}_app_vs_mapt_region_parameters.csv")
        df["protein"] = protein
        df["region_base"] = df["region"].map(strip_hemi)
        rhat_mask = df["active_any"].astype(bool).to_numpy()
        for parameter in PARAMETERS:
            rhat_mask &= pd.to_numeric(df[f"{parameter}_rhat_app"], errors="coerce").to_numpy() <= 1.05
            rhat_mask &= pd.to_numeric(df[f"{parameter}_rhat_mapt"], errors="coerce").to_numpy() <= 1.05
        df["active_rhat"] = rhat_mask
        for parameter in PARAMETERS:
            diff_col = f"{parameter}_diff_app_minus_mapt"
            if diff_col not in df.columns:
                df[diff_col] = df[f"{parameter}_app"] - df[f"{parameter}_mapt"]
        tables[protein] = df
    return tables


def add_delta_pca(df: pd.DataFrame, protein: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = df["active_rhat"].to_numpy()
    diff_cols = [f"{p}_diff_app_minus_mapt" for p in PARAMETERS]
    x = df.loc[mask, diff_cols].to_numpy(dtype=float)
    finite = np.isfinite(x).all(axis=1)
    idx = df.index[mask][finite]
    xz = np.column_stack([zscore(x[finite, j]) for j in range(x.shape[1])])
    pca = PCA(n_components=3).fit(xz)
    scores = xz @ pca.components_.T
    ab40 = df.loc[idx, "ab40_treatment_mean_prelimval"].to_numpy(dtype=float)
    orient_mask = np.isfinite(ab40) & np.isfinite(scores[:, 0])
    if orient_mask.sum() >= 3 and stats.pearsonr(scores[orient_mask, 0], ab40[orient_mask]).statistic < 0:
        pca.components_[0, :] *= -1
        scores[:, 0] *= -1

    out = df.copy()
    out["delta_pc1"] = np.nan
    out["delta_pc2"] = np.nan
    out["delta_pc3"] = np.nan
    out.loc[idx, ["delta_pc1", "delta_pc2", "delta_pc3"]] = scores

    rows = []
    for comp_idx in range(3):
        row = {
            "protein": protein,
            "component": f"PC{comp_idx + 1}",
            "explained_variance_ratio": pca.explained_variance_ratio_[comp_idx],
        }
        for j, parameter in enumerate(PARAMETERS):
            row[f"loading_delta_{parameter}"] = pca.components_[comp_idx, j]
        rows.append(row)
    return out, pd.DataFrame(rows)


def load_observation_means(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, na_values=["NA"])
    time_col = "mpi"
    region_cols = [c for c in raw.columns if c not in {"mouse", time_col}]
    mean = raw.groupby(time_col)[region_cols].mean(numeric_only=True).sort_index()
    se = raw.groupby(time_col)[region_cols].sem(numeric_only=True).sort_index()
    return mean, se


def plot_amyloid_parameter_grid(tables: dict[str, pd.DataFrame], figure_dir: Path) -> None:
    for amyloid in AMYLOIDS:
        amy_col = f"{amyloid}_treatment_mean_prelimval"
        fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.6), sharex=False)
        for row_idx, protein in enumerate(["syn", "tau"]):
            df = tables[protein]
            mask = df["active_rhat"].to_numpy()
            for col_idx, parameter in enumerate(PARAMETERS):
                ax = axes[row_idx, col_idx]
                y_col = f"{parameter}_diff_app_minus_mapt"
                x = df.loc[mask, amy_col].to_numpy(dtype=float)
                y = df.loc[mask, y_col].to_numpy(dtype=float)
                ax.scatter(x, y, s=28, alpha=0.72, color=PROTEIN_COLOR[protein], edgecolor="white", linewidth=0.35)
                finite = np.isfinite(x) & np.isfinite(y)
                if finite.sum() >= 3:
                    slope, intercept, r, p, _ = stats.linregress(x[finite], y[finite])
                    xx = np.linspace(np.nanmin(x[finite]), np.nanmax(x[finite]), 100)
                    ax.plot(xx, intercept + slope * xx, color="black", lw=1.8)
                    ax.text(0.04, 0.96, f"r={r:.2f}, p={p:.2g}\nn={finite.sum()}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
                ax.axhline(0, color="0.75", lw=0.8)
                ax.set_title(f"{PROTEIN_LABEL[protein]}: delta {parameter}")
                ax.set_xlabel(f"{amyloid.upper()} amyloid load")
                ax.set_ylabel("APP - MAPT")
                style_axis(ax)
        fig.suptitle(f"Amyloid load vs APP-induced dynamical parameter shifts ({amyloid.upper()})", y=0.98)
        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.90, wspace=0.34, hspace=0.52)
        save_figure(fig, figure_dir / f"{amyloid}_parameter_shift_grid")


def plot_delta_pca(tables: dict[str, pd.DataFrame], pca_summary: pd.DataFrame, figure_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 8.7))
    for row_idx, protein in enumerate(["syn", "tau"]):
        df = tables[protein]
        summary = pca_summary[(pca_summary["protein"] == protein) & (pca_summary["component"] == "PC1")].iloc[0]
        ax = axes[row_idx, 0]
        loads = [summary[f"loading_delta_{p}"] for p in PARAMETERS]
        ax.bar(PARAMETERS, loads, color=PROTEIN_COLOR[protein], alpha=0.86)
        ax.axhline(0, color="0.35", lw=0.8)
        ax.set_ylim(-1, 1)
        ax.set_title(f"{PROTEIN_LABEL[protein]} delta-PC1 loadings")
        ax.set_ylabel("loading")
        style_axis(ax)

        for col_idx, amyloid in enumerate(AMYLOIDS, start=1):
            ax = axes[row_idx, col_idx]
            amy_col = f"{amyloid}_treatment_mean_prelimval"
            sub = df[df["active_rhat"]].copy()
            x = sub[amy_col].to_numpy(dtype=float)
            y = sub["delta_pc1"].to_numpy(dtype=float)
            ax.scatter(x, y, s=30, color=PROTEIN_COLOR[protein], alpha=0.72, edgecolor="white", linewidth=0.35)
            finite = np.isfinite(x) & np.isfinite(y)
            if finite.sum() >= 3:
                slope, intercept, r, p, _ = stats.linregress(x[finite], y[finite])
                xx = np.linspace(np.nanmin(x[finite]), np.nanmax(x[finite]), 100)
                ax.plot(xx, intercept + slope * xx, color="black", lw=1.8)
                ax.text(0.04, 0.96, f"r={r:.2f}, p={p:.2g}\nn={finite.sum()}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
            ax.axhline(0, color="0.75", lw=0.8)
            ax.set_title(f"delta-PC1 vs {amyloid.upper()}")
            ax.set_xlabel(f"{amyloid.upper()} amyloid load")
            ax.set_ylabel("amyloid sensitivity axis")
            style_axis(ax)
    fig.suptitle("Low-dimensional APP effect on local RF dynamics", y=0.98)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.90, wspace=0.34, hspace=0.52)
    save_figure(fig, figure_dir / "delta_dynamics_pca")


def plot_top_region_bars(tables: dict[str, pd.DataFrame], figure_dir: Path, n: int = 12) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.5))
    for ax, protein in zip(axes, ["syn", "tau"]):
        sub = tables[protein][tables[protein]["active_rhat"]].copy()
        sub = sub[np.isfinite(sub["delta_pc1"])].copy()
        top = pd.concat([sub.nsmallest(n // 2, "delta_pc1"), sub.nlargest(n // 2, "delta_pc1")]).sort_values("delta_pc1")
        colors = ["#777777" if v < 0 else PROTEIN_COLOR[protein] for v in top["delta_pc1"]]
        ax.barh(top["region"], top["delta_pc1"], color=colors, alpha=0.88)
        ax.axvline(0, color="black", lw=0.9)
        ax.set_title(f"{PROTEIN_LABEL[protein]} regions with strongest APP-effect axis scores")
        ax.set_xlabel("delta-PC1")
        style_axis(ax)
    fig.tight_layout()
    save_figure(fig, figure_dir / "top_region_delta_pc1_bars")


def plot_syn_tau_comparison(tables: dict[str, pd.DataFrame], figure_dir: Path) -> pd.DataFrame:
    syn = tables["syn"][["region", "active_rhat"] + [f"{p}_diff_app_minus_mapt" for p in PARAMETERS] + ["delta_pc1"]].copy()
    tau = tables["tau"][["region", "active_rhat"] + [f"{p}_diff_app_minus_mapt" for p in PARAMETERS] + ["delta_pc1"]].copy()
    merged = syn.merge(tau, on="region", suffixes=("_syn", "_tau"))
    merged["active_rhat_both"] = merged["active_rhat_syn"] & merged["active_rhat_tau"]
    rows = []
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.8))
    for ax, parameter in zip(axes, PARAMETERS + ["delta_pc1"]):
        x_col = f"{parameter}_diff_app_minus_mapt_syn" if parameter in PARAMETERS else "delta_pc1_syn"
        y_col = f"{parameter}_diff_app_minus_mapt_tau" if parameter in PARAMETERS else "delta_pc1_tau"
        sub = merged[merged["active_rhat_both"]]
        x = sub[x_col].to_numpy(dtype=float)
        y = sub[y_col].to_numpy(dtype=float)
        ax.scatter(x, y, s=28, color="#2E2E2E", alpha=0.68, edgecolor="white", linewidth=0.3)
        finite = np.isfinite(x) & np.isfinite(y)
        if finite.sum() >= 3:
            slope, intercept, r, p, _ = stats.linregress(x[finite], y[finite])
            xx = np.linspace(np.nanmin(x[finite]), np.nanmax(x[finite]), 100)
            ax.plot(xx, intercept + slope * xx, color="#C43616", lw=1.8)
            ax.text(0.04, 0.96, f"r={r:.2f}, p={p:.2g}\nn={finite.sum()}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
            rows.append({"comparison": parameter, "n_regions": int(finite.sum()), "pearson_r": r, "pearson_p": p})
        ax.axhline(0, color="0.78", lw=0.8)
        ax.axvline(0, color="0.78", lw=0.8)
        ax.set_title(parameter)
        ax.set_xlabel("Syn APP - MAPT")
        ax.set_ylabel("Tau APP - MAPT")
        style_axis(ax)
    fig.suptitle("Are amyloid-induced shifts shared by synuclein and tau?", y=1.04)
    save_figure(fig, figure_dir / "syn_tau_shift_comparison")
    return pd.DataFrame(rows)


def plot_trajectory_examples(tables: dict[str, pd.DataFrame], data_dir: Path, run_root: Path, figure_dir: Path, n_regions: int = 4) -> None:
    for protein in ["syn", "tau"]:
        df = tables[protein][tables[protein]["active_rhat"]].copy()
        df = df[np.isfinite(df["delta_pc1"])].nlargest(n_regions, "delta_pc1")
        app_mean, app_se = load_observation_means(data_dir / f"{protein}_pathology_app.csv")
        mapt_mean, mapt_se = load_observation_means(data_dir / f"{protein}_pathology_mapt.csv")
        app_pred = pd.read_csv(run_root / f"copath_{protein}_app" / "predictions_train.csv").set_index("region")
        mapt_pred = pd.read_csv(run_root / f"copath_{protein}_mapt" / "predictions_train.csv").set_index("region")
        fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), squeeze=False)
        for ax, (_, row) in zip(axes.ravel(), df.iterrows()):
            region = row["region"]
            times = app_mean.index.to_numpy(dtype=float)
            ax.errorbar(times, app_mean[region], yerr=app_se[region], fmt="o", color=PROTEIN_COLOR[protein], label="APP obs")
            ax.errorbar(times, mapt_mean[region], yerr=mapt_se[region], fmt="o", color="#777777", label="MAPT obs")
            pred_times = np.array([float(c) for c in app_pred.columns], dtype=float)
            ax.plot(pred_times, app_pred.loc[region].to_numpy(dtype=float), color=PROTEIN_COLOR[protein], lw=2.2, label="APP fit")
            ax.plot(pred_times, mapt_pred.loc[region].to_numpy(dtype=float), color="#333333", lw=2.2, label="MAPT fit")
            ax.set_title(f"{region}: delta-PC1={row['delta_pc1']:.2f}")
            ax.set_xlabel("MPI")
            ax.set_ylabel("pathology")
            style_axis(ax)
        axes[0, 0].legend(frameon=False, fontsize=8)
        fig.suptitle(f"{PROTEIN_LABEL[protein]} trajectories in high amyloid-sensitive regions", y=1.02)
        save_figure(fig, figure_dir / f"{protein}_high_sensitivity_trajectory_examples")


def load_expression(path: Path) -> tuple[pd.DataFrame, list[str]]:
    expr = pd.read_csv(path)
    expr = expr.rename(columns={expr.columns[0]: "region_base"})
    gene_cols = [c for c in expr.columns if c != "region_base"]
    return expr, gene_cols


def gene_axis_analysis(tables: dict[str, pd.DataFrame], expression_path: Path, out_dir: Path, figure_dir: Path, gmt_path: Path) -> pd.DataFrame:
    expr, gene_cols = load_expression(expression_path)
    summary_rows = []
    for protein, df in tables.items():
        sub = df[df["active_rhat"]][["region", "region_base", "delta_pc1"]].dropna()
        merged = sub.merge(expr, on="region_base", how="inner")
        rows = []
        eta = merged["delta_pc1"].to_numpy(dtype=float)
        for gene in gene_cols:
            y = pd.to_numeric(merged[gene], errors="coerce").to_numpy(dtype=float)
            mask = np.isfinite(y) & np.isfinite(eta)
            if mask.sum() < 8:
                continue
            r, p = stats.pearsonr(y[mask], eta[mask])
            rows.append({"gene": gene, "r": r, "p": p, "n_regions": int(mask.sum())})
        corr = pd.DataFrame(rows).sort_values("r", ascending=False)
        corr["p_fdr"] = bh_fdr(corr["p"].to_numpy())
        protein_dir = out_dir / protein
        protein_dir.mkdir(parents=True, exist_ok=True)
        corr.to_csv(protein_dir / "gene_delta_pc1_correlations.csv", index=False)

        ranked = corr[["gene", "r"]].dropna().drop_duplicates("gene")
        ranked.to_csv(protein_dir / "ranked_genes_delta_pc1.rnk", sep="\t", header=False, index=False)
        pre_res = gp.prerank(
            rnk=ranked,
            gene_sets=str(gmt_path),
            min_size=10,
            max_size=500,
            permutation_num=1000,
            outdir=str(protein_dir / "gseapy"),
            seed=7,
            verbose=False,
        )
        gsea = pre_res.res2d.copy()
        gsea.to_csv(protein_dir / "gsea_delta_pc1_all.csv", index=False)
        gsea["protein"] = protein
        summary_rows.append(gsea)
    all_gsea = pd.concat(summary_rows, ignore_index=True)
    all_gsea.to_csv(out_dir / "gsea_delta_pc1_all_proteins.csv", index=False)
    plot_gsea_dotplot(all_gsea, figure_dir)
    return all_gsea


def plot_gsea_dotplot(gsea: pd.DataFrame, figure_dir: Path, top_n: int = 12) -> None:
    term_col = "Term"
    nes_col = "NES"
    fdr_col = "FDR q-val"
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=False)
    for ax, protein in zip(axes, ["syn", "tau"]):
        sub = gsea[gsea["protein"] == protein].copy()
        sub[fdr_col] = pd.to_numeric(sub[fdr_col], errors="coerce")
        sub[nes_col] = pd.to_numeric(sub[nes_col], errors="coerce")
        sub = sub.sort_values(fdr_col).head(top_n).sort_values(nes_col)
        y = np.arange(len(sub))
        sizes = 30 + 180 * np.clip(-np.log10(sub[fdr_col].fillna(1).to_numpy()), 0, 5) / 5
        sc = ax.scatter(sub[nes_col], y, s=sizes, c=sub[nes_col], cmap="coolwarm", vmin=-3, vmax=3, edgecolor="white", linewidth=0.4)
        ax.axvline(0, color="0.75", lw=0.8)
        ax.set_yticks(y)
        labels = [truncate_label(str(t).replace("KEGG_2019_Mouse__", "").replace("_", " ")) for t in sub[term_col]]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("NES")
        ax.set_title(PROTEIN_LABEL[protein])
        style_axis(ax)
    fig.colorbar(sc, ax=axes, pad=0.02, label="NES")
    fig.suptitle("Pathways aligned with amyloid sensitivity axis", y=1.02)
    save_figure(fig, figure_dir / "gene_pathway_delta_pc1_gsea")


def celltype_axis_analysis(
    tables: dict[str, pd.DataFrame],
    celltype_path: Path,
    out_dir: Path,
    figure_dir: Path,
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    cell = pd.read_csv(celltype_path)
    frac_cols = [c for c in cell.columns if c.startswith("frac_") and c != "frac_row_sum"]
    frac = cell[frac_cols].astype(float).to_numpy()
    frac = np.where(np.isfinite(frac), frac, 0)
    clr = clr_transform(frac)
    clr_df = pd.DataFrame(clr, columns=[c.replace("frac_", "clr_") for c in frac_cols])
    ma_cols = [c for c in MA_CELLTYPES if c in frac_cols]
    ma_idx = [frac_cols.index(c) for c in ma_cols]
    cell = pd.concat([cell[["node", "region_base"]], clr_df], axis=1)
    cell["monoaminergic_score"] = clr[:, ma_idx].mean(axis=1)
    clr_cols = list(clr_df.columns)

    rows = []
    mono_rows = []
    joint_rows = []
    outcome_labels = ["delta_pc1"] + [f"{p}_diff_app_minus_mapt" for p in PARAMETERS]
    for protein, df in tables.items():
        sub = df[df["active_rhat"]][["region", "delta_pc1"] + [f"{p}_diff_app_minus_mapt" for p in PARAMETERS]].dropna(subset=["delta_pc1"])
        merged = sub.merge(cell, left_on="region", right_on="node", how="inner")
        clusters = merged["region_base"].astype(str).to_numpy()
        for outcome_idx, outcome in enumerate(outcome_labels):
            y = merged[outcome].to_numpy(dtype=float)
            for cell_type in clr_cols:
                x = merged[cell_type].to_numpy(dtype=float)
                st = corr_stats(x, y)
                p_perm = cluster_permutation_p(x, y, clusters, float(st["spearman_r"]), n_perm, seed + outcome_idx * 100 + len(rows))
                rows.append({"protein": protein, "outcome": outcome, "cell_type": cell_type.replace("clr_", ""), **st, "p_perm": p_perm})

            ma = merged["monoaminergic_score"].to_numpy(dtype=float)
            st = corr_stats(ma, y)
            p_perm = cluster_permutation_p(ma, y, clusters, float(st["spearman_r"]), n_perm, seed + 10_000 + outcome_idx)
            mono_rows.append(
                {
                    "protein": protein,
                    "outcome": outcome,
                    "score": "monoaminergic",
                    "components": ";".join(ma_cols),
                    **st,
                    "p_perm": p_perm,
                }
            )
        joint = merged[["region", "region_base", "delta_pc1"] + [f"{p}_diff_app_minus_mapt" for p in PARAMETERS]].copy()
        joint["protein"] = protein
        joint["monoaminergic_score"] = merged["monoaminergic_score"]
        joint_rows.append(joint)
    out = pd.DataFrame(rows)
    out["p_fdr"] = out.groupby(["protein", "outcome"])["p_perm"].transform(lambda s: bh_fdr(s.to_numpy()))
    out.to_csv(out_dir / "celltype_delta_axis_correlations.csv", index=False)
    mono = pd.DataFrame(mono_rows)
    mono["p_fdr"] = bh_fdr(mono["p_perm"].to_numpy())
    mono.to_csv(out_dir / "monoaminergic_delta_axis_stats.csv", index=False)
    pd.concat(joint_rows, ignore_index=True).to_csv(out_dir / "celltype_delta_axis_joint_table.csv", index=False)
    plot_celltype_heatmap(out, figure_dir)
    plot_celltype_bars(out, figure_dir)
    plot_monoaminergic_scores(mono, pd.concat(joint_rows, ignore_index=True), figure_dir)
    return out


def plot_celltype_heatmap(corr: pd.DataFrame, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)
    outcomes = ["delta_pc1"] + [f"{p}_diff_app_minus_mapt" for p in PARAMETERS]
    clean = {"delta_pc1": "delta-PC1", "alpha_diff_app_minus_mapt": "delta alpha", "beta_diff_app_minus_mapt": "delta beta", "gamma_diff_app_minus_mapt": "delta gamma"}
    for ax, protein in zip(axes, ["syn", "tau"]):
        sub = corr[corr["protein"] == protein]
        mat = sub.pivot(index="cell_type", columns="outcome", values="spearman_r").reindex(columns=outcomes)
        im = ax.imshow(mat.to_numpy(), cmap="coolwarm", vmin=-0.6, vmax=0.6, aspect="auto")
        ax.set_xticks(np.arange(len(outcomes)))
        ax.set_xticklabels([clean[o] for o in outcomes], rotation=35, ha="right")
        ax.set_yticks(np.arange(len(mat.index)))
        ax.set_yticklabels([c.replace("frac_", "") for c in mat.index])
        ax.set_title(PROTEIN_LABEL[protein])
    fig.colorbar(im, ax=axes, pad=0.02, label="Spearman rho")
    fig.suptitle("Cell-type composition associated with amyloid-sensitive dynamics", y=1.02)
    save_figure(fig, figure_dir / "celltype_delta_axis_heatmap")


def plot_celltype_bars(corr: pd.DataFrame, figure_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 7.0), sharex=True)
    outcomes = ["delta_pc1", "gamma_diff_app_minus_mapt"]
    clean = {"delta_pc1": "delta-PC1", "gamma_diff_app_minus_mapt": "delta gamma"}
    for row_idx, protein in enumerate(["syn", "tau"]):
        for col_idx, outcome in enumerate(outcomes):
            ax = axes[row_idx, col_idx]
            sub = corr[(corr["protein"] == protein) & (corr["outcome"] == outcome)].sort_values("spearman_r")
            labels = sub["cell_type"].to_list()
            colors = ["#C43616" if x < 0 else "#2C6DB2" for x in sub["spearman_r"]]
            ax.barh(labels, sub["spearman_r"], color=colors, alpha=0.88)
            ax.axvline(0, color="black", lw=0.8)
            ax.set_title(f"{PROTEIN_LABEL[protein]} {clean[outcome]}")
            ax.set_xlabel("Spearman rho")
            style_axis(ax)
    fig.tight_layout()
    save_figure(fig, figure_dir / "celltype_delta_axis_bars")


def plot_monoaminergic_scores(mono: pd.DataFrame, joint: pd.DataFrame, figure_dir: Path) -> None:
    outcomes = ["delta_pc1", "alpha_diff_app_minus_mapt", "beta_diff_app_minus_mapt", "gamma_diff_app_minus_mapt"]
    clean = {
        "delta_pc1": "delta-PC1",
        "alpha_diff_app_minus_mapt": "delta alpha",
        "beta_diff_app_minus_mapt": "delta beta",
        "gamma_diff_app_minus_mapt": "delta gamma",
    }
    fig, axes = plt.subplots(2, 4, figsize=(15.2, 6.8), sharey=True)
    for row_idx, protein in enumerate(["syn", "tau"]):
        for col_idx, outcome in enumerate(outcomes):
            ax = axes[row_idx, col_idx]
            sub = joint[joint["protein"] == protein]
            stats_row = mono[(mono["protein"] == protein) & (mono["outcome"] == outcome)].iloc[0]
            x = sub[outcome].to_numpy(dtype=float)
            y = sub["monoaminergic_score"].to_numpy(dtype=float)
            ax.scatter(x, y, s=34, color=PROTEIN_COLOR[protein], alpha=0.72, edgecolor="white", linewidth=0.35)
            finite = np.isfinite(x) & np.isfinite(y)
            if finite.sum() >= 3:
                slope, intercept, *_ = stats.linregress(x[finite], y[finite])
                xx = np.linspace(np.nanmin(x[finite]), np.nanmax(x[finite]), 100)
                ax.plot(xx, intercept + slope * xx, color="black", lw=1.5)
            ax.axvline(0, color="0.78", lw=0.8)
            ax.text(
                0.04,
                0.96,
                rf"$\rho={stats_row['spearman_r']:.2f}$" + f"\nperm p={stats_row['p_perm']:.3g}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9,
            )
            ax.set_title(f"{PROTEIN_LABEL[protein]} {clean[outcome]}")
            ax.set_xlabel(clean[outcome])
            if col_idx == 0:
                ax.set_ylabel("monoaminergic score")
            style_axis(ax)
    fig.suptitle("Monoaminergic score vs APP-induced local RF shifts", y=1.02)
    fig.tight_layout()
    save_figure(fig, figure_dir / "monoaminergic_score_vs_parameter_shifts")


def write_correlations_and_regressions(tables: dict[str, pd.DataFrame], out_dir: Path) -> None:
    corr_rows = []
    reg_rows = []
    for protein, df in tables.items():
        for amyloid in AMYLOIDS:
            amy_col = f"{amyloid}_treatment_mean_prelimval"
            mask = df["active_rhat"].to_numpy()
            for parameter in PARAMETERS + ["delta_pc1"]:
                y_col = f"{parameter}_diff_app_minus_mapt" if parameter in PARAMETERS else parameter
                st = corr_stats(df.loc[mask, amy_col], df.loc[mask, y_col])
                corr_rows.append({"protein": protein, "amyloid": amyloid, "outcome": parameter, **st})
            reg_rows.append(regression_table(df, protein, amyloid))
    corr = pd.DataFrame(corr_rows)
    corr["pearson_p_fdr"] = bh_fdr(corr["pearson_p"].to_numpy())
    corr["spearman_p_fdr"] = bh_fdr(corr["spearman_p"].to_numpy())
    corr.to_csv(out_dir / "amyloid_parameter_shift_and_pc1_correlations.csv", index=False)
    pd.concat(reg_rows, ignore_index=True).to_csv(out_dir / "amyloid_shift_regressions_adjusted_for_baseline.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-dir", default="paper-copath/results/region_rf_condition_comparison")
    parser.add_argument("--data-dir", default="paper-copath/data")
    parser.add_argument("--run-root", default="runs/region_rf")
    parser.add_argument("--expression", default="paper-rf/data/transcriptomics/avg_Pangea_exp.csv")
    parser.add_argument("--cell-types", default="paper-rf/data/cell_types/connectome_celltype.csv")
    parser.add_argument("--gene-sets", default="paper-rf/results/enrichment/striatum/gseapy/gene_sets.gmt")
    parser.add_argument("--out-dir", default="paper-copath/results/amyloid_sensitivity")
    parser.add_argument("--figure-dir", default="paper-copath/figures/amyloid_sensitivity")
    parser.add_argument("--n-celltype-perm", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    figure_dir = Path(args.figure_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    tables = load_comparison_tables(Path(args.comparison_dir))
    pca_rows = []
    for protein in ["syn", "tau"]:
        tables[protein], pca_summary = add_delta_pca(tables[protein], protein)
        tables[protein].to_csv(out_dir / f"{protein}_amyloid_sensitivity_region_table.csv", index=False)
        pca_rows.append(pca_summary)
    pca_summary = pd.concat(pca_rows, ignore_index=True)
    pca_summary.to_csv(out_dir / "delta_parameter_pca_summary.csv", index=False)

    write_correlations_and_regressions(tables, out_dir)
    plot_amyloid_parameter_grid(tables, figure_dir)
    plot_delta_pca(tables, pca_summary, figure_dir)
    plot_top_region_bars(tables, figure_dir)
    plot_syn_tau_comparison(tables, figure_dir).to_csv(out_dir / "syn_tau_shift_correlations.csv", index=False)
    plot_trajectory_examples(tables, Path(args.data_dir), Path(args.run_root), figure_dir)
    gene_axis_analysis(tables, Path(args.expression), out_dir / "transcriptomics", figure_dir, Path(args.gene_sets))
    celltype_axis_analysis(tables, Path(args.cell_types), out_dir, figure_dir, args.n_celltype_perm, args.seed)
    print(out_dir)
    print(figure_dir)


if __name__ == "__main__":
    main()
