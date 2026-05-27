#!/usr/bin/env python3
"""Gene coefficient PCA for independent REGION-RF striatum/hippocampus fits."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
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


def load_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path)
    expr = expr.rename(columns={expr.columns[0]: "region_base"})
    expr["region_base"] = expr["region_base"].astype(str)
    return expr


def load_region_rf_params(
    run_dir: Path,
    parameters: list[str],
    beta_min: float,
    rhat_threshold: float | None,
) -> pd.DataFrame:
    summary = pd.read_csv(run_dir / "region_rf_summary.csv")
    params = summary[["region", *parameters, "observed_peak_mean"]].copy()
    params["region_base"] = params["region"].map(strip_hemi)
    params["hemi"] = params["region"].astype(str).str[0]
    params["active"] = pd.to_numeric(params["observed_peak_mean"], errors="coerce") > 0

    if rhat_threshold is not None:
        posterior = pd.read_csv(run_dir / "region_rf_posterior_summary_long.csv")
        rhats = posterior[posterior["parameter"].isin(parameters)][["region", "parameter", "rhat"]]
        rhats = rhats.pivot(index="region", columns="parameter", values="rhat").reset_index()
        rhats = rhats.rename(columns={p: f"{p}_rhat" for p in parameters})
        params = params.merge(rhats, on="region", how="left")
    else:
        for p in parameters:
            params[f"{p}_rhat"] = np.nan

    mask = np.isfinite(params[parameters].to_numpy(dtype=float)).all(axis=1)
    mask &= params["active"].to_numpy()
    mask &= params["beta"].to_numpy() > beta_min
    if rhat_threshold is not None:
        for p in parameters:
            mask &= params[f"{p}_rhat"].to_numpy() <= rhat_threshold
    params = params.loc[mask].copy()
    for p in parameters:
        params[f"z_{p}"] = zscore(params[p].to_numpy())
    return params


def align_expression(params: pd.DataFrame, expression: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    merged = params.merge(expression, on="region_base", how="inner")
    meta_cols = set(params.columns)
    gene_cols = [c for c in merged.columns if c not in meta_cols]
    return merged, gene_cols


def regress_gene_coefficients(df: pd.DataFrame, gene_cols: list[str], parameters: list[str]) -> pd.DataFrame:
    x = np.column_stack([df[f"z_{p}"].to_numpy() for p in parameters])
    rows = []
    for gene in gene_cols:
        y = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if mask.sum() < 5:
            continue
        yz = zscore(y[mask])
        coef, *_ = np.linalg.lstsq(x[mask], yz, rcond=None)
        rows.append((gene, *coef, int(mask.sum())))
    return pd.DataFrame(rows, columns=["gene", *[f"coef_{p}" for p in parameters], "n_used"])


def orient_pc1(loadings: np.ndarray, parameters: list[str]) -> np.ndarray:
    gamma_idx = parameters.index("gamma")
    return -loadings if loadings[gamma_idx] < 0 else loadings


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


def save_figure(fig: plt.Figure, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def coef_col(parameter: str) -> str:
    return f"coef_{parameter}"


def plot_gene_coefficients(result_dir: Path, figure_dir: Path, label: str, parameters: list[str]) -> None:
    coefs = pd.read_csv(result_dir / "gene_parameter_coefficients.csv")
    corr = pd.read_csv(result_dir / "gene_eta_correlations.csv")[["gene", "r"]]
    pca_summary = pd.read_csv(result_dir / "pca_summary.csv")
    pc1 = pca_summary.iloc[0]
    coef_cols = [coef_col(p) for p in parameters]
    pca_fit = PCA(n_components=min(3, len(parameters))).fit(coefs[coef_cols].to_numpy())
    gamma_idx = parameters.index("gamma")
    if pca_fit.components_[0, gamma_idx] < 0:
        pca_fit.components_[0, :] *= -1
    df = coefs.merge(corr, on="gene", how="left")

    if len(parameters) == 2:
        pairs = [(parameters[0], parameters[1])]
        figsize = (4.2, 3.8)
    else:
        pairs = [("beta", "gamma"), ("alpha", "gamma"), ("alpha", "beta")]
        pairs = [(x, y) for x, y in pairs if x in parameters and y in parameters]
        figsize = (4.4 * len(pairs), 3.8)

    fig, axes = plt.subplots(1, len(pairs), figsize=figsize, squeeze=False)
    scatter = None
    for ax, (xpar, ypar) in zip(axes[0], pairs):
        xcol = coef_col(xpar)
        ycol = coef_col(ypar)
        scatter = ax.scatter(df[xcol], df[ycol], c=df["r"], cmap="coolwarm", s=5, alpha=0.65, rasterized=True)
        scale = max(df[xcol].abs().quantile(0.995), df[ycol].abs().quantile(0.995))
        xi = parameters.index(xpar)
        yi = parameters.index(ypar)
        ax.plot([0, pca_fit.components_[0, xi] * scale], [0, pca_fit.components_[0, yi] * scale], color="black", lw=2.0)
        if pca_fit.n_components_ > 1:
            ax.plot([0, pca_fit.components_[1, xi] * scale], [0, pca_fit.components_[1, yi] * scale], color="0.35", lw=1.5, ls="--")
        ax.axhline(0, color="0.75", lw=0.7)
        ax.axvline(0, color="0.75", lw=0.7)
        ax.set_xlabel(rf"coefficient for $z({xpar})$")
        ax.set_ylabel(rf"coefficient for $z({ypar})$")
        ax.set_title(f"{xpar} vs {ypar}")
        style_axis(ax)
    explained = "\n".join(
        f"PC{i + 1}: {100 * v:.1f}%" for i, v in enumerate(pca_fit.explained_variance_ratio_[:3])
    )
    axes[0, 0].text(
        0.04,
        0.96,
        explained,
        transform=axes[0, 0].transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 2.5},
    )
    fig.suptitle(f"{label}: gene coefficient PCA", y=1.03)
    cbar = fig.colorbar(scatter, ax=list(axes[0]), pad=0.02)
    cbar.set_label(r"corr(gene, $\eta$)")
    save_figure(fig, figure_dir / "pca_gene_coefficients")


def plot_beta_gamma_eta(result_dir: Path, figure_dir: Path, label: str, parameters: list[str]) -> None:
    if "beta" not in parameters or "gamma" not in parameters:
        return
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
        "$\\eta="
        + "+".join([rf"{pca[f'loading_{p}']:.2f}z({p})" for p in parameters])
        + "$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r"$\eta$")
    style_axis(ax)
    save_figure(fig, figure_dir / "beta_gamma_colored_by_eta")


def corr_text(x: pd.Series, y: pd.Series) -> str:
    mask = np.isfinite(x.to_numpy(dtype=float)) & np.isfinite(y.to_numpy(dtype=float))
    if mask.sum() < 3:
        return f"n={mask.sum()}"
    pearson = stats.pearsonr(x[mask], y[mask])
    spearman = stats.spearmanr(x[mask], y[mask])
    return f"r={pearson.statistic:.2f}, p={pearson.pvalue:.2g}\nrho={spearman.statistic:.2f}, p={spearman.pvalue:.2g}\nn={mask.sum()}"


def plot_parameter_pairs(result_dir: Path, figure_dir: Path, label: str, parameters: list[str]) -> None:
    pairs = [("beta", "gamma"), ("alpha", "gamma"), ("alpha", "beta")]
    pairs = [(x, y) for x, y in pairs if x in parameters and y in parameters]
    if not pairs:
        return
    df = pd.read_csv(result_dir / "filtered_parameters.csv")
    fig, axes = plt.subplots(1, len(pairs), figsize=(4.2 * len(pairs), 3.6), squeeze=False)
    for ax, (xpar, ypar) in zip(axes[0], pairs):
        ax.scatter(df[xpar], df[ypar], s=34, color="#0047AB", alpha=0.7, edgecolor="white", linewidth=0.35)
        ax.axhline(0, color="0.75", lw=0.7)
        ax.axvline(0, color="0.75", lw=0.7)
        ax.set_xlabel(xpar)
        ax.set_ylabel(ypar)
        ax.set_title(f"{ypar} vs {xpar}")
        ax.text(
            0.04,
            0.96,
            corr_text(df[xpar], df[ypar]),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 2.5},
        )
        style_axis(ax)
    fig.suptitle(f"{label}: regional parameter pairs", y=1.03)
    save_figure(fig, figure_dir / "regional_parameter_pairs")


def run_dataset(
    name: str,
    run_dir: Path,
    expression: pd.DataFrame,
    out_root: Path,
    figure_root: Path,
    parameters: list[str],
    beta_min: float,
    rhat_threshold: float | None,
) -> dict[str, float | str | int]:
    params = load_region_rf_params(run_dir, parameters, beta_min, rhat_threshold)
    merged, gene_cols = align_expression(params, expression)
    coefs = regress_gene_coefficients(merged, gene_cols, parameters)

    coef_cols = [coef_col(p) for p in parameters]
    n_components = min(3, len(parameters))
    pca = PCA(n_components=n_components)
    pca.fit(coefs[coef_cols].to_numpy())
    pc1 = orient_pc1(pca.components_[0].copy(), parameters)
    if not np.allclose(pc1, pca.components_[0]):
        pca.components_[0, :] *= -1
    zmat = np.column_stack([merged[f"z_{p}"].to_numpy() for p in parameters])
    eta = zmat @ pc1
    merged["eta"] = eta
    coefs["pc1_score"] = coefs[coef_cols].to_numpy() @ pc1

    out_dir = out_root / name
    fig_dir = figure_root / name
    out_dir.mkdir(parents=True, exist_ok=True)
    params.to_csv(out_dir / "filtered_parameters.csv", index=False)
    merged[["region", "region_base", "hemi", *parameters, *[f"z_{p}" for p in parameters], "eta"]].to_csv(
        out_dir / "region_axis.csv",
        index=False,
    )
    coefs.to_csv(out_dir / "gene_parameter_coefficients.csv", index=False)
    gene_correlations(merged, gene_cols, eta).to_csv(out_dir / "gene_eta_correlations.csv", index=False)
    pca_rows = []
    for idx in range(n_components):
        component = pca.components_[idx].copy()
        row = {
            "component": f"PC{idx + 1}",
            "explained_variance_ratio": pca.explained_variance_ratio_[idx],
            "n_regions": len(merged),
            "n_genes": len(coefs),
        }
        for pidx, parameter in enumerate(parameters):
            row[f"loading_{parameter}"] = component[pidx]
        pca_rows.append(row)
    pd.DataFrame(pca_rows).to_csv(out_dir / "pca_summary.csv", index=False)
    plot_gene_coefficients(out_dir, fig_dir, name.capitalize(), parameters)
    plot_beta_gamma_eta(out_dir, fig_dir, name.capitalize(), parameters)
    plot_parameter_pairs(out_dir, fig_dir, name.capitalize(), parameters)
    result: dict[str, float | str | int] = {
        "dataset": name,
        "pc1_explained": pca.explained_variance_ratio_[0],
        "n_regions": len(merged),
        "n_genes": len(coefs),
    }
    if n_components > 1:
        result["pc2_explained"] = pca.explained_variance_ratio_[1]
    if n_components > 2:
        result["pc3_explained"] = pca.explained_variance_ratio_[2]
    for pidx, parameter in enumerate(parameters):
        result[f"loading_{parameter}"] = pc1[pidx]
    return {
        **result,
    }


def compare_directions(summary: pd.DataFrame, parameters: list[str]) -> pd.DataFrame:
    if not {"striatum", "hippocampus"}.issubset(set(summary["dataset"])):
        return pd.DataFrame()
    striatum = summary[summary["dataset"] == "striatum"].iloc[0]
    hippocampus = summary[summary["dataset"] == "hippocampus"].iloc[0]
    v_str = np.array([striatum[f"loading_{p}"] for p in parameters], dtype=float)
    v_hip = np.array([hippocampus[f"loading_{p}"] for p in parameters], dtype=float)
    cosine = float(np.clip(np.dot(v_str, v_hip), -1, 1))
    row = {
        "pc1_cosine_similarity": cosine,
        "pc1_absolute_angle_degrees": float(np.degrees(np.arccos(abs(cosine)))),
        "striatum_pc1_explained": striatum.pc1_explained,
        "hippocampus_pc1_explained": hippocampus.pc1_explained,
        "striatum_n_regions": striatum.n_regions,
        "hippocampus_n_regions": hippocampus.n_regions,
    }
    for parameter in parameters:
        row[f"striatum_loading_{parameter}"] = striatum[f"loading_{parameter}"]
        row[f"hippocampus_loading_{parameter}"] = hippocampus[f"loading_{parameter}"]
    return pd.DataFrame([row])


def plot_direction_comparison(
    summary: pd.DataFrame, comparison: pd.DataFrame, out_path: Path, parameters: list[str]
) -> None:
    if comparison.empty:
        return
    if "beta" not in parameters or "gamma" not in parameters:
        return
    row = comparison.iloc[0]
    colors = {"striatum": "#0047AB", "hippocampus": "#C43616"}
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    ax.axhline(0, color="0.85", lw=1)
    ax.axvline(0, color="0.85", lw=1)
    circle = plt.Circle((0, 0), 1, color="0.9", fill=False, lw=1)
    ax.add_patch(circle)
    for _, s in summary.iterrows():
        ax.arrow(
            0,
            0,
            s["loading_beta"],
            s["loading_gamma"],
            head_width=0.035,
            length_includes_head=True,
            color=colors.get(s["dataset"], "0.25"),
            lw=2.5,
        )
        ax.text(
            1.08 * s["loading_beta"],
            1.08 * s["loading_gamma"],
            f"{s['dataset'].capitalize()}\nPC1 {100*s['pc1_explained']:.1f}%",
            color=colors.get(s["dataset"], "0.25"),
            ha="center",
            va="center",
            fontsize=9,
        )
    ax.set_title(f"REGION-RF PC1 directions\nangle={row['pc1_absolute_angle_degrees']:.1f} deg")
    ax.set_xlabel("PC1 loading on z(beta)")
    ax.set_ylabel("PC1 loading on z(gamma)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.1, 1.05)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    save_figure(fig, out_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", default="paper-rf/data/transcriptomics/avg_Pangea_exp.csv")
    parser.add_argument("--striatum-run", default="runs/region_rf/striatum")
    parser.add_argument("--hippocampus-run", default="runs/region_rf/hippocampus")
    parser.add_argument("--out-dir", default="paper-rf/results/region_rf_gene_pca")
    parser.add_argument("--figure-dir", default="paper-rf/figures/region_rf_gene_pca")
    parser.add_argument("--parameters", default="beta,gamma")
    parser.add_argument("--beta-min", type=float, default=0.0)
    parser.add_argument("--no-rhat-filter", action="store_true")
    parser.add_argument("--rhat-threshold", type=float, default=1.05)
    args = parser.parse_args()

    expression = load_expression(Path(args.expression))
    out_root = Path(args.out_dir)
    figure_root = Path(args.figure_dir)
    parameters = [p.strip() for p in args.parameters.split(",") if p.strip()]
    if "gamma" not in parameters:
        raise ValueError("--parameters must include gamma so PC1 can use the paper fall-positive orientation")
    rhat_threshold = None if args.no_rhat_filter else args.rhat_threshold
    rows = [
        run_dataset(
            "striatum",
            Path(args.striatum_run),
            expression,
            out_root,
            figure_root,
            parameters,
            args.beta_min,
            rhat_threshold,
        ),
        run_dataset(
            "hippocampus",
            Path(args.hippocampus_run),
            expression,
            out_root,
            figure_root,
            parameters,
            args.beta_min,
            rhat_threshold,
        ),
    ]
    summary = pd.DataFrame(rows)
    out_root.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_root / "region_rf_pca_summary.csv", index=False)
    comparison = compare_directions(summary, parameters)
    comparison.to_csv(out_root / "pc1_direction_comparison.csv", index=False)
    plot_direction_comparison(summary, comparison, figure_root / "pc1_direction_comparison", parameters)
    print(out_root)


if __name__ == "__main__":
    main()
