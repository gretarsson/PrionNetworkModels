#!/usr/bin/env python3
"""Curated TCA/OxPhos machinery correlations with the vulnerability axis."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "paper-rf" / "analyses" / "plotting"))

from plot_manuscript_figure_6_7_panels import setup_style  # noqa: E402


RESULTS_DIR = ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "tca_oxphos_machinery"
FIGURES_DIR = ROOT / "paper-rf" / "figures" / "energy_metabolism"
EXPRESSION = ROOT / "paper-rf" / "data" / "transcriptomics" / "avg_Pangea_exp.csv"
N_PERMUTATIONS = 10000

DATASETS = {
    "striatum": {
        "label": "Striatum",
        "plot_label": "Striatal injection",
        "corr": ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "striatum" / "transcriptomics" / "gene_eta_correlations.csv",
        "axis": ROOT / "paper-rf" / "results" / "pooled_z" / "all" / "striatum" / "transcriptomics" / "region_axis.csv",
        "color": "#8b1e3f",
    },
    "hippocampus": {
        "label": "Hippocampus",
        "plot_label": "Hippocampal injection",
        "corr": ROOT
        / "paper-rf"
        / "results"
        / "pooled_z"
        / "all"
        / "hippocampus_C3_C4"
        / "transcriptomics"
        / "gene_eta_correlations.csv",
        "axis": ROOT
        / "paper-rf"
        / "results"
        / "pooled_z"
        / "all"
        / "hippocampus_C3_C4"
        / "transcriptomics"
        / "region_axis.csv",
        "color": "#00897b",
    },
}

MODULE_ORDER = [
    "TCA cycle",
    "Complex I",
    "Complex II",
    "Complex III",
    "Complex IV",
    "Complex V",
]

MODULE_COLORS = {
    "TCA cycle": "#1b9e77",
    "Complex I": "#274c77",
    "Complex II": "#386fa4",
    "Complex III": "#4f8fc0",
    "Complex IV": "#6aaed6",
    "Complex V": "#9ecae1",
}


def rows(module: str, genes: list[str], role: str, tier: str = "machinery", expected_effect: str = "positive_capacity") -> list[dict]:
    return [
        {
            "gene": gene,
            "module": module,
            "tier": tier,
            "expected_effect": expected_effect,
            "biochemical_role": role,
        }
        for gene in genes
    ]


def machinery_catalog() -> pd.DataFrame:
    entries: list[dict] = []
    entries += rows(
        "TCA cycle",
        ["Cs", "Aco2", "Idh3a", "Idh3b", "Idh3g", "Ogdh", "Ogdhl", "Dlst", "Dld", "Suclg1", "Sucla2", "Suclg2", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Fh1", "Mdh2"],
        "Core TCA-cycle enzymes carrying carbon flux and generating reducing equivalents.",
    )
    entries += rows(
        "Complex I",
        [
            "Ndufa1",
            "Ndufa2",
            "Ndufa3",
            "Ndufa4",
            "Ndufa5",
            "Ndufa6",
            "Ndufa7",
            "Ndufa8",
            "Ndufa9",
            "Ndufa10",
            "Ndufa11",
            "Ndufa12",
            "Ndufa13",
            "Ndufab1",
            "Ndufb1",
            "Ndufb2",
            "Ndufb3",
            "Ndufb4",
            "Ndufb5",
            "Ndufb6",
            "Ndufb7",
            "Ndufb8",
            "Ndufb9",
            "Ndufb10",
            "Ndufb11",
            "Ndufc1",
            "Ndufc2",
            "Ndufs1",
            "Ndufs2",
            "Ndufs3",
            "Ndufs4",
            "Ndufs5",
            "Ndufs6",
            "Ndufs7",
            "Ndufs8",
            "Ndufv1",
            "Ndufv2",
            "Ndufv3",
        ],
        "NADH:ubiquinone oxidoreductase subunits; NADH oxidation and proton pumping.",
    )
    entries += rows(
        "Complex II",
        ["Sdha", "Sdhb", "Sdhc", "Sdhd"],
        "Succinate dehydrogenase; oxidizes succinate and feeds electrons into ubiquinone.",
    )
    entries += rows(
        "Complex III",
        ["Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrb", "Uqcrh", "Uqcrq", "Uqcr10", "Uqcr11", "Cyc1"],
        "Cytochrome bc1 complex and cytochrome c; transfers electrons from ubiquinol toward cytochrome c.",
    )
    entries += rows(
        "Complex IV",
        ["Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l", "Cox7b", "Cox7b2", "Cox7c", "Cox8a", "Cox8b", "Cox8c"],
        "Cytochrome c oxidase nuclear subunits; terminal electron transfer to oxygen and proton pumping.",
    )
    entries += rows(
        "Complex V",
        ["Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5j", "Atp5j2", "Atp5k", "Atp5l", "Atp5md", "Atp5mpl", "Atp5o", "Atp5pb"],
        "ATP synthase subunits; converts proton-motive force into ATP.",
    )

    catalog = pd.DataFrame(entries).drop_duplicates(["gene", "module", "tier"])
    catalog["gene_key"] = catalog["gene"].str.upper()
    catalog["module"] = pd.Categorical(catalog["module"], MODULE_ORDER, ordered=True)
    return catalog.sort_values(["module", "gene"]).reset_index(drop=True)


def read_correlations() -> pd.DataFrame:
    frames = []
    for dataset, spec in DATASETS.items():
        corr = pd.read_csv(spec["corr"])
        corr["gene_key"] = corr["gene"].str.upper()
        corr["dataset"] = dataset
        corr["dataset_label"] = spec["label"]
        frames.append(corr)
    return pd.concat(frames, ignore_index=True)


def merge_catalog(catalog: pd.DataFrame, corr: pd.DataFrame) -> pd.DataFrame:
    merged = catalog.merge(corr, on="gene_key", how="left", suffixes=("_catalog", ""))
    merged["gene_symbol_observed"] = merged["gene"].fillna(merged["gene_catalog"])
    merged["present"] = merged["r"].notna()
    return merged.drop(columns=["gene"]).rename(columns={"gene_catalog": "gene"})


def summarize_modules(merged: pd.DataFrame) -> pd.DataFrame:
    rows_out = []
    for (dataset, dataset_label, module, tier), sub in merged.groupby(["dataset", "dataset_label", "module", "tier"], observed=True):
        vals = sub["r"].dropna().to_numpy()
        if vals.size == 0:
            continue
        n_pos = int((vals > 0).sum())
        try:
            sign_p = stats.binomtest(n_pos, vals.size, 0.5, alternative="greater").pvalue
        except Exception:
            sign_p = np.nan
        t_stat, t_p_two_sided = stats.ttest_1samp(vals, popmean=0.0, nan_policy="omit")
        t_p_greater = t_p_two_sided / 2 if np.isfinite(t_stat) and t_stat > 0 else 1 - (t_p_two_sided / 2)
        wilcoxon_p = stats.wilcoxon(vals, alternative="greater").pvalue if vals.size >= 3 and np.any(vals != 0) else np.nan
        rows_out.append(
            {
                "dataset": dataset,
                "dataset_label": dataset_label,
                "module": module,
                "tier": tier,
                "n_present": int(vals.size),
                "n_positive": n_pos,
                "fraction_positive": n_pos / vals.size,
                "mean_r": float(np.mean(vals)),
                "median_r": float(np.median(vals)),
                "min_r": float(np.min(vals)),
                "max_r": float(np.max(vals)),
                "one_sample_t": float(t_stat),
                "ttest_p_greater_than_zero": float(t_p_greater),
                "sign_test_p_greater_than_half_positive": sign_p,
                "wilcoxon_p_greater_than_zero": wilcoxon_p,
            }
        )
    out = pd.DataFrame(rows_out)
    out["module"] = pd.Categorical(out["module"], MODULE_ORDER, ordered=True)
    return out.sort_values(["dataset", "tier", "module"])


def load_expression_matrix() -> pd.DataFrame:
    expr = pd.read_csv(EXPRESSION)
    expr = expr.rename(columns={expr.columns[0]: "region_base"})
    expr["region_base"] = expr["region_base"].astype(str)
    return expr


def pearson_correlations_against_eta(df: pd.DataFrame, genes: list[str], eta: np.ndarray) -> np.ndarray:
    y = df[genes].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    x = np.asarray(eta, dtype=float)
    out = np.full(y.shape[1], np.nan)
    for j in range(y.shape[1]):
        mask = np.isfinite(x) & np.isfinite(y[:, j])
        if mask.sum() < 5:
            continue
        out[j] = stats.pearsonr(y[mask, j], x[mask]).statistic
    return out


def prepare_expression_z(df: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    y = df[genes].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    valid_rows = np.isfinite(y).all(axis=1)
    y = y[valid_rows]
    y = y - y.mean(axis=0, keepdims=True)
    sd = y.std(axis=0, ddof=1, keepdims=True)
    keep_cols = np.isfinite(sd.ravel()) & (sd.ravel() > 0)
    y = y[:, keep_cols] / sd[:, keep_cols]
    return y, valid_rows


def pearson_corrs_from_z(y_z: np.ndarray, eta: np.ndarray) -> np.ndarray:
    x = np.asarray(eta, dtype=float)
    x = x - x.mean()
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.full(y_z.shape[1], np.nan)
    x = x / sd
    return (y_z.T @ x) / (len(x) - 1)


def permutation_module_tests(catalog: pd.DataFrame, observed_long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = load_expression_matrix()
    rng = np.random.default_rng(20260709)
    summary_rows = []
    null_rows = []

    for dataset, spec in DATASETS.items():
        axis = pd.read_csv(spec["axis"])
        merged = axis.merge(expr, on="region_base", how="inner")
        eta = merged["eta"].to_numpy(dtype=float)

        for module in MODULE_ORDER:
            genes = catalog.loc[catalog["module"].astype(str).eq(module), "gene"].tolist()
            genes = [gene for gene in genes if gene in merged.columns]
            if not genes:
                continue

            observed = observed_long[
                (observed_long["dataset"] == dataset)
                & (observed_long["module"].astype(str) == module)
                & (observed_long["gene"].isin(genes))
            ]["r"].dropna()
            observed_mean = float(observed.mean())
            y_z, valid_rows = prepare_expression_z(merged, genes)
            eta_valid = eta[valid_rows]

            null = np.empty(N_PERMUTATIONS, dtype=float)
            for i in range(N_PERMUTATIONS):
                perm_eta = rng.permutation(eta_valid)
                null[i] = float(np.nanmean(pearson_corrs_from_z(y_z, perm_eta)))

            p_greater = (1 + int(np.sum(null >= observed_mean))) / (N_PERMUTATIONS + 1)
            summary_rows.append(
                {
                    "dataset": dataset,
                    "dataset_label": spec["label"],
                    "module": module,
                    "n_genes": len(genes),
                    "observed_mean_r": observed_mean,
                    "permutation_mean": float(np.mean(null)),
                    "permutation_sd": float(np.std(null, ddof=1)),
                    "permutation_p_greater_than_observed": p_greater,
                    "n_permutations": N_PERMUTATIONS,
                }
            )
            null_rows.extend(
                {
                    "dataset": dataset,
                    "module": module,
                    "permutation": i,
                    "mean_r": value,
                }
                for i, value in enumerate(null)
            )

    summary = pd.DataFrame(summary_rows)
    summary["module"] = pd.Categorical(summary["module"], MODULE_ORDER, ordered=True)
    nulls = pd.DataFrame(null_rows)
    nulls["module"] = pd.Categorical(nulls["module"], MODULE_ORDER, ordered=True)
    return summary.sort_values(["dataset", "module"]), nulls.sort_values(["dataset", "module", "permutation"])


def paired_dataset_table(merged: pd.DataFrame) -> pd.DataFrame:
    cols = ["gene", "module", "tier", "expected_effect", "biochemical_role"]
    wide = (
        merged.pivot_table(index=cols, columns="dataset", values=["r", "p_un", "p_fdr"], aggfunc="first")
        .reset_index()
    )
    wide.columns = ["_".join(c).rstrip("_") if isinstance(c, tuple) else c for c in wide.columns]
    if "r_striatum" in wide and "r_hippocampus" in wide:
        wide["same_direction"] = np.sign(wide["r_striatum"]) == np.sign(wide["r_hippocampus"])
        wide["both_positive"] = (wide["r_striatum"] > 0) & (wide["r_hippocampus"] > 0)
        wide["mean_r_across_datasets"] = wide[["r_striatum", "r_hippocampus"]].mean(axis=1)
    wide["module"] = pd.Categorical(wide["module"], MODULE_ORDER, ordered=True)
    return wide.sort_values(["module", "gene"])


def significance_stars(value: float) -> str:
    if not np.isfinite(value):
        return ""
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return ""


def bold_star_suffix(stars: str) -> str:
    if not stars:
        return ""
    return " (" + "\N{HEAVY ASTERISK}" * len(stars) + ")"


def plot_module_distributions(merged: pd.DataFrame, summary: pd.DataFrame, out: Path) -> None:
    plot_df = merged[(merged["tier"] == "machinery") & merged["present"]].copy()
    modules = [m for m in MODULE_ORDER if (plot_df["module"].astype(str) == m).any()]
    rng = np.random.default_rng(7)
    fig, axes = plt.subplots(2, 1, figsize=(5.1, 9.8), sharex=True, sharey=True)
    for ax, (dataset, spec) in zip(axes, DATASETS.items()):
        sub = plot_df[plot_df["dataset"] == dataset]
        test_sub = summary[(summary["dataset"] == dataset) & (summary["tier"] == "machinery")]
        for i, module in enumerate(modules):
            vals = sub.loc[sub["module"].astype(str) == module, "r"].to_numpy()
            if vals.size == 0:
                continue
            x = i + rng.normal(0, 0.055, vals.size)
            ax.scatter(x, vals, s=58, alpha=0.76, color=MODULE_COLORS[module], linewidths=0)
            med = np.median(vals)
            ax.plot([i - 0.26, i + 0.26], [med, med], color="black", lw=3.6, solid_capstyle="round")
            p_rows = test_sub[test_sub["module"].astype(str) == module]
            if not p_rows.empty:
                stars = significance_stars(float(p_rows["ttest_p_greater_than_zero"].iloc[0]))
                if stars:
                    ax.text(
                        i,
                        0.388,
                        stars,
                        ha="center",
                        va="center",
                        fontsize=17,
                        fontweight="bold",
                        clip_on=False,
                    )
        ax.axhline(0, color="0.45", lw=1.35, ls=(0, (4, 4)))
        ax.set_title(spec.get("plot_label", spec["label"]), fontsize=20, pad=12)
        ax.set_xticks(range(len(modules)))
        ax.set_xticklabels(modules, rotation=28, ha="right", fontsize=15)
        ax.tick_params(axis="y", labelsize=15)
        ax.set_xlabel("")
    for ax in axes:
        ax.set_ylabel("Pearson r with vulnerability axis", fontsize=17)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.12)
    fig.savefig(out.with_suffix(".png"), dpi=250, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def plot_dataset_scatter(wide: pd.DataFrame, out: Path) -> None:
    plot_df = wide[wide["tier"].eq("machinery") & wide["r_striatum"].notna() & wide["r_hippocampus"].notna()].copy()
    fig, ax = plt.subplots(figsize=(8.4, 7.25))
    for module in MODULE_ORDER:
        sub = plot_df[plot_df["module"].astype(str) == module]
        if sub.empty:
            continue
        sign_agreement = 100 * (np.sign(sub["r_striatum"]) == np.sign(sub["r_hippocampus"])).mean()
        ax.scatter(
            sub["r_striatum"],
            sub["r_hippocampus"],
            s=76,
            alpha=0.78,
            label=f"{module} ({sign_agreement:.0f}% sign agreement)",
            color=MODULE_COLORS[module],
            linewidths=0,
        )
    lim = max(abs(plot_df[["r_striatum", "r_hippocampus"]].min().min()), abs(plot_df[["r_striatum", "r_hippocampus"]].max().max())) * 1.08
    ax.axhline(0, color="0.5", lw=1.3, ls=(0, (4, 4)))
    ax.axvline(0, color="0.5", lw=1.3, ls=(0, (4, 4)))
    rho, p = stats.spearmanr(plot_df["r_striatum"], plot_df["r_hippocampus"])
    stars = significance_stars(float(p))
    stars_suffix = bold_star_suffix(stars)
    ax.text(
        0.04,
        0.96,
        rf"Spearman $\rho$={rho:.2f}, p={p:.2g}{stars_suffix}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17.5,
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"Striatum gene correlation with $\eta$", fontsize=20.0)
    ax.set_ylabel(r"Hippocampus gene correlation with $\eta$", fontsize=20.0, labelpad=8)
    ax.tick_params(axis="both", labelsize=17.0)
    ax.legend(frameon=False, fontsize=15.0, ncol=1, loc="lower right", markerscale=1.2)
    fig.subplots_adjust(left=0.27, right=0.96, bottom=0.24, top=0.91)
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=250)
    plt.close(fig)


def plot_heatmap(wide: pd.DataFrame, out: Path) -> None:
    plot_df = wide[wide["tier"].eq("machinery") & wide["r_striatum"].notna() & wide["r_hippocampus"].notna()].copy()
    plot_df["module"] = pd.Categorical(plot_df["module"], MODULE_ORDER, ordered=True)
    plot_df = plot_df.sort_values(["module", "mean_r_across_datasets", "gene"], ascending=[True, False, True])
    data = plot_df[["r_striatum", "r_hippocampus"]].to_numpy()
    vmax = np.nanmax(np.abs(data))
    fig_h = max(7.1, 0.112 * len(plot_df))
    fig = plt.figure(figsize=(3.75, fig_h))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.48, 0.10], wspace=0.035)
    ax = fig.add_subplot(gs[0, 0])
    label_ax = fig.add_subplot(gs[0, 1], sharey=ax)
    cax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(data, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Striatum", "Hippocampus"], fontsize=8.7)
    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["gene"], fontsize=6.2)
    ax.tick_params(axis="x", length=0)
    boundaries = []
    current = None
    for i, module in enumerate(plot_df["module"].astype(str)):
        if module != current:
            boundaries.append((i, module))
            current = module
    for i, module in boundaries:
        if i > 0:
            ax.axhline(i - 0.5, color="black", lw=0.55)
            label_ax.axhline(i - 0.5, color="black", lw=0.55)
        label_ax.text(0.03, i, module, va="top", ha="left", fontsize=10.5, color=MODULE_COLORS.get(module, "black"))
    label_ax.set_xlim(0, 1)
    label_ax.set_xticks([])
    label_ax.grid(False)
    label_ax.tick_params(left=False, labelleft=False)
    for spine in label_ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Pearson r", fontsize=12)
    cbar.ax.tick_params(labelsize=10.5)
    cbar.ax.grid(False)
    fig.subplots_adjust(left=0.23, right=0.97, bottom=0.045, top=0.995)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=250, bbox_inches="tight")
    plt.close(fig)


def plot_permutation_nulls(
    permutation_nulls: pd.DataFrame,
    permutation_summary: pd.DataFrame,
    module_summary: pd.DataFrame,
    out: Path,
) -> None:
    fig, axes = plt.subplots(2, len(MODULE_ORDER), figsize=(15.0, 5.8), sharex=False, sharey=False)
    for row_idx, (dataset, spec) in enumerate(DATASETS.items()):
        for col_idx, module in enumerate(MODULE_ORDER):
            ax = axes[row_idx, col_idx]
            null = permutation_nulls[
                (permutation_nulls["dataset"] == dataset)
                & (permutation_nulls["module"].astype(str) == module)
            ]["mean_r"].to_numpy(dtype=float)
            obs_row = permutation_summary[
                (permutation_summary["dataset"] == dataset)
                & (permutation_summary["module"].astype(str) == module)
            ]
            test_row = module_summary[
                (module_summary["dataset"] == dataset)
                & (module_summary["module"].astype(str) == module)
                & (module_summary["tier"] == "machinery")
            ]
            if null.size == 0 or obs_row.empty:
                ax.axis("off")
                continue
            obs = float(obs_row["observed_mean_r"].iloc[0])
            perm_p = float(obs_row["permutation_p_greater_than_observed"].iloc[0])
            t_p = float(test_row["ttest_p_greater_than_zero"].iloc[0]) if not test_row.empty else np.nan
            ax.hist(null, bins=36, color="0.78", edgecolor="white", linewidth=0.35)
            ax.axvline(obs, color=MODULE_COLORS[module], lw=3.0)
            if row_idx == 0:
                ax.set_title(module, fontsize=13.2, color=MODULE_COLORS[module], pad=7)
            if col_idx == 0:
                ax.set_ylabel(spec["label"], fontsize=13.5)
            ax.text(
                0.03,
                0.91,
                f"t-test p={format_p(t_p)}\nperm p={format_p(perm_p)}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10.6,
                linespacing=1.18,
            )
            ax.tick_params(axis="both", labelsize=9.6)
    fig.supxlabel("Mean Pearson r under eta permutation", y=0.02, fontsize=12.5)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=250, bbox_inches="tight")
    plt.close(fig)


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 1e-4:
        return f"{value:.2e}"
    if value < 0.01:
        return f"{value:.4f}".rstrip("0").rstrip(".")
    return f"{value:.3g}"


def write_stats_text(summary: pd.DataFrame, out: Path) -> None:
    lines = [
        "Energy metabolism module analysis",
        "=================================",
        "",
        "Curated modules: TCA cycle and OxPhos complexes I-V.",
        "Gene-wise statistics are Pearson correlations between regional gene expression and vulnerability axis eta.",
        "T-test: one-sided one-sample t-test of gene-wise correlations against zero, alternative mean r > 0.",
        "Permutation test: eta labels permuted across regions; module mean Pearson r recomputed for each permutation.",
        "",
    ]

    for dataset in DATASETS:
        label = DATASETS[dataset]["label"]
        lines += [label, "-" * len(label)]
        sub = summary[(summary["dataset"] == dataset) & (summary["tier"] == "machinery")].copy()
        sub["module"] = pd.Categorical(sub["module"], MODULE_ORDER, ordered=True)
        for _, row in sub.sort_values("module").iterrows():
            lines.append(
                f"{row['module']}: "
                f"n={int(row['n_present'])}, "
                f"positive={int(row['n_positive'])}/{int(row['n_present'])} ({row['fraction_positive']:.3f}), "
                f"mean r={row['mean_r']:.3f}, "
                f"median r={row['median_r']:.3f}, "
                f"one-sided t p={format_p(float(row['ttest_p_greater_than_zero']))}, "
                f"eta-permutation p={format_p(float(row['permutation_p_greater_than_observed']))} "
                f"({int(row['n_permutations'])} permutations)"
            )
        lines.append("")

    lines += [
        "Interpretation note",
        "-------------------",
        "The striatal dataset shows a robust positive shift for TCA/OxPhos machinery under both the simple t-test and eta-permutation null. The hippocampal dataset shows directionally consistent positive shifts by the simple t-test, but these do not pass the more conservative eta-permutation test.",
        "",
    ]
    out.write_text("\n".join(lines))


def main() -> None:
    setup_style()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    catalog = machinery_catalog()
    corr = read_correlations()
    merged = merge_catalog(catalog, corr)
    summary = summarize_modules(merged)
    wide = paired_dataset_table(merged)
    permutation_summary, permutation_nulls = permutation_module_tests(catalog, merged)
    combined_summary = summary.merge(
        permutation_summary[
            [
                "dataset",
                "module",
                "observed_mean_r",
                "permutation_mean",
                "permutation_sd",
                "permutation_p_greater_than_observed",
                "n_permutations",
            ]
        ],
        on=["dataset", "module"],
        how="left",
    )

    catalog.to_csv(RESULTS_DIR / "curated_tca_oxphos_gene_catalog.csv", index=False)
    merged.to_csv(RESULTS_DIR / "curated_tca_oxphos_gene_eta_correlations_long.csv", index=False)
    wide.to_csv(RESULTS_DIR / "curated_tca_oxphos_gene_eta_correlations_by_gene.csv", index=False)
    combined_summary.to_csv(RESULTS_DIR / "curated_tca_oxphos_module_summary.csv", index=False)
    permutation_summary.to_csv(RESULTS_DIR / "curated_tca_oxphos_module_permutation_summary.csv", index=False)
    permutation_nulls.to_csv(RESULTS_DIR / "curated_tca_oxphos_module_permutation_nulls.csv", index=False)

    plot_module_distributions(merged, summary, FIGURES_DIR / "module_correlation_distributions.pdf")
    plot_dataset_scatter(wide, FIGURES_DIR / "striatal_vs_hippocampal_gene_correlations.pdf")
    plot_heatmap(wide, FIGURES_DIR / "gene_correlation_heatmap.pdf")
    plot_permutation_nulls(permutation_nulls, permutation_summary, summary, FIGURES_DIR / "module_permutation_nulls.pdf")
    write_stats_text(combined_summary, FIGURES_DIR / "energy_metabolism_stats.txt")

    print(RESULTS_DIR)
    print(FIGURES_DIR)
    print(combined_summary[combined_summary["module"].astype(str).isin(MODULE_ORDER) & combined_summary["tier"].eq("machinery")].to_string(index=False))


if __name__ == "__main__":
    main()
