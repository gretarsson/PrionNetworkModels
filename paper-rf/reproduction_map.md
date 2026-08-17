# Manuscript Reproduction Map

This map connects manuscript claims and figures to reproducible inputs and scripts in
this repository.

## Core Modeling Figures

These are generated from `runs/` bundles using the Julia package and plotting scripts.

| Manuscript item | Content | Reproduction source |
| --- | --- | --- |
| Fig. 1 / overview | Experimental design and conceptual model schematic | Mostly manually composed figure; model components documented in `src/models.jl` |
| Fig. 2 / transport | Transport mechanism model comparison | translated `runs/` bundles, `paper-rf/analyses/figures/export_run_bundle_waic.jl`, and `rebuild_figure2.py` |
| Fig. 3 / model comparison | DIFF, DIFF-R, DIFF-RF predictive performance | `scripts/plot_run.jl`, `rebuild_figure3.py`, selected merged run bundles |
| Fig. 4 / null models | Connectivity and seed nulls | translated null WAIC bundles, `export_run_bundle_waic.jl`, and `rebuild_figure4.py` |
| Fig. 5 / out-of-sample | Leave-final-timepoint-out prediction | imported T-1 run bundles and `rebuild_figure5.py` |
| Fig. S diagnostics | Rhat and chain diagnostics | `scripts/plot_run.jl` diagnostics outputs |
| Fig. S posteriors | Prior/posterior and beta/gamma relationships | `scripts/plot_run.jl` diagnostics plus posterior summaries |

## Biological Interpretation Figures

These are now routed through `paper-rf/analyses/`.

| Manuscript item | Content | Reproduction source |
| --- | --- | --- |
| Main PCA figure | Gene-level beta/gamma coefficient PCA, eta axis, gene ranking | `paper-rf/analyses/transcriptomics/gene_parameter_pca.py` |
| GSEA panels | KEGG pre-ranked enrichment along eta | `paper-rf/analyses/transcriptomics/run_gsea.py` |
| Cell-type panels | CLR cell-type associations and monoaminergic score | `paper-rf/analyses/cell_types/cell_type_axis_associations.py` |
| Hippocampus PCA figure | Replication of eta in hippocampal seeding | Same transcriptomics/cell-type scripts with hippocampus parameter tables |
| Cross-dataset gene comparison | PC direction and gene eta correlation preservation | `paper-rf/analyses/transcriptomics/compare_axes.py` |
| Fig. S PCA filtering | Full biological analyses repeated across all, beta-positive, and updated filters | `paper-rf/run_paper_analyses.sh` filter-level outputs |
| Fig. S hippocampal enrichment | Additional hippocampal category/NES panels | `run_gsea.py` plus downstream plotting refinement |
| Fig. S beta/gamma comparison | Filter-specific gene-coefficient comparisons and one shared regional-parameter striatum/hippocampus comparison | `paper-rf/analyses/plotting/plot_biological_figures.py` |

## Current Selected Run Bundles

The current working paper choices are:

| Dataset | Run bundle | Notes |
| --- | --- | --- |
| Striatum DIFF-RF | `runs/striatum_DIFF-RF_RETRO` | Striatal DIFF-RF posterior bundle |
| Hippocampus DIFF-RF | `runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_C3_C4` | Retained hippocampal posterior mode using striatal global posterior priors |
| Hippocampus normal DIFF-RF | `runs/hippocampus_DIFF-RF_RETRO_C1_C2_C3` | Useful comparison run without posterior-derived priors |

## Porting Status

Done:

- curated transcriptomic input copied into `paper-rf/data/transcriptomics/`
- curated cell-type input copied into `paper-rf/data/cell_types/`
- posterior beta/gamma table export from run bundles
- gene coefficient PCA and eta ranking
- pre-ranked GSEA wrapper
- cell-type and monoaminergic association script
- striatum/hippocampus axis comparison script
- manuscript-style biological plotting script
- end-to-end paper wrapper that regenerates tables and biological panels
- manuscript all-region vulnerability-axis analyses and appendix sensitivity panels
- striatum/hippocampus gene-coefficient and regional-parameter comparison panels

Still to refine:

- full 1000-permutation GSEA outputs for final paper archives
