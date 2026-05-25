# Manuscript Reproduction Map

This map connects manuscript claims and figures to reproducible inputs and scripts in
this repository.

## Core Modeling Figures

These are generated from `runs/` bundles using the Julia package and plotting scripts.

| Manuscript item | Content | Reproduction source |
| --- | --- | --- |
| Fig. 1 / overview | Experimental design and conceptual model schematic | Mostly manually composed figure; model components documented in `src/models.jl` |
| Fig. 2 / transport | Transport mechanism model comparison | Future paper wrapper around transport configs and WAIC summaries |
| Fig. 3 / model comparison | DIFF, DIFF-R, DIFF-RF predictive performance | `scripts/plot_run.jl`, selected merged run bundles |
| Fig. 4 / null models | Connectivity and seed nulls | Future paper wrapper for null configs |
| Fig. 5 / out-of-sample | Leave-final-timepoint-out prediction | Future paper wrapper for holdout configs |
| Fig. S diagnostics | Rhat and chain diagnostics | `scripts/plot_run.jl` diagnostics outputs |
| Fig. S posteriors | Prior/posterior and beta/gamma relationships | `scripts/plot_run.jl` diagnostics plus posterior summaries |

## Biological Interpretation Figures

These are now routed through `paper/analyses/`.

| Manuscript item | Content | Reproduction source |
| --- | --- | --- |
| Main PCA figure | Gene-level beta/gamma coefficient PCA, eta axis, gene ranking | `paper/analyses/transcriptomics/gene_parameter_pca.py` |
| GSEA panels | KEGG pre-ranked enrichment along eta | `paper/analyses/transcriptomics/run_gsea.py` |
| Cell-type panels | CLR cell-type associations and monoaminergic score | `paper/analyses/cell_types/cell_type_axis_associations.py` |
| Hippocampus PCA figure | Replication of eta in hippocampal seeding | Same transcriptomics/cell-type scripts with hippocampus parameter tables |
| Cross-dataset gene comparison | PC direction and gene eta correlation preservation | `paper/analyses/transcriptomics/compare_axes.py` |
| Fig. S PCA filtering | Robustness to filtering choices | `gene_parameter_pca.py --no-update-filter` and `--beta-min` variants |
| Fig. S hippocampal enrichment | Additional hippocampal category/NES panels | `run_gsea.py` plus downstream plotting refinement |

## Current Selected Run Bundles

The current working paper choices are:

| Dataset | Run bundle | Notes |
| --- | --- | --- |
| Striatum DIFF-RF | `runs/striatum_DIFF-RF_RETRO_C1_C3_C4` | Selected chains 1, 3, 4; chain 2 excluded based on lower/different mode |
| Hippocampus DIFF-RF | `runs/hippocampus_DIFF-RF_RETRO_striatum-global-priors_partial_C1_C4` | Partial posterior-prior merge until remaining chains are available |
| Hippocampus normal DIFF-RF | `runs/hippocampus_DIFF-RF_RETRO_C1_C2_C3` | Useful comparison run, not the manuscript posterior-prior configuration |

## Porting Status

Done:

- curated transcriptomic input copied into `paper/data/transcriptomics/`
- curated cell-type input copied into `paper/data/cell_types/`
- posterior beta/gamma table export from run bundles
- gene coefficient PCA and eta ranking
- pre-ranked GSEA wrapper
- cell-type and monoaminergic association script
- striatum/hippocampus axis comparison script
- manuscript-style biological plotting script
- end-to-end paper wrapper that regenerates tables and biological panels

Still to refine:

- WAIC/model-comparison paper wrappers
- null-model and holdout paper wrappers
- full 1000-permutation GSEA outputs for final paper archives
