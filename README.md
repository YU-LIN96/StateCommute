# StateCommute

Pathway-level characterization of disease-associated cell states.

The workflow scores gene sets per cell with **irGSEA (AUCell)**, factorizes the
resulting pathway × cell matrix with **NMF** to recover functional modules,
clusters cells into **states**, and collapses **redundant pathways per state** by
gene-set overlap similarity. A bridge function links the NMF modules to the
overlap-based reduction, using NMF loadings to pick a representative pathway for
each state.

## Installation

Core dependencies are on CRAN. Several stage-specific packages live on
Bioconductor / GitHub and are only needed for the stages that use them:

```r
# core (CRAN)
install.packages(c("dplyr", "ggplot2", "igraph", "magrittr", "Matrix", "scales"))

# scoring / NMF / Seurat side (as needed)
install.packages(c("Seurat", "RcppML", "msigdbr"))
# irGSEA and ComplexHeatmap are on Bioconductor / GitHub:
# BiocManager::install(c("AUCell", "UCell", "ComplexHeatmap"))
# remotes::install_github("chuiqin/irGSEA")

# then install this package
# remotes::install_local("StateCommute")
```

## Pipeline

| Stage | Function(s) |
|-------|-------------|
| Gene sets | `prepare_genesets_msigdb()`, `filter_genesets_to_features()`, `standardize_genesets()` |
| Scoring | `score_pathways_irgsea()`, `get_pathway_cell_matrix()` |
| NMF | `run_nmf()`, `embed_nmf()`, `nmf_top_pathways()`, `nmf_factor_by_state()` |
| Clustering | `cluster_states()` |
| QC (optional) | `flag_contamination()`, `filter_contaminated()` |
| Similarity | `build_geneset_sparse_matrix()`, `compute_similarity()` |
| Redundancy | `cluster_terms_hclust()`, `cluster_terms_graph()`, `choose_representative_terms()` |
| Reduction | `reduce_pathway_terms()`, `build_reduced_terms_simple()` |
| **Bridge** | `reduce_pathways_for_states()` |
| Visualization | `plot_reduction_sensitivity()`, `plot_group_size_distribution()`, `plot_similarity_heatmap()`, `plot_dendrogram_cutoff()`, `plot_redundancy_graph()`, `plot_treemap_simple()` |
| Orchestration | `run_statecommute()` |

## Quick start

```r
library(StateCommute)

# 1. Gene sets -> scoring -> NMF -> clustering
gs  <- prepare_genesets_msigdb()
gs  <- filter_genesets_to_features(gs, rownames(obj[["RNA"]]))
gs  <- standardize_genesets(gs)

obj <- score_pathways_irgsea(obj, gs, methods = "AUCell")
mat <- get_pathway_cell_matrix(obj)
nm  <- run_nmf(mat, k = 20, seed = 1234)
obj <- embed_nmf(obj, nm)
obj <- cluster_states(obj, resolution = 0.3)

# 2. Annotate clusters into states (manual), storing e.g. obj$state_anno2, then:
bridge <- reduce_pathways_for_states(
  obj, genesets = gs, nmf_model = nm,
  group.by = "state_anno2", n_top = 30,
  similarity_method = "overlap", similarity_threshold = 0.7
)

bridge$representative_table                        # per-state non-redundant pathways
plot_similarity_heatmap(bridge$results[["Antiviral"]])

# Or run stages 1-clustering in one call (reduction after you annotate states):
res <- run_statecommute(obj, k = 20, resolution = 0.3)
```

## Notes

- **AUCell only for NMF.** `get_pathway_cell_matrix()` pulls the AUCell assay
  because NMF needs a non-negative matrix; `run_nmf()` errors on negatives.
- **msigdbr API.** Uses the current `collection` / `subcollection` /
  `db_species` arguments.
- **RcppML version.** `run_nmf()` extracts `w`/`h` whether RcppML returns a list
  (older) or an S4 object (newer).
```
