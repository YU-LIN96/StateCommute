#' Score pathways per cell with irGSEA.
#'
#' For downstream NMF only the AUCell scores are used (they are non-negative).
#' Other methods may be requested for comparison but are optional.
#'
#' @param object      A Seurat object.
#' @param genesets    Named list of gene sets (already feature-filtered).
#' @param assay,slot  Assay/slot to score on (default "RNA"/"data").
#' @param methods     irGSEA methods; AUCell must be included for NMF.
#' @param seeds       RNG seed.
#' @param ncores      Number of cores for irGSEA.
#' @param min.cells,min.feature  Passed to \code{irGSEA::irGSEA.score}.
#' @param kcdf        Kernel for ssgsea-type methods (default "Gaussian").
#' @param ...         Further arguments forwarded to \code{irGSEA::irGSEA.score}.
#' @return The Seurat object with added score assays (e.g. AUCell).
#' @export
score_pathways_irgsea <- function(
    object, genesets,
    assay = "RNA", slot = "data",
    methods = c("AUCell"),
    seeds = 1234, ncores = 1,
    min.cells = 3, min.feature = 0,
    kcdf = "Gaussian", ...
) {
  .require("irGSEA")
  if (!"AUCell" %in% methods) {
    warning("AUCell not in `methods`; get_pathway_cell_matrix() will fail for NMF.")
  }

  irGSEA::irGSEA.score(
    object = object, assay = assay, slot = slot,
    seeds = seeds, ncores = ncores,
    min.cells = min.cells, min.feature = min.feature,
    custom = TRUE, geneset = genesets, msigdb = FALSE,
    geneid = "symbol", method = methods,
    aucell.MaxRank = NULL, ucell.MaxRank = NULL,
    kcdf = kcdf, ...
  )
}


#' Extract the pathway x cell matrix for NMF (non-negative AUCell scores).
#'
#' @param object A Seurat object carrying an AUCell assay.
#' @param assay  Score assay to extract (default "AUCell").
#' @param layer  Layer/slot to pull (default "data").
#' @return A (pathways x cells) matrix, typically a \code{dgCMatrix}.
#' @export
get_pathway_cell_matrix <- function(object, assay = "AUCell", layer = "data") {
  .require("Seurat")
  Seurat::GetAssayData(object, assay = assay, layer = layer)
}
