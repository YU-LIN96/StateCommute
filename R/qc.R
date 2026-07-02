#' Flag contaminated cells by DecontX score and neuronal signal.
#'
#' Adds a \code{con_level} factor ("low"/"high") and a logical \code{qc_flag}
#' marking cells that are both highly contaminated and express a neuronal marker
#' (e.g. \code{Snap25}), following the ambient-contamination logic used in the
#' microglia workflow. Cells are annotated, not removed; use
#' \code{\link{filter_contaminated}} to subset.
#'
#' @param object            A Seurat object.
#' @param contamination_col Metadata column with DecontX contamination scores.
#' @param threshold         Contamination cutoff (default 0.5).
#' @param neuronal_gene     Neuronal marker gene (default "Snap25").
#' @param neuron_max        Expression above this counts as neuronal (default 1).
#' @return The Seurat object with \code{con_level} and \code{qc_flag} added.
#' @export
flag_contamination <- function(object,
                               contamination_col = "contamination",
                               threshold         = 0.5,
                               neuronal_gene     = "Snap25",
                               neuron_max        = 1) {
  .require("Seurat")
  df     <- Seurat::FetchData(object, vars = c(contamination_col, neuronal_gene))
  contam <- df[[contamination_col]]
  neuro  <- df[[neuronal_gene]]

  object$con_level <- factor(ifelse(contam > threshold, "high", "low"),
                             levels = c("low", "high"))
  object$qc_flag <- (contam > threshold) & (neuro > neuron_max)
  object
}


#' Remove contaminated cells by DecontX score and neuronal signal.
#'
#' Keeps cells with contamination below \code{threshold} AND neuronal marker at
#' or below \code{neuron_max}.
#'
#' @inheritParams flag_contamination
#' @return A subset Seurat object.
#' @export
filter_contaminated <- function(object,
                                contamination_col = "contamination",
                                threshold         = 0.5,
                                neuronal_gene     = "Snap25",
                                neuron_max        = 1) {
  .require("Seurat")
  df   <- Seurat::FetchData(object, vars = c(contamination_col, neuronal_gene))
  keep <- (df[[contamination_col]] < threshold) & (df[[neuronal_gene]] <= neuron_max)
  object[, keep]
}
