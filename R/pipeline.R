#' Run the StateCommute workflow end to end.
#'
#' Runs gene set preparation (optional), irGSEA scoring, NMF, state clustering,
#' and — if a state annotation column is supplied — per-state pathway reduction.
#' State annotation itself is a manual step between clustering and reduction, so
#' when \code{state_col} is NULL the reduction stage is skipped and you can call
#' \code{\link{reduce_pathways_for_states}} afterwards.
#'
#' @param object               A Seurat object (RNA assay).
#' @param genesets             Named list of gene sets. If NULL, built with
#'   \code{\link{prepare_genesets_msigdb}}.
#' @param species,db_species   Passed to \code{\link{prepare_genesets_msigdb}}.
#' @param score_method         irGSEA method(s); must include "AUCell".
#' @param k                    NMF factors (default 20).
#' @param seed                 RNG seed.
#' @param resolution           Clustering resolution (default 0.3).
#' @param state_col            State annotation column; if NULL, reduction is skipped.
#' @param do_reduction         Attempt per-state reduction (default TRUE).
#' @param n_top                Top-loading pathways per factor for the bridge.
#' @param similarity_method,similarity_threshold,reduction_method Passed to the
#'   per-state reduction.
#' @param ncores               Cores for irGSEA scoring.
#' @param verbose              Print progress (default TRUE).
#' @return A list with \code{object}, \code{genesets}, \code{nmf_model}, and
#'   \code{reduction} (NULL if the reduction stage was skipped).
#' @export
run_statecommute <- function(object,
                             genesets   = NULL,
                             species    = "mouse",
                             db_species = "MM",
                             score_method = "AUCell",
                             k          = 20,
                             seed       = 1234,
                             resolution = 0.3,
                             state_col  = NULL,
                             do_reduction = TRUE,
                             n_top      = 30,
                             similarity_method   = "overlap",
                             similarity_threshold = 0.7,
                             reduction_method    = "hclust",
                             ncores     = 1,
                             verbose    = TRUE) {
  .require("Seurat")
  msg <- function(...) if (verbose) cat(...)

  # 1. Gene sets
  if (is.null(genesets)) {
    msg("Preparing MSigDB gene sets...\n")
    genesets <- prepare_genesets_msigdb(species = species, db_species = db_species)
  }
  genesets <- filter_genesets_to_features(genesets, rownames(object[["RNA"]]))
  genesets <- standardize_genesets(genesets)

  # 2. Score
  msg("Scoring pathways (irGSEA)...\n")
  object <- score_pathways_irgsea(object, genesets, methods = score_method,
                                  seeds = seed, ncores = ncores)

  # 3. NMF
  msg("Running NMF...\n")
  mat       <- get_pathway_cell_matrix(object)
  nmf_model <- run_nmf(mat, k = k, seed = seed)
  object    <- embed_nmf(object, nmf_model)

  # 4. Cluster
  msg("Clustering states on the NMF reduction...\n")
  object <- cluster_states(object, resolution = resolution)

  # 5. Optional per-state reduction
  reduction_out <- NULL
  if (do_reduction && !is.null(state_col)) {
    msg("Reducing redundant pathways per state...\n")
    reduction_out <- reduce_pathways_for_states(
      object, genesets = genesets, nmf_model = nmf_model,
      group.by = state_col, n_top = n_top,
      similarity_method = similarity_method,
      similarity_threshold = similarity_threshold,
      reduction_method = reduction_method
    )
  } else if (do_reduction) {
    msg("Skipping per-state reduction: no `state_col` supplied ",
        "(annotate clusters, then call reduce_pathways_for_states()).\n")
  }

  list(object = object, genesets = genesets,
       nmf_model = nmf_model, reduction = reduction_out)
}
