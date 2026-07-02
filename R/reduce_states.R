#' Reduce redundant pathways per state, bridging NMF modules to overlap reduction.
#'
#' NMF groups pathways by co-variation across cells; overlap reduction groups them
#' by shared genes. This bridge combines both. Per state it (1) assigns NMF factors
#' to states, (2) unions the top-loading pathways of that state's factors into a
#' candidate set, and (3) collapses gene-overlap redundancy with
#' \code{\link{reduce_pathway_terms}}, using NMF loading as the
#' representative-selection priority.
#'
#' @param object        A Seurat object with an NMF reduction and a state column.
#' @param genesets      Named list of gene sets used for scoring.
#' @param nmf_model     Output of \code{\link{run_nmf}}.
#' @param group.by      State metadata column (default "state_anno2").
#' @param reduction     NMF reduction name (default "nmf").
#' @param n_top         Top-loading pathways per factor to consider (default 30).
#' @param assign        Factor-to-state assignment: "argmax" (each factor to its
#'   most-active state) or "zscore" (threshold; factors may be shared or dropped).
#' @param zscore_thresh z-score cutoff when \code{assign = "zscore"} (default 1).
#' @param min_terms     If a state has fewer candidate pathways than this, skip
#'   reduction and keep them all as representatives (default 3).
#' @param use_loading_priority Pick representatives by NMF loading (default TRUE).
#' @param representative_rule  Used only when \code{use_loading_priority = FALSE}.
#' @param similarity_method,reduction_method,similarity_threshold Passed to
#'   \code{\link{reduce_pathway_terms}}.
#' @param ...           Further args forwarded to \code{\link{reduce_pathway_terms}}.
#' @return A list with \code{representative_table} (one row per state
#'   representative), \code{results} (per-state \code{reduce_pathway_terms}
#'   outputs), and \code{factor_state_map}.
#' @export
reduce_pathways_for_states <- function(
    object, genesets, nmf_model,
    group.by            = "state_anno2",
    reduction           = "nmf",
    n_top               = 30,
    assign              = c("argmax", "zscore"),
    zscore_thresh       = 1,
    min_terms           = 3,
    use_loading_priority = TRUE,
    representative_rule = "central_term",
    similarity_method   = "overlap",
    reduction_method    = "hclust",
    similarity_threshold = 0.7,
    ...
) {
  assign <- match.arg(assign)

  # --- Step 1: factor -> state assignment ---
  fbs    <- nmf_factor_by_state(object, group.by = group.by, reduction = reduction)
  states <- colnames(fbs)

  factor_state_map <- switch(assign,
    argmax = {
      winner <- states[apply(fbs, 1, which.max)]
      split(rownames(fbs), winner)                 # state -> factors
    },
    zscore = {
      z <- t(scale(t(fbs)))                        # row-wise z across states
      z[is.na(z)] <- -Inf                          # guard constant factors
      stats::setNames(
        lapply(states, function(s) rownames(fbs)[z[, s] >= zscore_thresh]),
        states)
    }
  )

  # --- Step 2: candidate pathways per state ---
  tops <- nmf_top_pathways(nmf_model, n = n_top)   # factor, rank, pathway, loading

  results  <- list()
  rep_rows <- list()
  skipped  <- character(0)

  for (s in states) {
    facs <- factor_state_map[[s]]
    if (length(facs) == 0) { skipped <- c(skipped, s); next }

    cand <- tops[tops$factor %in% facs, c("pathway", "loading")]
    cand <- cand[order(-cand$loading), ]
    cand <- cand[!duplicated(cand$pathway), ]                       # keep max loading
    cand <- cand[cand$pathway %in% names(genesets), , drop = FALSE] # must have gene set
    if (nrow(cand) == 0) { skipped <- c(skipped, s); next }

    gs_sub <- genesets[cand$pathway]
    prio   <- stats::setNames(cand$loading, cand$pathway)

    # --- Step 3a: too few to cluster -> keep all as representatives ---
    if (nrow(cand) < min_terms) {
      rep_rows[[s]] <- data.frame(
        state               = s,
        representative_term = cand$pathway,
        n_terms_in_group    = 1L,
        geneset_size        = lengths(gs_sub),
        selection_score     = NA_real_,
        loading             = cand$loading,
        stringsAsFactors    = FALSE
      )
      next
    }

    # --- Step 3b: gene-overlap redundancy reduction ---
    rr <- if (use_loading_priority) "user_priority" else representative_rule
    res <- reduce_pathway_terms(
      genesets             = gs_sub,
      similarity_method    = similarity_method,
      reduction_method     = reduction_method,
      similarity_threshold = similarity_threshold,
      representative_rule  = rr,
      term_priority        = if (use_loading_priority) prio else NULL,
      ...
    )
    results[[s]] <- res

    rt         <- res$representative_table
    rt$state   <- s
    rt$loading <- prio[rt$representative_term]
    rep_rows[[s]] <- rt[, c("state", "representative_term", "n_terms_in_group",
                            "geneset_size", "selection_score", "loading")]
  }

  if (length(skipped)) {
    message(sprintf("No factors/pathways assigned to state(s): %s",
                    paste(skipped, collapse = ", ")))
  }

  list(
    representative_table = dplyr::bind_rows(rep_rows),
    results              = results,
    factor_state_map     = factor_state_map
  )
}
