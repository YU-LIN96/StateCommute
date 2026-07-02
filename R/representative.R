#' Select one representative term per redundancy group.
#'
#' Rules:
#' \describe{
#'   \item{central_term}{Highest mean within-group similarity.}
#'   \item{max_variance}{Highest variance across cells (needs \code{pathway_matrix}).}
#'   \item{largest_geneset}{Largest gene set.}
#'   \item{highest_mean}{Highest mean value across cells (needs \code{pathway_matrix}).}
#'   \item{user_priority}{Highest value in \code{term_priority} (e.g. NMF loading).}
#' }
#' Ties are broken deterministically by larger gene set, then alphabetically.
#'
#' @param group_ids      Named integer vector: term -> group ID.
#' @param sim_mat        Full similarity matrix (terms x terms).
#' @param genesets       Standardised named list of gene sets.
#' @param pathway_matrix Optional numeric matrix (pathways x cells).
#' @param rule           Selection rule (see Details).
#' @param term_priority  Optional named numeric vector for "user_priority".
#' @return A data.frame: \code{term}, \code{group_id}, \code{is_representative},
#'   \code{selection_score}.
#' @export
choose_representative_terms <- function(group_ids,
                                        sim_mat,
                                        genesets,
                                        pathway_matrix  = NULL,
                                        rule = c("central_term",
                                                 "max_variance",
                                                 "largest_geneset",
                                                 "highest_mean",
                                                 "user_priority"),
                                        term_priority   = NULL) {
  rule        <- match.arg(rule)
  terms       <- names(group_ids)
  group_ids_v <- as.integer(group_ids)
  n_groups    <- max(group_ids_v)
  sizes       <- lengths(genesets[terms])

  if (rule %in% c("max_variance", "highest_mean") && is.null(pathway_matrix)) {
    stop(sprintf("rule = '%s' requires pathway_matrix.", rule))
  }

  if (!is.null(pathway_matrix)) {
    shared <- intersect(terms, rownames(pathway_matrix))
    if (length(shared) < length(terms)) {
      warning(sprintf("%d terms not found in pathway_matrix rownames.",
                      length(terms) - length(shared)))
    }
  }

  score_vec <- switch(rule,

    central_term = {
      # Computed per group inside the loop below.
      rep(NA_real_, length(terms))
    },

    max_variance = {
      pm <- pathway_matrix[intersect(terms, rownames(pathway_matrix)), ,
                           drop = FALSE]
      v  <- apply(pm, 1, stats::var, na.rm = TRUE)
      v[terms]
    },

    largest_geneset = {
      sizes
    },

    highest_mean = {
      pm <- pathway_matrix[intersect(terms, rownames(pathway_matrix)), ,
                           drop = FALSE]
      m  <- rowMeans(pm, na.rm = TRUE)
      m[terms]
    },

    user_priority = {
      if (is.null(term_priority)) {
        stop("term_priority must be supplied for rule = 'user_priority'.")
      }
      p <- term_priority[terms]
      p[is.na(p)] <- -Inf
      p
    }
  )
  names(score_vec) <- terms

  results <- lapply(seq_len(n_groups), function(g) {
    grp_terms <- terms[group_ids_v == g]

    if (length(grp_terms) == 1) {
      return(data.frame(
        term              = grp_terms,
        group_id          = g,
        is_representative = TRUE,
        selection_score   = NA_real_,
        stringsAsFactors  = FALSE
      ))
    }

    if (rule == "central_term") {
      sub_sim   <- sim_mat[grp_terms, grp_terms, drop = FALSE]
      mean_sims <- (rowSums(sub_sim) - 1) / (length(grp_terms) - 1)
      score_vec[grp_terms] <<- mean_sims
      best <- grp_terms[which.max(mean_sims)]
    } else {
      grp_scores <- score_vec[grp_terms]
      best <- grp_terms[
        order(-grp_scores, -sizes[grp_terms], grp_terms)[1]
      ]
    }

    data.frame(
      term              = grp_terms,
      group_id          = g,
      is_representative = grp_terms == best,
      selection_score   = score_vec[grp_terms],
      stringsAsFactors  = FALSE
    )
  })

  dplyr::bind_rows(results)
}
