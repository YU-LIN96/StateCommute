#' Reduce pathway terms by gene-set overlap similarity.
#'
#' Standardises gene sets, builds a sparse membership matrix, computes pairwise
#' similarity, groups redundant terms, and picks one representative per group.
#'
#' @param genesets               Named list of gene sets.
#' @param pathway_matrix         Optional pathway x cell numeric matrix.
#' @param similarity_method      "overlap" | "jaccard" | "dice" | "ochiai".
#' @param reduction_method       "hclust" | "graph".
#' @param similarity_threshold   Minimum similarity to be considered redundant.
#' @param hclust_method          Linkage for hclust (default "average").
#' @param community_detection    For graph mode: "components" | "louvain".
#' @param representative_rule    How to pick the representative term per group.
#' @param term_priority          Named numeric vector (for "user_priority").
#' @param min_geneset_size       Drop gene sets smaller than this.
#' @param return_similarity_matrix Include the full similarity matrix in output?
#' @param return_membership      Include the membership data.frame in output?
#' @return A named list: \code{selected_terms}, \code{representative_table},
#'   \code{clustering_result}, \code{params}, and optionally
#'   \code{membership_table} / \code{similarity_matrix}.
#' @export
reduce_pathway_terms <- function(
    genesets,
    pathway_matrix        = NULL,
    similarity_method     = c("overlap", "jaccard", "dice", "ochiai"),
    reduction_method      = c("hclust", "graph"),
    similarity_threshold  = 0.7,
    hclust_method         = "average",
    community_detection   = c("components", "louvain"),
    representative_rule   = c("central_term", "max_variance",
                              "largest_geneset", "highest_mean",
                              "user_priority"),
    term_priority         = NULL,
    min_geneset_size      = 5,
    return_similarity_matrix = TRUE,
    return_membership     = TRUE
) {
  similarity_method   <- match.arg(similarity_method)
  reduction_method    <- match.arg(reduction_method)
  community_detection <- match.arg(community_detection)
  representative_rule <- match.arg(representative_rule)

  # --- Step 1: Standardise gene sets ---
  cat("Standardising gene sets...\n")
  genesets <- standardize_genesets(genesets, min_geneset_size)
  n_terms  <- length(genesets)
  cat(sprintf("  %d gene sets retained.\n", n_terms))

  # --- Step 2: Build sparse membership matrix ---
  cat("Building sparse membership matrix...\n")
  M <- build_geneset_sparse_matrix(genesets)

  # --- Step 3: Compute pairwise similarity ---
  cat(sprintf("Computing %s similarity (%d x %d)...\n",
              similarity_method, n_terms, n_terms))
  sim_mat <- compute_similarity(M, method = similarity_method)

  # --- Step 4: Redundancy grouping ---
  cat(sprintf("Grouping terms via %s (threshold = %.2f)...\n",
              reduction_method, similarity_threshold))

  if (reduction_method == "hclust") {
    clust_result <- cluster_terms_hclust(sim_mat, similarity_threshold,
                                         hclust_method)
  } else {
    clust_result <- cluster_terms_graph(sim_mat, similarity_threshold,
                                        community_detection)
  }
  group_ids <- clust_result$group_ids

  # --- Step 5: Representative selection ---
  cat(sprintf("Selecting representatives by '%s'...\n", representative_rule))
  membership_df <- choose_representative_terms(
    group_ids      = group_ids,
    sim_mat        = sim_mat,
    genesets       = genesets,
    pathway_matrix = pathway_matrix,
    rule           = representative_rule,
    term_priority  = term_priority
  )

  # Attach gene set sizes
  membership_df$geneset_size <- lengths(genesets[membership_df$term])

  # --- Step 6: Build representative table ---
  selected_terms <- membership_df$term[membership_df$is_representative]

  rep_table <- membership_df %>%
    dplyr::filter(.data$is_representative) %>%
    dplyr::left_join(
      membership_df %>%
        dplyr::group_by(.data$group_id) %>%
        dplyr::summarise(n_terms_in_group = dplyr::n(), .groups = "drop"),
      by = "group_id"
    ) %>%
    dplyr::mutate(selection_rule = representative_rule) %>%
    dplyr::rename(representative_term = "term") %>%
    dplyr::select("representative_term", "group_id",
                  "n_terms_in_group", "geneset_size",
                  "selection_rule", "selection_score")

  member_list <- membership_df %>%
    dplyr::left_join(
      membership_df %>%
        dplyr::filter(.data$is_representative) %>%
        dplyr::select("group_id", representative_term = "term"),
      by = "group_id"
    ) %>%
    dplyr::group_by(.data$representative_term) %>%
    dplyr::summarise(member_terms = list(.data$term), .groups = "drop")

  rep_table <- rep_table %>%
    dplyr::left_join(member_list, by = "representative_term")

  cat(sprintf("Done. %d terms -> %d groups -> %d representatives.\n",
              n_terms, max(group_ids), length(selected_terms)))

  out <- list(
    selected_terms       = selected_terms,
    representative_table = rep_table,
    clustering_result    = clust_result,
    params = list(
      similarity_method    = similarity_method,
      reduction_method     = reduction_method,
      similarity_threshold = similarity_threshold,
      representative_rule  = representative_rule,
      min_geneset_size     = min_geneset_size
    )
  )

  if (return_membership)        out$membership_table  <- membership_df
  if (return_similarity_matrix) out$similarity_matrix <- sim_mat

  out
}


#' Build a compact reduced-terms table (parent/child) for treemap-style plots.
#'
#' A lighter-weight alternative to \code{\link{reduce_pathway_terms}} that works
#' from a pre-computed similarity matrix and returns a tidy parent/child table.
#'
#' @param terms              Character vector of terms to reduce.
#' @param sim_mat            Pre-computed similarity matrix.
#' @param genesets           Named list of gene sets.
#' @param threshold          Similarity threshold (default 0.3).
#' @param method             "graph" (default) or "hclust".
#' @param community_detection For graph mode: "components" or "louvain".
#' @param hclust_method      Linkage for hclust.
#' @param rep_rule           "central_term" or "largest_geneset".
#' @param score_vector       Optional named numeric vector of per-term scores.
#' @return A data.frame with \code{term}, \code{parent}, \code{parentTerm},
#'   \code{group}, \code{size}, \code{selection_score}, \code{score}.
#' @export
build_reduced_terms_simple <- function(
    terms,
    sim_mat,
    genesets,
    threshold = 0.3,
    method = c("graph", "hclust"),
    community_detection = c("components", "louvain"),
    hclust_method = "average",
    rep_rule = c("central_term", "largest_geneset"),
    score_vector = NULL
) {
  method              <- match.arg(method)
  community_detection <- match.arg(community_detection)
  rep_rule            <- match.arg(rep_rule)

  # keep shared terms only
  terms <- intersect(terms, rownames(sim_mat))
  terms <- intersect(terms, names(genesets))

  sim_sub      <- sim_mat[terms, terms, drop = FALSE]
  genesets_sub <- genesets[terms]

  group_ids <- if (method == "graph") {
    cluster_terms_graph(
      sim_mat = sim_sub,
      similarity_threshold = threshold,
      community_detection = community_detection
    )$group_ids
  } else {
    cluster_terms_hclust(
      sim_mat = sim_sub,
      similarity_threshold = threshold,
      hclust_method = hclust_method
    )$group_ids
  }

  group_ids <- group_ids[terms]

  rep_df <- choose_representative_terms(
    group_ids = group_ids,
    sim_mat   = sim_sub,
    genesets  = genesets_sub,
    rule      = rep_rule
  )

  rep_map <- rep_df %>%
    dplyr::filter(.data$is_representative) %>%
    dplyr::select("group_id", parent = "term")

  reduced_terms <- rep_df %>%
    dplyr::left_join(rep_map, by = "group_id") %>%
    dplyr::mutate(
      parentTerm = .data$parent,
      group = .data$group_id,
      size  = lengths(genesets_sub[.data$term])
    ) %>%
    dplyr::select("term", "parent", "parentTerm", "group", "size",
                  "selection_score")

  if (is.null(score_vector)) {
    reduced_terms$score <- ifelse(is.na(reduced_terms$selection_score),
                                  1, reduced_terms$selection_score)
  } else {
    score_vector        <- score_vector[terms]
    reduced_terms$score <- score_vector[reduced_terms$term]
  }

  reduced_terms
}


#' Treemap of a reduced-terms table.
#'
#' @param reduced_terms Output of \code{\link{build_reduced_terms_simple}}.
#' @param size          Column to size rectangles by ("size" or "score").
#' @param title         Plot title.
#' @return Invisibly, the treemap object (drawn as a side effect).
#' @export
plot_treemap_simple <- function(reduced_terms, size = "size", title = "") {
  .require("treemap")
  size <- match.arg(size, c("size", "score"))

  treemap::treemap(
    reduced_terms,
    index = c("parentTerm", "term"),
    vSize = size,
    type  = "index",
    title = title,
    palette = scales::hue_pal()(length(unique(reduced_terms$parent))),
    fontcolor.labels = c("#FFFFFFDD", "#00000080"),
    bg.labels  = 0,
    border.col = "#00000080"
  )
}
