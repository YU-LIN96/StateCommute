#' Group pathway terms by hierarchical clustering on the similarity matrix.
#'
#' Internally distance = 1 - similarity, and the tree is cut at height
#' 1 - \code{similarity_threshold}. Terms in the same cluster are considered
#' redundant enough to merge.
#'
#' @param sim_mat              Dense similarity matrix (terms x terms).
#' @param similarity_threshold Minimum similarity to be considered redundant (0-1).
#' @param hclust_method        Linkage passed to \code{stats::hclust} (default "average").
#' @return A list with \code{group_ids} (named integer vector), \code{hclust_obj},
#'   \code{cut_height}, \code{threshold}.
#' @export
cluster_terms_hclust <- function(sim_mat,
                                 similarity_threshold = 0.7,
                                 hclust_method        = "average") {
  dist_mat <- stats::as.dist(1 - sim_mat)
  hc       <- stats::hclust(dist_mat, method = hclust_method)

  cut_height <- 1 - similarity_threshold
  groups     <- stats::cutree(hc, h = cut_height)

  list(group_ids  = groups,
       hclust_obj = hc,
       cut_height = cut_height,
       threshold  = similarity_threshold)
}


#' Group pathway terms by connected components in a similarity graph.
#'
#' An edge connects two terms if their similarity is at least the threshold.
#' Connected components define redundancy groups; Louvain gives finer grouping.
#'
#' @param sim_mat              Dense similarity matrix.
#' @param similarity_threshold Edge threshold.
#' @param community_detection  "components" (default) or "louvain".
#' @return A list with \code{group_ids} (named integer vector), \code{graph_obj}
#'   (igraph), \code{threshold}, \code{method}.
#' @export
cluster_terms_graph <- function(sim_mat,
                                similarity_threshold = 0.7,
                                community_detection  = c("components",
                                                         "louvain")) {
  community_detection <- match.arg(community_detection)

  adj       <- sim_mat
  diag(adj) <- 0
  adj[adj < similarity_threshold] <- 0

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)

  if (community_detection == "components") {
    groups <- igraph::components(g)$membership
  } else {
    groups <- igraph::membership(igraph::cluster_louvain(g))
  }

  names(groups) <- rownames(sim_mat)
  list(group_ids = groups,
       graph_obj = g,
       threshold = similarity_threshold,
       method    = community_detection)
}
