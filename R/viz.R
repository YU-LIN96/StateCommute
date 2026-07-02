#' Plot the number of retained terms across a range of similarity thresholds.
#'
#' @param genesets          Named list of gene sets.
#' @param sim_mat           Pre-computed similarity matrix (optional speedup).
#' @param thresholds        Numeric vector of thresholds to evaluate.
#' @param similarity_method Used only if \code{sim_mat} is NULL.
#' @param reduction_method  "hclust" or "graph".
#' @param min_geneset_size  Passed to \code{\link{standardize_genesets}}.
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_reduction_sensitivity <- function(
    genesets,
    sim_mat           = NULL,
    thresholds        = seq(0.3, 0.95, by = 0.05),
    similarity_method = "overlap",
    reduction_method  = "hclust",
    min_geneset_size  = 5
) {
  genesets <- standardize_genesets(genesets, min_geneset_size)

  if (is.null(sim_mat)) {
    M       <- build_geneset_sparse_matrix(genesets)
    sim_mat <- compute_similarity(M, method = similarity_method)
  }

  scan_df <- summarize_threshold_scan(sim_mat, thresholds, reduction_method)

  ggplot(scan_df, aes(x = .data$threshold, y = .data$n_selected_terms)) +
    geom_line(linewidth = 1.1, color = "#1565C0") +
    geom_point(size = 3, color = "#1565C0") +
    geom_text(aes(label = .data$n_selected_terms), vjust = -0.8, size = 3.2,
              color = "grey30") +
    labs(title    = "Pathway Term Reduction Sensitivity",
         subtitle = sprintf("Method: %s | Similarity: %s",
                            reduction_method, similarity_method),
         x = "Similarity threshold",
         y = "Number of retained representative terms") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}


#' Histogram of redundancy group sizes.
#'
#' @param result Output of \code{\link{reduce_pathway_terms}}.
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_group_size_distribution <- function(result) {
  df <- result$membership_table %>%
    dplyr::group_by(.data$group_id) %>%
    dplyr::summarise(group_size = dplyr::n(), .groups = "drop")

  ggplot(df, aes(x = .data$group_size)) +
    geom_histogram(binwidth = 1, fill = "#42A5F5", color = "white") +
    labs(title    = "Redundancy Group Size Distribution",
         subtitle = sprintf("%d groups from %d terms  |  threshold = %.2f",
                            nrow(result$representative_table),
                            nrow(result$membership_table),
                            result$params$similarity_threshold),
         x = "Group size (number of terms per group)",
         y = "Count") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}


#' Heatmap of the pairwise similarity matrix.
#'
#' Uses \pkg{ComplexHeatmap} if available, else \pkg{pheatmap}, else base heatmap.
#'
#' @param result        Output of \code{\link{reduce_pathway_terms}}.
#' @param max_terms     Subsample if more terms than this (for readability).
#' @param annotate_reps Highlight representative terms.
#' @return A heatmap object (drawn as a side effect).
#' @export
plot_similarity_heatmap <- function(result,
                                    max_terms     = 100,
                                    annotate_reps = TRUE) {
  sim_mat <- result$similarity_matrix
  if (is.null(sim_mat)) {
    stop("similarity_matrix not found. Re-run with return_similarity_matrix = TRUE.")
  }

  n <- nrow(sim_mat)
  if (n > max_terms) {
    message(sprintf("Subsampling to %d terms for readability.", max_terms))
    idx     <- sort(sample(n, max_terms))
    sim_mat <- sim_mat[idx, idx]
  }

  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {

    col_fun <- circlize::colorRamp2(c(0, 0.5, 1),
                                    c("#FFFFFF", "#90CAF9", "#0D47A1"))

    if (annotate_reps && !is.null(result$membership_table)) {
      is_rep <- result$membership_table$is_representative[
        match(rownames(sim_mat), result$membership_table$term)]
      ra <- ComplexHeatmap::rowAnnotation(
        Representative = ifelse(is_rep, "yes", "no"),
        col = list(Representative = c(yes = "#E53935", no = "grey85")),
        show_annotation_name = FALSE
      )
    } else {
      ra <- NULL
    }

    ComplexHeatmap::Heatmap(
      sim_mat,
      name              = "Similarity",
      col               = col_fun,
      show_row_names    = n <= 60,
      show_column_names = FALSE,
      show_row_dend     = FALSE,
      show_column_dend  = FALSE,
      border            = TRUE,
      clustering_distance_rows    = function(x) stats::as.dist(1 - x),
      clustering_distance_columns = function(x) stats::as.dist(1 - x),
      right_annotation  = ra,
      column_title      = sprintf("Pairwise %s Similarity",
                                  result$params$similarity_method)
    )

  } else if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
      sim_mat,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      show_rownames = n <= 60,
      show_colnames = FALSE,
      color = grDevices::colorRampPalette(c("white", "#90CAF9", "#0D47A1"))(100),
      main  = sprintf("Pairwise %s Similarity",
                      result$params$similarity_method)
    )
  } else {
    stats::heatmap(sim_mat, symm = TRUE,
            col  = grDevices::colorRampPalette(c("white", "#90CAF9", "#0D47A1"))(100),
            main = "Pairwise Similarity")
  }
}


#' Dendrogram with a cutoff line at the similarity threshold (hclust mode only).
#'
#' @param result Output of \code{\link{reduce_pathway_terms}} run with
#'   \code{reduction_method = "hclust"}.
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_dendrogram_cutoff <- function(result) {
  if (result$params$reduction_method != "hclust") {
    stop("Dendrogram plot only available for reduction_method = 'hclust'.")
  }
  .require("ggdendro")

  hc         <- result$clustering_result$hclust_obj
  cut_height <- result$clustering_result$cut_height
  threshold  <- result$params$similarity_threshold

  ddata <- ggdendro::dendro_data(stats::as.dendrogram(hc), type = "rectangle")

  ggplot(ggdendro::segment(ddata)) +
    geom_segment(aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                 color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = cut_height,
               color = "#E53935", linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = 1, y = cut_height + 0.01,
             label = sprintf("threshold = %.2f  (cut height = %.2f)",
                             threshold, cut_height),
             hjust = 0, color = "#E53935", size = 3.5) +
    scale_y_reverse() +
    labs(title    = "Term Redundancy Dendrogram",
         subtitle = sprintf("Distance = 1 - %s similarity  |  linkage = %s",
                            result$params$similarity_method,
                            result$clustering_result$hclust_obj$method),
         x = NULL, y = "Distance (1 - similarity)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x  = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold"))
}


#' Visualise the similarity graph, coloured by redundancy group (graph mode only).
#'
#' Representative terms are drawn larger and labelled.
#'
#' @param result    Output of \code{\link{reduce_pathway_terms}} run with
#'   \code{reduction_method = "graph"}.
#' @param max_nodes Subsample nodes if the graph is very large.
#' @return NULL; the plot is drawn as a side effect.
#' @export
plot_redundancy_graph <- function(result, max_nodes = 80) {
  if (result$params$reduction_method != "graph") {
    stop("Graph plot only available for reduction_method = 'graph'.")
  }

  g   <- result$clustering_result$graph_obj
  mem <- result$membership_table

  if (igraph::vcount(g) > max_nodes) {
    message(sprintf("Subsampling to %d nodes for readability.", max_nodes))
    keep <- sort(sample(igraph::vcount(g), max_nodes))
    g    <- igraph::induced_subgraph(g, keep)
    mem  <- mem[mem$term %in% igraph::V(g)$name, ]
  }

  igraph::V(g)$group  <- mem$group_id[match(igraph::V(g)$name, mem$term)]
  igraph::V(g)$is_rep <- mem$is_representative[match(igraph::V(g)$name, mem$term)]
  igraph::V(g)$size   <- ifelse(igraph::V(g)$is_rep, 10, 5)
  igraph::V(g)$label  <- ifelse(igraph::V(g)$is_rep, igraph::V(g)$name, "")

  n_groups <- max(igraph::V(g)$group, na.rm = TRUE)
  pal      <- scales::hue_pal()(n_groups)
  igraph::V(g)$color <- pal[igraph::V(g)$group]

  set.seed(42)
  plot(g,
       vertex.size        = igraph::V(g)$size,
       vertex.color       = igraph::V(g)$color,
       vertex.label       = igraph::V(g)$label,
       vertex.label.cex   = 0.6,
       vertex.label.color = "black",
       vertex.frame.color = ifelse(igraph::V(g)$is_rep, "black", NA),
       edge.color         = "grey70",
       edge.width         = 0.5,
       layout             = igraph::layout_with_fr(g),
       main               = sprintf(
         "Redundancy Graph  |  threshold = %.2f  |  %d groups",
         result$params$similarity_threshold, n_groups))
  invisible(NULL)
}


#' For a range of thresholds, compute grouping statistics.
#'
#' @param sim_mat          Pre-computed similarity matrix.
#' @param thresholds       Numeric vector of thresholds.
#' @param reduction_method "hclust" or "graph".
#' @return A data.frame of per-threshold group statistics.
#' @export
summarize_threshold_scan <- function(sim_mat,
                                     thresholds       = seq(0.3, 0.95, by = 0.05),
                                     reduction_method = "hclust") {
  do.call(rbind, lapply(thresholds, function(thr) {
    grp <- if (reduction_method == "hclust") {
      cluster_terms_hclust(sim_mat, thr)$group_ids
    } else {
      cluster_terms_graph(sim_mat, thr)$group_ids
    }
    sizes <- table(grp)
    data.frame(
      threshold         = thr,
      n_groups          = length(sizes),
      n_selected_terms  = length(sizes),   # one rep per group
      mean_group_size   = round(mean(sizes), 2),
      median_group_size = stats::median(sizes),
      max_group_size    = max(sizes)
    )
  }))
}
