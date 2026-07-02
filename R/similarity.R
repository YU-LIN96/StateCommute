#' Build a sparse binary gene-set x gene membership matrix.
#'
#' Rows are pathways, columns are the union of genes across all sets.
#'
#' @param genesets Standardised named list of gene sets.
#' @return A \code{dgCMatrix} of dimension (n_sets x n_genes).
#' @export
build_geneset_sparse_matrix <- function(genesets) {
  all_genes <- sort(unique(unlist(genesets)))
  n_sets    <- length(genesets)
  n_genes   <- length(all_genes)
  gene_idx  <- stats::setNames(seq_along(all_genes), all_genes)

  set_indices  <- rep(seq_along(genesets), lengths(genesets))
  gene_indices <- gene_idx[unlist(genesets)]

  Matrix::sparseMatrix(
    i    = set_indices,
    j    = gene_indices,
    x    = 1,
    dims = c(n_sets, n_genes),
    dimnames = list(names(genesets), all_genes)
  )
}


# --- Individual metrics (internal) ---------------------------------------------
# Each accepts a sparse binary matrix M (sets x genes) and returns a dense
# symmetric similarity matrix (sets x sets).
#
# Metric choice: the overlap coefficient is preferred when one set is largely
# nested inside another (a small specific pathway inside a broad one). Jaccard
# penalises both sets for non-overlapping members and can under-estimate
# redundancy between nested sets; overlap only penalises by the smaller set size,
# making it more sensitive to containment.

#' Overlap coefficient: |A n B| / min(|A|, |B|). Internal.
#' @keywords internal
#' @noRd
compute_overlap_similarity <- function(M) {
  M            <- methods::as(M, "dgCMatrix")
  intersection <- as.matrix(Matrix::tcrossprod(M))
  sizes        <- Matrix::rowSums(M)
  min_size     <- outer(sizes, sizes, pmin)
  sim          <- intersection / min_size
  diag(sim)    <- 1
  sim
}

#' Jaccard similarity: |A n B| / |A u B|. Internal.
#' @keywords internal
#' @noRd
compute_jaccard_similarity <- function(M) {
  M            <- methods::as(M, "dgCMatrix")
  intersection <- as.matrix(Matrix::tcrossprod(M))
  sizes        <- Matrix::rowSums(M)
  union_size   <- outer(sizes, sizes, "+") - intersection
  sim          <- intersection / union_size
  diag(sim)    <- 1
  sim
}

#' Sorensen-Dice coefficient: 2|A n B| / (|A| + |B|). Internal.
#' @keywords internal
#' @noRd
compute_dice_similarity <- function(M) {
  M            <- methods::as(M, "dgCMatrix")
  intersection <- as.matrix(Matrix::tcrossprod(M))
  sizes        <- Matrix::rowSums(M)
  sum_size     <- outer(sizes, sizes, "+")
  sim          <- 2 * intersection / sum_size
  diag(sim)    <- 1
  sim
}

#' Ochiai (binary cosine): |A n B| / sqrt(|A| * |B|). Internal.
#' @keywords internal
#' @noRd
compute_ochiai_similarity <- function(M) {
  M            <- methods::as(M, "dgCMatrix")
  intersection <- as.matrix(Matrix::tcrossprod(M))
  sizes        <- Matrix::rowSums(M)
  geom_mean    <- sqrt(outer(sizes, sizes, "*"))
  sim          <- intersection / geom_mean
  diag(sim)    <- 1
  sim
}

#' Dispatch pairwise gene-set similarity by method name.
#'
#' @param M      Sparse binary membership matrix (sets x genes).
#' @param method One of "overlap", "jaccard", "dice", "ochiai".
#' @return A dense symmetric similarity matrix.
#' @export
compute_similarity <- function(M,
                               method = c("overlap", "jaccard",
                                          "dice", "ochiai")) {
  method <- match.arg(method)
  switch(method,
         overlap = compute_overlap_similarity(M),
         jaccard = compute_jaccard_similarity(M),
         dice    = compute_dice_similarity(M),
         ochiai  = compute_ochiai_similarity(M)
  )
}
