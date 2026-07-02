#' Run NMF on a pathway x cell matrix (e.g. AUCell scores).
#'
#' Factorizes \eqn{A \approx W H}, where \code{w} = pathways x k (loadings) and
#' \code{h} = k x cells. RcppML does not preserve dimnames, so they are captured
#' from \code{mat} and returned for later re-attachment.
#'
#' @param mat  Non-negative matrix (pathways x cells); a \code{dgCMatrix} is fine.
#' @param k    Number of factors / functional modules (default 20).
#' @param seed RNG seed for reproducibility (default 1234).
#' @param ...  Passed through to \code{RcppML::nmf} (e.g. tol, maxit).
#' @return A list with \code{w}, \code{h}, \code{k}, \code{features}, \code{cells}
#'   and the raw \code{model}.
#' @export
run_nmf <- function(mat, k = 20, seed = 1234, ...) {
  .require("RcppML")
  if (min(mat, na.rm = TRUE) < 0) {
    stop("NMF requires a non-negative matrix (use the AUCell assay).", call. = FALSE)
  }

  set.seed(seed)
  m <- RcppML::nmf(mat, k = k, ...)

  list(
    w        = .nmf_extract(m, "w"),   # pathways x k
    h        = .nmf_extract(m, "h"),   # k x cells
    k        = k,
    features = rownames(mat),
    cells    = colnames(mat),
    model    = m
  )
}


#' Embed an NMF result into a Seurat object as a DimReduc.
#'
#' @param object         A Seurat object.
#' @param nmf_model      Output of \code{\link{run_nmf}}.
#' @param reduction_name Name of the reduction slot (default "nmf").
#' @param key            DimReduc key, must end with "_" (default "NMF_").
#' @param assay          Assay to associate (default "AUCell").
#' @return The Seurat object with the reduction added.
#' @export
embed_nmf <- function(object, nmf_model,
                      reduction_name = "nmf",
                      key            = "NMF_",
                      assay          = "AUCell") {
  .require("Seurat")
  k     <- nmf_model$k
  cells <- nmf_model$cells %||% colnames(object)

  cell_embeddings <- t(nmf_model$h)                      # cells x k
  colnames(cell_embeddings) <- paste0("NMF_", seq_len(k))
  rownames(cell_embeddings) <- cells

  feature_loadings <- nmf_model$w                         # pathways x k
  colnames(feature_loadings) <- paste0("NMF_", seq_len(k))
  rownames(feature_loadings) <- nmf_model$features

  object[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = cell_embeddings,
    loadings   = feature_loadings,
    assay      = assay,
    key        = key
  )
  object
}


#' Top-loading pathways per NMF factor (interpret each module).
#'
#' @param nmf_model Output of \code{\link{run_nmf}}.
#' @param n         Number of top pathways per factor (default 30).
#' @return A long data.frame with columns \code{factor}, \code{rank},
#'   \code{pathway}, \code{loading}.
#' @export
nmf_top_pathways <- function(nmf_model, n = 30) {
  w <- nmf_model$w
  rownames(w) <- nmf_model$features
  k <- nmf_model$k

  do.call(rbind, lapply(seq_len(k), function(j) {
    ord <- order(w[, j], decreasing = TRUE)[seq_len(min(n, nrow(w)))]
    data.frame(
      factor  = paste0("NMF_", j),
      rank    = seq_along(ord),
      pathway = rownames(w)[ord],
      loading = w[ord, j],
      stringsAsFactors = FALSE
    )
  }))
}


#' Mean NMF embedding per group (which module is high in which state).
#'
#' @param object    A Seurat object with an NMF reduction.
#' @param group.by  Metadata column (e.g. "state_anno2").
#' @param reduction Reduction name (default "nmf").
#' @return A matrix of factor x group mean embeddings.
#' @export
nmf_factor_by_state <- function(object, group.by = "state_anno2",
                                reduction = "nmf") {
  .require("Seurat")
  emb <- Seurat::Embeddings(object, reduction = reduction)   # cells x k
  grp <- object[[group.by, drop = TRUE]]
  t(apply(emb, 2, function(col) tapply(col, grp, mean, na.rm = TRUE)))
}


#' Build an SNN graph on the NMF reduction, cluster, and optionally run UMAP.
#'
#' @param object     A Seurat object with an NMF reduction.
#' @param reduction  Reduction to use (default "nmf").
#' @param dims       Dims to use (default \code{1:k} inferred from the reduction).
#' @param resolution Clustering resolution (default 0.3).
#' @param graph.name SNN graph name (default "NMF_snn").
#' @param run_umap   Also compute a UMAP on the reduction (default TRUE).
#' @return The Seurat object with clusters (and \code{nmf_umap}).
#' @export
cluster_states <- function(object, reduction = "nmf",
                           dims = NULL, resolution = 0.3,
                           graph.name = "NMF_snn", run_umap = TRUE) {
  .require("Seurat")
  if (is.null(dims)) {
    dims <- seq_len(ncol(Seurat::Embeddings(object, reduction)))
  }

  object <- Seurat::FindNeighbors(object, reduction = reduction,
                                  dims = dims, graph.name = graph.name)
  object <- Seurat::FindClusters(object, resolution = resolution,
                                 graph.name = graph.name, verbose = FALSE)
  if (run_umap) {
    object <- Seurat::RunUMAP(object, reduction = reduction, dims = dims,
                              reduction.name = "nmf_umap",
                              reduction.key  = "NMFUMAP_")
  }
  object
}
