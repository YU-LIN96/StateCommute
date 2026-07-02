#' Build mouse-native MSigDB gene sets for pathway scoring.
#'
#' Uses the current \pkg{msigdbr} API (\code{collection} / \code{subcollection} /
#' \code{db_species}). Reactome / WikiPathways are pulled from the M2 collection;
#' GO:BP is size-filtered because it is otherwise dominated by very broad terms.
#'
#' @param species     Species passed to msigdbr (default "mouse").
#' @param db_species  Underlying database species (default "MM" for mouse-native).
#' @param specs       List of collection specs. Each element is a list with
#'                    \code{$collection} and optional \code{$subcollection}.
#' @param go_bp_size  Length-2 numeric giving the min/max set size filter applied
#'                    only to the M5:GO:BP collection (default \code{c(20, 200)}).
#' @return Named list of character vectors (gene-deduped, name-deduped).
#' @export
prepare_genesets_msigdb <- function(
    species    = "mouse",
    db_species = "MM",
    specs      = list(
      list(collection = "MH"),
      list(collection = "M2", subcollection = "CP:REACTOME"),
      list(collection = "M2", subcollection = "CP:WIKIPATHWAYS"),
      list(collection = "M7"),
      list(collection = "M5", subcollection = "GO:BP")
    ),
    go_bp_size = c(20, 200)
) {
  .require("msigdbr")

  gs_list <- list()
  for (s in specs) {
    args <- list(species = species, db_species = db_species,
                 collection = s$collection)
    if (!is.null(s$subcollection)) args$subcollection <- s$subcollection
    df <- do.call(msigdbr::msigdbr, args)

    # Size-filter GO:BP only.
    if (identical(s$collection, "M5") &&
        identical(s$subcollection, "GO:BP")) {
      df <- df %>%
        dplyr::group_by(.data$gs_name) %>%
        dplyr::filter(dplyr::n() >= go_bp_size[1] &
                        dplyr::n() <= go_bp_size[2]) %>%
        dplyr::ungroup()
    }
    gs_list <- c(gs_list, split(df$gene_symbol, df$gs_name))
  }

  # Flatten + dedupe genes, then dedupe set names.
  gs_list <- lapply(gs_list, function(x) as.character(unique(unlist(x))))
  names(gs_list) <- make.unique(names(gs_list))
  gs_list
}


#' Standardise a named list of gene sets (object-agnostic).
#'
#' Removes duplicated genes within each set, drops empty/NULL sets, filters sets
#' smaller than \code{min_geneset_size}, and errors on unnamed sets.
#'
#' @param genesets         Named list of character vectors.
#' @param min_geneset_size Integer; minimum gene set size to retain.
#' @return Cleaned named list.
#' @export
standardize_genesets <- function(genesets, min_geneset_size = 5) {
  stopifnot(is.list(genesets))
  if (is.null(names(genesets)) || any(names(genesets) == "")) {
    stop("All elements of genesets must be named.")
  }
  genesets <- lapply(genesets, function(g) unique(g[!is.na(g)]))
  sizes    <- lengths(genesets)
  too_small <- sizes < min_geneset_size

  if (any(too_small)) {
    message(sprintf("Removing %d gene sets smaller than %d genes.",
                    sum(too_small), min_geneset_size))
  }

  genesets <- genesets[!too_small]
  if (length(genesets) == 0) stop("No gene sets remain after filtering.")
  genesets
}


#' Restrict gene sets to features present in an expression matrix (object-aware).
#'
#' The object-aware complement to \code{\link{standardize_genesets}}: intersects
#' each set with a feature universe (e.g. \code{rownames(obj[["RNA"]])}), then
#' re-applies a minimum size filter. This mirrors the cleaning step used before
#' irGSEA scoring so that scores are not diluted by genes absent from the data.
#'
#' @param genesets         Named list of gene sets.
#' @param features         Character vector of the feature universe.
#' @param min_geneset_size Drop sets smaller than this AFTER intersection.
#' @return Filtered named list.
#' @export
filter_genesets_to_features <- function(genesets, features,
                                        min_geneset_size = 5) {
  gs <- lapply(genesets, intersect, y = features)
  gs <- gs[lengths(gs) >= min_geneset_size]
  if (length(gs) == 0) stop("No gene sets remain after feature filtering.")
  message(sprintf("Feature-filtered: %d -> %d gene sets.",
                  length(genesets), length(gs)))
  gs
}
