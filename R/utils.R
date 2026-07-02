#' NULL-coalescing operator (internal).
#' @param a,b Values; returns \code{a} unless it is NULL, else \code{b}.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Require a Suggested package or stop with an informative message (internal).
#' @param pkg Package name.
#' @keywords internal
#' @noRd
.require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this function. Please install it.",
                 pkg), call. = FALSE)
  }
  invisible(TRUE)
}

#' Extract a component from an NMF result, whether it is a list or an S4 object
#' (internal). RcppML returns a list in older versions and an S4 object in newer
#' ones; this accessor handles both.
#' @param m    NMF result object.
#' @param name Component name ("w" or "h").
#' @keywords internal
#' @noRd
.nmf_extract <- function(m, name) {
  out <- tryCatch(m[[name]], error = function(e) NULL)
  if (is.null(out))
    out <- tryCatch(methods::slot(m, name), error = function(e) NULL)
  if (is.null(out))
    stop(sprintf("Could not extract '%s' from the NMF result.", name), call. = FALSE)
  out
}

#' Pipe operator
#'
#' Re-exported from \pkg{magrittr}. See \code{\link[magrittr]{\%>\%}}.
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @usage lhs \%>\% rhs
#' @param lhs A value or the result of a previous pipe.
#' @param rhs A function call using the pipe semantics.
#' @return The result of calling \code{rhs(lhs)}.
NULL
