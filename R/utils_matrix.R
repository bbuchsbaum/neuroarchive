#' Column-wise cumulative sums
#'
#' Computes column-wise cumulative sums of a numeric matrix. If the
#' `matrixStats` package is installed, `matrixStats::colCumsums` is used
#' for efficiency. Otherwise a simple column loop is performed.
#'
#' @param x Numeric matrix.
#' @return Matrix of the same dimensions as `x` containing cumulative sums
#'   down each column.
#' @keywords internal
.col_cumsums <- function(x) {
  stopifnot(is.matrix(x))
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colCumsums(x)
  } else {
    res <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for (j in seq_len(ncol(x))) {
      res[, j] <- cumsum(x[, j])
    }
    res
  }
}
