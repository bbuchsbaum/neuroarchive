#' Temporal Reconstruction Functions
#'
#' This file contains all temporal reconstruction functions that reconstruct
#' data from temporal coefficients using various basis types.

#' Reconstruct data from temporal coefficients
#'
#' Mirrors `temporal_project` for the inverse operation. The default
#' simply multiplies the basis by the coefficient matrix.
#' @keywords internal
temporal_reconstruct <- function(kind, basis, coeff, ...) {
  stopifnot(is.character(kind), length(kind) == 1)
  obj <- structure(kind, class = c(kind, "character"))
  UseMethod("temporal_reconstruct", obj)
}

#' @export
temporal_reconstruct.default <- function(kind, basis, coeff, ...) {
  basis %*% coeff
}

#' @export
temporal_reconstruct.bspline <- function(kind, basis, coeff, ...) {
  basis %*% coeff
} 