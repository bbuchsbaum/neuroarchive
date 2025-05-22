#' Minimum input dimensionality for a transform
#'
#' Provides the minimum number of dimensions required by a transform's
#' forward step. Packages can define methods for their own transforms.
#' The default requirement is 3 dimensions.
#'
#' @param type Character transform type.
#' @keywords internal
transform_min_dims <- function(type) {
  UseMethod("transform_min_dims", type)
}

#' @export
#' @keywords internal
transform_min_dims.default <- function(type) {
  3L
}

#' @export
#' @keywords internal
transform_min_dims.quant <- function(type) {
  1L
}
