#' Minimum input dimensionality for a transform
#'
#' Provides the minimum number of dimensions required by a transform's
#' forward step. Packages can define methods for their own transforms.
#' The default requirement is 3 dimensions.
#'
#' @param type Character transform type.
#' @export
transform_min_dims <- function(type) {
  UseMethod("transform_min_dims", type)
}

#' @export
transform_min_dims.default <- function(type) {
  3L
}

#' @export
transform_min_dims.quant <- function(type) {
  1L
}

#' @export
transform_min_dims.basis <- function(type) {
  2L
}

#' @export
transform_min_dims.embed <- function(type) {
  2L # Embed also typically works on 2D (time x features) or (time x voxels)
}
