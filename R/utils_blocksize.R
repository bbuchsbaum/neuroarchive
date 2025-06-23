#' Automatic Block Size for Spatial Slabs
#'
#' Determines spatial slab dimensions for block-wise processing. The
#' returned slab is chosen so that the estimated memory footprint does
#' not exceed `target_slab_bytes`. The depth (Z dimension) is reduced
#' first, then the Y and X dimensions as needed.
#'
#' @param spatial_dims Integer vector of length 3 giving the X, Y, Z
#'   dimensions.
#' @param element_size_bytes Size in bytes of a single element.
#' @param target_slab_bytes Target maximum size of a slab in bytes.
#' @return A list with `slab_dims` (integer vector of length 3) and
#'   `iterate_slabs` (number of slabs along each dimension).
#' @keywords internal
auto_block_size <- function(spatial_dims, element_size_bytes,
                            target_slab_bytes = NULL) {
  if (is.null(target_slab_bytes)) {
    target_slab_bytes <- lna_options("memory.target_slab_bytes")[[1]]
  }
  stopifnot(is.numeric(spatial_dims), length(spatial_dims) == 3)
  stopifnot(is.numeric(element_size_bytes), length(element_size_bytes) == 1)
  dims <- pmax(as.integer(spatial_dims), 1L)
  size <- as.numeric(element_size_bytes)

  slab <- dims
  slab[3] <- min(dims[3], 1L)

  bytes <- prod(slab) * size
  while (bytes > target_slab_bytes && slab[2] > 1) {
    slab[2] <- ceiling(slab[2] / 2)
    bytes <- prod(slab) * size
  }
  while (bytes > target_slab_bytes && slab[1] > 1) {
    slab[1] <- ceiling(slab[1] / 2)
    bytes <- prod(slab) * size
  }

  slab <- pmax(slab, 1L)
  iterate <- ceiling(dims / slab)
  list(slab_dims = slab, iterate_slabs = iterate)
}

