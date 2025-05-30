#' Placeholder synthesis for Haar octwave coefficients
#'
#' This helper is a stub and simply returns a zero matrix of the
#' appropriate dimensions. It will be replaced by a real implementation
#' of the inverse lifting scheme.
#' @keywords internal
perform_haar_lift_synthesis <- function(coeff_list, mask_3d_array, levels) {
  n_vox <- sum(mask_3d_array)
  if (is.null(coeff_list$root)) return(matrix(0, 0, n_vox))
  n_time <- nrow(coeff_list$root)
  matrix(0, nrow = n_time, ncol = n_vox)
}

#' Morton-ordered voxel indices from a 3D mask
#'
#' Returns the 1D linear indices of voxels inside a logical 3D mask
#' ordered by the Morton (Z-order) curve. A seed can be supplied to
#' deterministically break ties if the Morton code generation yields
#' identical codes for some voxels.
#'
#' @param mask_3d_array Logical 3D array representing the mask.
#' @param z_order_seed Integer seed used for tie-breaking. Default is 42.
#' @return Integer vector of 1-based linear indices in Morton order.
#' @keywords internal
get_morton_ordered_indices <- function(mask_3d_array, z_order_seed = 42L) {
  if (!is.array(mask_3d_array) || length(dim(mask_3d_array)) != 3L) {
    abort_lna("mask_3d_array must be a 3D array",
              .subclass = "lna_error_validation",
              location = "get_morton_ordered_indices:mask")
  }

  mask_logical <- as.logical(mask_3d_array)
  vox_coords <- which(mask_logical, arr.ind = TRUE)
  if (nrow(vox_coords) == 0L) {
    return(integer(0))
  }

  dims <- dim(mask_logical)
  x <- vox_coords[, 1] - 1L
  y <- vox_coords[, 2] - 1L
  z <- vox_coords[, 3] - 1L

  max_dim <- max(dims)
  bits <- ceiling(log2(max_dim))
  codes <- integer(nrow(vox_coords))
  for (b in seq_len(bits)) {
    shift <- b - 1L
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(x, shift), 1L), 3L * shift))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(y, shift), 1L), 3L * shift + 1L))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(z, shift), 1L), 3L * shift + 2L))
  }

  set.seed(z_order_seed)
  tie_break <- runif(length(codes))
  ordering <- order(codes, tie_break)
  ordered_coords <- vox_coords[ordering, , drop = FALSE]

  linear_idx <- (ordered_coords[, 1] - 1L) +
    (ordered_coords[, 2] - 1L) * dims[1L] +
    (ordered_coords[, 3] - 1L) * dims[1L] * dims[2L] + 1L

  as.integer(linear_idx)
}

#' SHA1 hash for a sequence of Morton-ordered indices
#'
#' Computes a SHA1 hash of the integer vector representing Morton-ordered
#' voxel indices. The resulting string is prefixed with ``"sha1:"`` to
#' mirror other hash fields within the package.
#'
#' @param ordered_indices_vector Integer vector of indices from
#'   `get_morton_ordered_indices`.
#' @return Character scalar ``"sha1:<hash>"``.
#' @keywords internal
morton_indices_to_hash <- function(ordered_indices_vector) {
  hash_val <- digest::digest(ordered_indices_vector, algo = "sha1",
                             serialize = TRUE)
  paste0("sha1:", hash_val)
}
