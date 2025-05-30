#' Perform inverse Haar lifting synthesis
#'
#' Reconstructs a Time x MaskedVoxels matrix from Haar octwave
#' coefficients as produced by `perform_haar_lift_analysis`.
#'
#' @param coeff_list List with elements `root` (matrix `Time x 1`) and
#'   `detail` (list of length `levels`, each a `Time x Nvox` matrix).
#' @param mask_3d_array Logical 3D array defining the voxel mask.
#' @param levels Number of decomposition levels used.
#' @param z_order_seed Seed used for Morton ordering (must match the one
#'   used during analysis). Default is 42.
#' @return Matrix of reconstructed values with dimensions `Time x
#'   MaskedVoxels`.
#' @keywords internal
perform_haar_lift_synthesis <- function(coeff_list, mask_3d_array, levels,
                                         z_order_seed = 42L) {
  if (is.null(coeff_list$root) || !is.matrix(coeff_list$root)) {
    abort_lna("coeff_list$root must be a matrix",
              .subclass = "lna_error_validation",
              location = "perform_haar_lift_synthesis:root")
  }

  if (!is.list(coeff_list$detail) || length(coeff_list$detail) < levels) {
    abort_lna("coeff_list$detail must be a list of length 'levels'",
              .subclass = "lna_error_validation",
              location = "perform_haar_lift_synthesis:detail")
  }

  morton_idx <- get_morton_ordered_indices(mask_3d_array, z_order_seed)
  mask_linear <- which(as.logical(mask_3d_array))
  perm <- match(morton_idx, mask_linear)
  inv_perm <- order(perm)

  full_order <- get_morton_ordered_indices(array(TRUE, dim(mask_3d_array)),
                                           z_order_seed)
  mask_flat <- as.logical(mask_3d_array)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- precompute_haar_scalings(mask_3d_array, levels)

  n_time <- nrow(coeff_list$root)
  n_vox <- length(morton_idx)
  reco_morton <- matrix(0, nrow = n_time, ncol = n_vox)

  use_rcpp <- isTRUE(getOption("lna.hwt.use_rcpp", TRUE)) &&
    exists("inverse_lift_rcpp")

  for (tt in seq_len(n_time)) {
    detail_vecs <- lapply(seq_len(levels),
                          function(lvl) coeff_list$detail[[lvl]][tt, ])
    if (use_rcpp) {
      reco_morton[tt, ] <- inverse_lift_rcpp(
        coeff_list$root[tt, 1],
        detail_vecs,
        mask_flat_morton,
        dim(mask_3d_array),
        levels,
        scalings
      )
    } else {
      reco_morton[tt, ] <- inverse_lift_R(
        coeff_list$root[tt, 1],
        detail_vecs,
        mask_flat_morton,
        dim(mask_3d_array),
        levels,
        scalings
      )
    }
  }

  reco_morton[, inv_perm, drop = FALSE]
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
  use_rcpp <- isTRUE(getOption("lna.hwt.use_rcpp", TRUE)) &&
    exists("morton_indices_to_hash_rcpp")
  if (use_rcpp) {
    return(morton_indices_to_hash_rcpp(as.integer(ordered_indices_vector)))
  }
  hash_val <- digest::digest(ordered_indices_vector, algo = "sha1",
                             serialize = TRUE)
  paste0("sha1:", hash_val)
}

#' Precompute scaling factors for Haar lifting
#'
#' Returns a list with one element per decomposition level. Each element
#' contains numeric vectors `sqrt_nvalid` and `sqrt_nvalid_div_8`
#' representing the $\sqrt{n_{valid}}$ and $\sqrt{n_{valid}/8}$ factors
#' for all conceptual $2 \times 2 \times 2$ blocks in Morton order.
#' @keywords internal
precompute_haar_scalings <- function(mask_3d_array, levels) {
  if (!is.array(mask_3d_array) || length(dim(mask_3d_array)) != 3L) {
    abort_lna("mask_3d_array must be a 3D array",
              .subclass = "lna_error_validation",
              location = "precompute_haar_scalings:mask")
  }

  mask_logical <- as.logical(mask_3d_array)
  scalings <- vector("list", levels)
  current <- mask_logical
  for (lvl in seq_len(levels)) {
    dims <- dim(current)
    x_seq <- seq(1L, dims[1], by = 2L)
    y_seq <- seq(1L, dims[2], by = 2L)
    z_seq <- seq(1L, dims[3], by = 2L)
    counts <- integer(length(x_seq) * length(y_seq) * length(z_seq))
    idx <- 1L
    for (x in x_seq) {
      for (y in y_seq) {
        for (z in z_seq) {
          block <- current[
            x:min(x + 1L, dims[1]),
            y:min(y + 1L, dims[2]),
            z:min(z + 1L, dims[3])
          ]
          counts[idx] <- sum(block)
          idx <- idx + 1L
        }
      }
    }
    scalings[[lvl]] <- list(
      sqrt_nvalid = sqrt(as.numeric(counts)),
      sqrt_nvalid_div_8 = sqrt(as.numeric(counts) / 8)
    )
    current <- array(counts > 0,
                     dim = c(length(x_seq), length(y_seq), length(z_seq)))
  }
  scalings
}

#' Morton codes of valid finest-level blocks
#'
#' Identifies all 2x2x2 blocks at the finest decomposition level that
#' contain at least one in-mask voxel and returns their Morton codes as
#' 0-based unsigned integers.
#'
#' @param mask_3d_array Logical 3D array defining the voxel mask.
#' @return Integer vector of Morton codes for valid finest blocks.
#' @keywords internal
get_valid_finest_blocks <- function(mask_3d_array) {
  if (!is.array(mask_3d_array) || length(dim(mask_3d_array)) != 3L) {
    abort_lna("mask_3d_array must be a 3D array",
              .subclass = "lna_error_validation",
              location = "get_valid_finest_blocks:mask")
  }

  mask_logical <- as.logical(mask_3d_array)
  dims <- dim(mask_logical)
  x_seq <- seq(1L, dims[1], by = 2L)
  y_seq <- seq(1L, dims[2], by = 2L)
  z_seq <- seq(1L, dims[3], by = 2L)
  block_mask <- array(FALSE,
                      dim = c(length(x_seq), length(y_seq), length(z_seq)))

  for (i in seq_along(x_seq)) {
    for (j in seq_along(y_seq)) {
      for (k in seq_along(z_seq)) {
        block <- mask_logical[
          x_seq[i]:min(x_seq[i] + 1L, dims[1]),
          y_seq[j]:min(y_seq[j] + 1L, dims[2]),
          z_seq[k]:min(z_seq[k] + 1L, dims[3])
        ]
        block_mask[i, j, k] <- any(block)
      }
    }
  }

  if (!any(block_mask)) return(integer(0))

  coords <- which(block_mask, arr.ind = TRUE)
  x <- coords[, 1] - 1L
  y <- coords[, 2] - 1L
  z <- coords[, 3] - 1L

  max_dim <- max(dim(block_mask))
  bits <- ceiling(log2(max_dim))
  codes <- integer(nrow(coords))
  for (b in seq_len(bits)) {
    shift <- b - 1L
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(x, shift), 1L), 3L * shift))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(y, shift), 1L),
                                3L * shift + 1L))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(z, shift), 1L),
                                3L * shift + 2L))
  }

  sort(as.integer(codes))
}

#' Pure-R forward Haar lifting implementation
#'
#' @keywords internal
forward_lift_R <- function(data_masked_morton_ordered,
                           mask_flat_morton_ordered,
                           mask_dims,
                           levels,
                           scaling_factors_per_level) {
  current <- data_masked_morton_ordered
  detail_list <- vector("list", levels)
  for (lvl in seq_len(levels)) {
    sc <- scaling_factors_per_level[[lvl]]
    counts <- as.integer(round(sc$sqrt_nvalid^2))
    next_data <- numeric(length(counts))
    dvec <- numeric(length(current))
    idx_in <- 1L
    idx_out <- 1L
    for (b in seq_along(counts)) {
      nv <- counts[b]
      if (nv > 0) {
        vals <- current[idx_in:(idx_in + nv - 1L)]
        avg <- mean(vals)
        next_data[b] <- avg * sc$sqrt_nvalid[b]
        dvec[idx_out:(idx_out + nv - 1L)] <- (vals - avg) * sc$sqrt_nvalid_div_8[b]
        idx_in <- idx_in + nv
        idx_out <- idx_out + nv
      }
    }
    detail_list[[lvl]] <- dvec
    current <- next_data
  }
  list(root_coeff = current[1], detail_coeffs_by_level = detail_list)
}

#' Pure-R inverse Haar lifting implementation
#'
#' @keywords internal
inverse_lift_R <- function(root_coeff,
                           detail_coeffs_by_level,
                           mask_flat_morton_ordered,
                           mask_dims,
                           levels,
                           scaling_factors_per_level) {
  current <- root_coeff
  for (lvl in rev(seq_len(levels))) {
    sc <- scaling_factors_per_level[[lvl]]
    counts <- as.integer(round(sc$sqrt_nvalid^2))
    next_data <- numeric(sum(counts))
    idx_out <- 1L
    idx_det <- 1L
    for (b in seq_along(counts)) {
      nv <- counts[b]
      if (nv > 0) {
        avg <- current[b] / sc$sqrt_nvalid[b]
        vals <- detail_coeffs_by_level[[lvl]][idx_det:(idx_det + nv - 1L)] /
          sc$sqrt_nvalid_div_8[b] + avg
        next_data[idx_out:(idx_out + nv - 1L)] <- vals
        idx_out <- idx_out + nv
        idx_det <- idx_det + nv
      }
    }
    current <- next_data
  }
  current
}

#' Perform forward Haar lifting analysis
#'
#' Orchestrates the recursive lifting for a full Time x MaskedVoxels
#' matrix. Columns of `data_matrix_T_x_Nmask` are reordered to Morton
#' order and each time point processed with `forward_lift_rcpp`.
#' @keywords internal
perform_haar_lift_analysis <- function(data_matrix_T_x_Nmask, mask_3d_array,
                                       levels, z_order_seed = 42L) {
  if (!is.matrix(data_matrix_T_x_Nmask)) {
    abort_lna("data_matrix_T_x_Nmask must be a matrix",
              .subclass = "lna_error_validation",
              location = "perform_haar_lift_analysis:data")
  }

  morton_idx <- get_morton_ordered_indices(mask_3d_array, z_order_seed)
  mask_linear <- which(as.logical(mask_3d_array))
  perm <- match(morton_idx, mask_linear)
  data_morton <- data_matrix_T_x_Nmask[, perm, drop = FALSE]

  full_order <- get_morton_ordered_indices(array(TRUE, dim(mask_3d_array)),
                                           z_order_seed)
  mask_flat <- as.logical(mask_3d_array)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- precompute_haar_scalings(mask_3d_array, levels)

  n_time <- nrow(data_morton)
  n_vox <- length(morton_idx)
  root_coeffs <- matrix(0, nrow = n_time, ncol = 1)
  detail_coeffs <- lapply(seq_len(levels),
                          function(i) matrix(0, nrow = n_time, ncol = n_vox))

  use_rcpp <- isTRUE(getOption("lna.hwt.use_rcpp", TRUE)) &&
    exists("forward_lift_rcpp")

  for (tt in seq_len(n_time)) {
    if (use_rcpp) {
      res <- forward_lift_rcpp(
        data_morton[tt, ],
        mask_flat_morton,
        dim(mask_3d_array),
        levels,
        scalings
      )
    } else {
      res <- forward_lift_R(
        data_morton[tt, ],
        mask_flat_morton,
        dim(mask_3d_array),
        levels,
        scalings
      )
    }
    root_coeffs[tt, 1] <- res$root_coeff
    for (lvl in seq_len(levels)) {
      detail_coeffs[[lvl]][tt, ] <- res$detail_coeffs_by_level[[lvl]]
    }
  }

  list(root = root_coeffs, detail = detail_coeffs)
}

#' Morton code from block coordinates
#'
#' Helper used for ROI streaming. Encodes 0-based block coordinates
#' `(x, y, z)` into an integer Morton code using `bits` bits per axis.
#' @keywords internal
.encode_morton3d <- function(x, y, z, bits) {
  code <- 0L
  for (b in seq_len(bits)) {
    shift <- b - 1L
    code <- bitwOr(code, bitwShiftL(bitwAnd(bitwShiftR(x, shift), 1L), 3L * shift))
    code <- bitwOr(code, bitwShiftL(bitwAnd(bitwShiftR(y, shift), 1L), 3L * shift + 1L))
    code <- bitwOr(code, bitwShiftL(bitwAnd(bitwShiftR(z, shift), 1L), 3L * shift + 2L))
  }
  as.integer(code)
}

#' Block map for Haar octwave coefficients
#'
#' Builds a mapping from block Morton codes to the column index ranges of
#' detail coefficients for each decomposition level. The mapping mirrors the
#' enumeration order used by `perform_haar_lift_analysis` and therefore allows
#' subsetting coefficient matrices during inverse reconstruction.
#'
#' @return List of length `levels` with elements `code`, `count` and `start`
#'   (0-based start indices).
#' @keywords internal
.compute_block_map <- function(mask_3d_array, levels) {
  mask_logical <- as.logical(mask_3d_array)
  mapping <- vector("list", levels)
  current <- mask_logical
  for (lvl in seq_len(levels)) {
    dims <- dim(current)
    x_seq <- seq(1L, dims[1], by = 2L)
    y_seq <- seq(1L, dims[2], by = 2L)
    z_seq <- seq(1L, dims[3], by = 2L)
    nbx <- length(x_seq); nby <- length(y_seq); nbz <- length(z_seq)
    nblocks <- nbx * nby * nbz
    counts <- integer(nblocks)
    codes <- integer(nblocks)
    bits <- ceiling(log2(max(nbx, nby, nbz)))
    idx <- 1L
    for (i in seq_along(x_seq)) {
      for (j in seq_along(y_seq)) {
        for (k in seq_along(z_seq)) {
          block <- current[
            x_seq[i]:min(x_seq[i] + 1L, dims[1]),
            y_seq[j]:min(y_seq[j] + 1L, dims[2]),
            z_seq[k]:min(z_seq[k] + 1L, dims[3])
          ]
          counts[idx] <- sum(block)
          codes[idx] <- .encode_morton3d(i - 1L, j - 1L, k - 1L, bits)
          idx <- idx + 1L
        }
      }
    }
    start <- c(0L, cumsum(counts)[-length(counts)])
    mapping[[lvl]] <- list(code = codes, count = counts, start = start)
    current <- array(counts > 0L, dim = c(nbx, nby, nbz))
  }
  mapping
}

#' Determine detail coefficient indices required for an ROI
#'
#' Given an ROI mask, the original mask and the number of decomposition levels,
#' returns a list of column index vectors for each level identifying which
#' detail coefficients are necessary for reconstruction. Indices are 1-based.
#' @keywords internal
get_roi_detail_indices <- function(roi_mask, mask_3d_array, levels) {
  if (!is.array(roi_mask) || !identical(dim(roi_mask), dim(mask_3d_array))) {
    abort_lna("roi_mask must match mask dimensions",
              .subclass = "lna_error_validation",
              location = "get_roi_detail_indices:roi")
  }

  mapping <- .compute_block_map(mask_3d_array, levels)
  roi_codes_finest <- get_valid_finest_blocks(roi_mask)
  if (length(roi_codes_finest) == 0L) {
    return(lapply(seq_len(levels), function(x) integer()))
  }

  codes_by_level <- vector("list", levels)
  codes_by_level[[levels]] <- roi_codes_finest
  if (levels > 1) {
    for (lvl in seq_len(levels - 1)) {
      shift <- 3L * (levels - lvl)
      codes_by_level[[lvl]] <- unique(roi_codes_finest %/% (2L ^ shift))
    }
  }

  indices_by_level <- vector("list", levels)
  for (lvl in seq_len(levels)) {
    codes <- codes_by_level[[lvl]]
    map <- mapping[[lvl]]
    pos <- match(codes, map$code, nomatch = 0L)
    keep <- pos > 0L & map$count[pos] > 0L
    start <- map$start[pos[keep]]
    cnt <- map$count[pos[keep]]
    idx <- unlist(mapply(function(s, c) seq.int(s + 1L, s + c), start, cnt,
                         SIMPLIFY = FALSE), use.names = FALSE)
    indices_by_level[[lvl]] <- idx
  }
  indices_by_level
}
