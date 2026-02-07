#' Coerce object to a dense matrix
#'
#' Provides a simple S3 generic used by write adapters to obtain
#' a time-by-voxel matrix representation. Packages can implement
#' methods for their own classes.
#'
#' @param obj Object to coerce.
#' @return A matrix.
#' @keywords internal
#' @export
as_dense_mat <- function(obj) {
  UseMethod("as_dense_mat")
}

#' @export
#' @keywords internal
as_dense_mat.default <- function(obj) {
  if (is.matrix(obj)) {
    obj
  } else {
    as.matrix(obj)
  }
}

#' @export
#' @keywords internal
as_dense_mat.array <- function(obj) {
  d <- dim(obj)
  if (length(d) <= 2) {
    return(as.matrix(obj))
  }
  time_dim <- d[length(d)]
  vox_dim <- prod(d[-length(d)])
  mat <- matrix(as.numeric(aperm(obj, c(length(d), seq_len(length(d) - 1)))),
                nrow = time_dim, ncol = vox_dim)
  mat
}

#' Coerce object to a 4-D array
#'
#' Provides a simple S3 generic used by write adapters when
#' volumetric data is required.
#'
#' @param obj Object to coerce.
#' @return A 4-D array.
#' @keywords internal
#' @export
as_vol4d <- function(obj) {
  UseMethod("as_vol4d")
}

#' @export
#' @keywords internal
as_vol4d.default <- function(obj) {
  if (is.array(obj) && length(dim(obj)) == 4) {
    obj
  } else {
    abort_lna(
      "cannot coerce object to 4-D array",
      .subclass = "lna_error_validation",
      location = "as_vol4d.default"
    )
  }
}

#' Convert data to a MaskedVoxels x Time matrix
#'
#' Helper used by spatial transforms. Accepts common representations
#' such as 4-D arrays (X\times Y\times Z\times Time), matrices in either
#' voxel-by-time or time-by-voxel orientation, and `neuroim2::NeuroVec`
#' objects when the \code{neuroim2} package is available. The result is
#' always a numeric matrix with rows corresponding to voxels inside the
#' supplied mask and columns corresponding to time points.
#'
#' @param x Input data object.
#' @param mask_3d_array Logical 3-D array defining the voxel mask.
#' @return Numeric matrix of dimension MaskedVoxels x Time.
#' @keywords internal
convert_to_masked_vox_time_matrix <- function(x, mask_3d_array) {
  if (is.null(mask_3d_array)) {
    abort_lna("mask_3d_array required",
              .subclass = "lna_error_validation",
              location = "convert_to_masked_vox_time_matrix:mask")
  }

  mask_idx <- which(as.logical(mask_3d_array))
  nvox <- prod(dim(mask_3d_array))
  nmask <- length(mask_idx)

  # neuroim2 object support if package is available
  if (inherits(x, "NeuroVec")) {
    if (requireNamespace("neuroim2", quietly = TRUE)) {
      mat <- as.matrix(neuroim2::series(x, mask_3d_array))
      return(t(mat))
    } else {
      abort_lna("neuroim2 package required for NeuroVec input",
                .subclass = "lna_error_validation",
                location = "convert_to_masked_vox_time_matrix:neuroim2")
    }
  }

  if (is.array(x) && length(dim(x)) == 4) {
    d <- dim(x)
    spatial_n <- prod(d[1:3])
    if (!identical(as.integer(spatial_n), as.integer(nvox))) {
      abort_lna(
        "mask_3d_array dimensions do not match x spatial dimensions",
        .subclass = "lna_error_validation",
        location = "convert_to_masked_vox_time_matrix:mask_dims"
      )
    }

    time_n <- d[4]
    if (nmask == 0) {
      return(matrix(numeric(0), nrow = 0, ncol = time_n))
    }

    # Low-allocation extraction: avoid building a full (Voxels x Time) matrix.
    offsets <- spatial_n * rep(0:(time_n - 1L), each = nmask)
    idx_vec <- rep.int(mask_idx, times = time_n) + offsets
    out <- matrix(as.numeric(x[idx_vec]), nrow = nmask, ncol = time_n)
    return(out)
  }

  mat <- as_dense_mat(x)

  if (nrow(mat) == nvox) {
    return(mat[mask_idx, , drop = FALSE])
  }

  if (ncol(mat) == nvox) {
    return(t(mat)[mask_idx, , drop = FALSE])
  }

  # Already masked: Voxel x Time or Time x Voxel
  if (nrow(mat) == nmask) {
    return(mat)
  }
  if (ncol(mat) == nmask) {
    return(t(mat))
  }

  # Fallback: return as-is
  mat
}
