# Bridge: shared-basis LNA archive -> fmrilatent LatentNeuroVec
#
# For archives written with embed.shared_basis, reconstructs a
# LatentNeuroVec where basis = per-subject coefficients and
# loadings = the shared spatial basis.

#' Coerce a shared-basis LNA archive to a LatentNeuroVec
#'
#' Reads a `DataHandle` produced by `read_lna()` on a shared-basis
#' archive and constructs a `fmrilatent::LatentNeuroVec`. The shared
#' basis is resolved via the standard resolution chain.
#'
#' @param x A `DataHandle` from `read_lna()` on an archive written
#'   with the `embed.shared_basis` transform.
#' @param ... Additional arguments (unused).
#' @return A `fmrilatent::LatentNeuroVec`.
#' @export
as_latent_neurovec <- function(x, ...) {
  UseMethod("as_latent_neurovec")
}

#' @export
as_latent_neurovec.DataHandle <- function(x, ...) {
  if (!requireNamespace("fmrilatent", quietly = TRUE)) {
    abort_lna("fmrilatent package is required for as_latent_neurovec()",
              .subclass = "lna_error_missing_package",
              location = "as_latent_neurovec.DataHandle")
  }
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    abort_lna("neuroim2 package is required for as_latent_neurovec()",
              .subclass = "lna_error_missing_package",
              location = "as_latent_neurovec.DataHandle")
  }

  # Find the embed.shared_basis descriptor
  desc <- .find_shared_basis_descriptor(x$h5)
  if (is.null(desc)) {
    abort_lna(
      "No embed.shared_basis transform found in this archive",
      .subclass = "lna_error_descriptor",
      location = "as_latent_neurovec.DataHandle"
    )
  }

  p <- desc$params %||% list()

  # Resolve the shared basis
  basis_data <- read_shared_basis(
    resolve_basis(
      basis_id = p$basis_id,
      checksum = p$basis_checksum,
      locator = p$basis_locator,
      archive_path = x$h5$filename
    )
  )
  loadings <- basis_data$loadings  # voxels x atoms

  # Get coefficients from the stash
  coeff <- x$stash$input
  if (is.null(coeff)) {
    # Try the coefficients key if input was not produced
    coeff <- x$stash$coefficients
  }
  if (is.null(coeff)) {
    abort_lna("No coefficient data found in stash",
              .subclass = "lna_error_missing_data",
              location = "as_latent_neurovec.DataHandle")
  }
  coeff <- as_dense_mat(coeff)  # time x atoms

  # Read offset if available
  offset <- numeric(0)
  offset_entries <- Filter(function(d) identical(d$role, "offset"), desc$datasets)
  if (length(offset_entries) > 0) {
    root <- x$h5[["/"]]
    offset <- tryCatch(
      h5_read(root, offset_entries[[1]]$path),
      error = function(e) numeric(0)
    )
    safe_h5_close(root)
  }

  # Build NeuroSpace and mask
  header <- .read_lna_header(x$h5, run_id = x$current_run_id)
  spatial_dims <- header$dims[1:3]
  spacing <- header$spacing %||% rep(1, 3)
  origin <- header$origin %||% rep(0, 3)
  n_time <- nrow(coeff)

  space4d <- neuroim2::NeuroSpace(c(spatial_dims, n_time),
                                   spacing = spacing, origin = origin)

  mask_arr <- basis_data$mask
  dim(mask_arr) <- spatial_dims
  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_arr,
    neuroim2::NeuroSpace(spatial_dims, spacing = spacing, origin = origin)
  )

  # Construct LatentNeuroVec: basis = coefficients (time x atoms),
  # loadings = shared spatial dictionary (voxels x atoms)
  fmrilatent::LatentNeuroVec(
    basis = Matrix::Matrix(coeff, sparse = FALSE),
    loadings = if (inherits(loadings, "Matrix")) loadings else Matrix::Matrix(loadings, sparse = FALSE),
    space = space4d,
    mask = mask_vol,
    offset = if (length(offset) > 0) as.numeric(offset) else numeric(0),
    label = p$basis_id %||% "",
    meta = list(
      family = "shared_basis",
      basis_id = p$basis_id,
      basis_kind = p$basis_kind %||% basis_data$basis_kind
    )
  )
}

#' @export
as_latent_neurovec.lna_reader <- function(x, ...) {
  handle <- x$data(...)
  as_latent_neurovec(handle)
}

# ---------------------------------------------------------------------------
# Internal: find embed.shared_basis descriptor in an HDF5 file
# ---------------------------------------------------------------------------

.find_shared_basis_descriptor <- function(h5) {
  if (!h5$exists("transforms")) return(NULL)
  tf_group <- h5[["transforms"]]
  on.exit(safe_h5_close(tf_group), add = TRUE)

  transforms <- discover_transforms(tf_group)
  if (nrow(transforms) == 0) return(NULL)

  idx <- which(transforms$type == "embed.shared_basis")
  if (length(idx) == 0) return(NULL)

  read_json_descriptor(tf_group, transforms$name[idx[1]])
}
