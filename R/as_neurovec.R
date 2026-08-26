# Coercion bridge: DataHandle / lna_reader -> neuroim2 DenseNeuroVec
#
# Reconstructs a spatially-aware NeuroVec from LNA read results by
# recovering the NeuroSpace from run-scoped header attributes (or the legacy
# single-run global header) and reshaping the flat stash matrix back to 4D.

#' Coerce an LNA read result to a NeuroVec
#'
#' Converts the flat matrix in a `DataHandle` (from `read_lna()`) into a
#' `neuroim2::DenseNeuroVec` with full spatial metadata. The spatial
#' information is recovered from `/header/runs/<run_id>` when present, with
#' `/header/global` used only as the legacy single-run compatibility path.
#'
#' @param x A `DataHandle` (from `read_lna()`) or an `lna_reader`
#'   (from `read_lna(..., lazy = TRUE)`), or a file path to an LNA file.
#' @param assume_identity Logical. If `FALSE` (the default), missing spatial
#'   dimensions or a missing affine are errors. Set to `TRUE` only for generic
#'   non-fMRI archives whose geometry is intentionally unknown; the returned
#'   object is then marked with `lna_geometry_provenance`.
#' @param ... Additional arguments. For `DataHandle`, pass `file` to
#'   specify the LNA file path if the HDF5 handle is closed.
#' @return A `neuroim2::DenseNeuroVec` or `neuroim2::SparseNeuroVec`.
#' @export
as_neurovec <- function(x, assume_identity = FALSE, ...) {
  UseMethod("as_neurovec")
}

#' @export
as_neurovec.DataHandle <- function(x, assume_identity = FALSE, file = NULL, ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    abort_lna("neuroim2 package is required for as_neurovec()",
              .subclass = "lna_error_missing_package",
              location = "as_neurovec.DataHandle")
  }

  header <- .read_lna_header_safe(x$h5, file, run_id = x$current_run_id)
  .as_neurovec_from_header(
    x,
    header,
    assume_identity = assume_identity,
    location = "as_neurovec.DataHandle"
  )
}

#' @export
as_neurovec.lna_reader <- function(x, assume_identity = FALSE, ...) {
  handle <- x$data(...)
  # The reader keeps h5 open, so header reading works via its handle
  .as_neurovec_with_h5(handle, x$h5, assume_identity = assume_identity)
}

#' @export
as_neurovec.character <- function(x, assume_identity = FALSE, ...) {
  handle <- read_lna(x, ...)
  as_neurovec(handle, assume_identity = assume_identity, file = x)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build DenseNeuroVec from a DataHandle using an open H5 for header
#' @keywords internal
.as_neurovec_with_h5 <- function(handle, h5, assume_identity = FALSE) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    abort_lna("neuroim2 package is required",
              .subclass = "lna_error_missing_package",
              location = ".as_neurovec_with_h5")
  }

  header <- .read_lna_header(h5, run_id = handle$current_run_id)
  .as_neurovec_from_header(
    handle,
    header,
    assume_identity = assume_identity,
    location = ".as_neurovec_with_h5"
  )
}

#' Resolve stored spatial identity into 3D and 4D NeuroSpace objects
#' @keywords internal
.lna_neurospaces_from_header <- function(header, inferred_dims, n_time,
                                         assume_identity = FALSE,
                                         location = "as_neurovec") {
  stored_dims <- header$dims
  has_dims <- !is.null(stored_dims) && length(stored_dims) >= 3L

  if (has_dims) {
    spatial_dims <- suppressWarnings(as.integer(stored_dims[1:3]))
    valid_dims <- is.numeric(stored_dims) &&
      length(spatial_dims) == 3L &&
      all(is.finite(stored_dims[1:3])) &&
      all(stored_dims[1:3] == spatial_dims) &&
      all(spatial_dims > 0L)
    if (!isTRUE(valid_dims)) {
      abort_lna(
        "Stored spatial dimensions must be three positive integers",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
    if (!is.null(inferred_dims) && length(inferred_dims) >= 3L &&
        !identical(spatial_dims, as.integer(inferred_dims[1:3]))) {
      abort_lna(
        "Stored spatial dimensions do not match reconstructed data",
        stored_dims = spatial_dims,
        data_dims = as.integer(inferred_dims[1:3]),
        .subclass = "lna_error_geometry",
        location = location
      )
    }
  } else if (isTRUE(assume_identity)) {
    if (is.null(inferred_dims) || length(inferred_dims) < 3L) {
      abort_lna(
        "Spatial dimensions cannot be inferred for identity geometry",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
    spatial_dims <- as.integer(inferred_dims[1:3])
    valid_dims <- all(is.finite(inferred_dims[1:3])) &&
      all(inferred_dims[1:3] == spatial_dims) &&
      all(spatial_dims > 0L)
    if (!isTRUE(valid_dims)) {
      abort_lna(
        "Inferred spatial dimensions must be three positive integers",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
  } else {
    abort_lna(
      "LNA spatial dimensions are missing; fMRI reconstruction requires stored geometry",
      .subclass = "lna_error_geometry",
      location = location
    )
  }

  stored_transform <- header$transform
  has_transform <- !is.null(stored_transform)
  if (has_transform) {
    if (!is.numeric(stored_transform) ||
        !identical(dim(stored_transform), c(4L, 4L)) ||
        any(!is.finite(stored_transform))) {
      abort_lna(
        "Stored affine must be a finite numeric 4 x 4 matrix",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
    affine <- stored_transform
    provenance <- list(
      schema_version = "1.0.0",
      source = "lna_header",
      assumed_identity = FALSE
    )
  } else if (isTRUE(assume_identity)) {
    affine <- diag(4L)
    provenance <- list(
      schema_version = "1.0.0",
      source = "explicit_identity_assumption",
      assumed_identity = TRUE
    )
  } else {
    abort_lna(
      "LNA affine is missing; fMRI reconstruction will not fabricate identity geometry",
      .subclass = "lna_error_geometry",
      location = location
    )
  }

  n_time_int <- suppressWarnings(as.integer(n_time))
  if (length(n_time_int) != 1L || !is.finite(n_time) ||
      n_time != n_time_int || n_time_int < 1L) {
    abort_lna(
      "Reconstructed time dimension must be one positive integer",
      .subclass = "lna_error_geometry",
      location = location
    )
  }

  spaces <- tryCatch(
    list(
      space3d = neuroim2::NeuroSpace(spatial_dims, trans = affine),
      space4d = neuroim2::NeuroSpace(c(spatial_dims, n_time_int),
                                     trans = affine)
    ),
    error = function(e) {
      abort_lna(
        "Stored affine cannot define a valid NeuroSpace",
        .subclass = "lna_error_geometry",
        location = location,
        parent = e
      )
    }
  )
  c(spaces, list(spatial_dims = spatial_dims, provenance = provenance))
}

#' Build a NeuroVec from a DataHandle and resolved header
#' @keywords internal
.as_neurovec_from_header <- function(handle, header, assume_identity = FALSE,
                                     location = "as_neurovec") {
  raw_data <- handle$stash$input
  if (is.null(raw_data)) {
    abort_lna(
      "No reconstructed data in stash (key 'input' is NULL)",
      .subclass = "lna_error_missing_data",
      location = location
    )
  }

  if (is.array(raw_data) && length(dim(raw_data)) == 4L) {
    data_dims <- dim(raw_data)
    resolved <- .lna_neurospaces_from_header(
      header,
      inferred_dims = data_dims[1:3],
      n_time = data_dims[4],
      assume_identity = assume_identity,
      location = location
    )
    result <- neuroim2::DenseNeuroVec(raw_data, resolved$space4d)
    attr(result, "lna_geometry_provenance") <- resolved$provenance
    return(result)
  }

  mat <- as_dense_mat(raw_data)
  roi_mask <- handle$subset$roi_mask %||% handle$subset$roi
  archive_mask <- handle$mask_info$mask
  mask_arr <- archive_mask

  if (!is.null(roi_mask)) {
    if (!is.null(archive_mask)) {
      archive_values <- as.logical(archive_mask)
      roi_values <- as.logical(roi_mask)
      if (length(archive_values) != length(roi_values)) {
        abort_lna(
          "ROI mask length does not match the stored archive mask",
          archive_mask_length = length(archive_values),
          roi_mask_length = length(roi_values),
          .subclass = "lna_error_geometry",
          location = location
        )
      }
      mask_dims <- dim(archive_mask) %||% dim(roi_mask)
      mask_arr <- archive_values & roi_values
      if (!is.null(mask_dims)) dim(mask_arr) <- mask_dims
    } else {
      mask_arr <- roi_mask
    }
  }

  if (!is.null(mask_arr)) {
    n_roi <- sum(as.logical(mask_arr))
    if (nrow(mat) == n_roi && ncol(mat) != n_roi) {
      mat <- t(mat)
    }
  }

  inferred_dims <- NULL
  if (!is.null(mask_arr) && !is.null(dim(mask_arr))) {
    inferred_dims <- dim(mask_arr)
  } else if (isTRUE(assume_identity)) {
    inferred_dims <- c(ncol(mat), 1L, 1L)
  }

  resolved <- .lna_neurospaces_from_header(
    header,
    inferred_dims = inferred_dims,
    n_time = nrow(mat),
    assume_identity = assume_identity,
    location = location
  )
  spatial_dims <- resolved$spatial_dims
  n_time <- nrow(mat)
  n_vox <- ncol(mat)
  n_spatial <- prod(spatial_dims)

  if (!is.null(mask_arr)) {
    if (length(mask_arr) != n_spatial) {
      abort_lna(
        "Stored mask length does not match spatial dimensions",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
    mask_arr <- array(as.logical(mask_arr), dim = spatial_dims)
  }

  if (n_vox < n_spatial) {
    if (is.null(mask_arr) || sum(mask_arr) != n_vox) {
      abort_lna(
        "Masked reconstruction requires a mask matching the stored voxel count",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
    mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, resolved$space3d)
    result <- neuroim2::SparseNeuroVec(
      data = t(mat),
      space = resolved$space4d,
      mask = mask_vol
    )
  } else {
    if (n_vox != n_spatial) {
      abort_lna(
        "Reconstructed voxel count does not match spatial dimensions",
        .subclass = "lna_error_geometry",
        location = location
      )
    }
    arr4d <- array(
      as.numeric(aperm(
        array(mat, dim = c(n_time, spatial_dims)),
        c(2, 3, 4, 1)
      )),
      dim = c(spatial_dims, n_time)
    )
    result <- neuroim2::DenseNeuroVec(arr4d, resolved$space4d)
  }

  attr(result, "lna_geometry_provenance") <- resolved$provenance
  result
}

#' Read the LNA header from an open HDF5 handle
#' @keywords internal
.read_lna_header_group <- function(group) {
  header <- list()
  fields <- c(
    "geometry_schema_version", "run_id", "dims", "spacing", "origin",
    "transform", "spatial_units", "voxel_ordering", "n_volumes", "grid_id"
  )
  for (nm in fields) {
    if (group$attr_exists(nm)) {
      attr_obj <- group$attr_open(nm)
      header[[nm]] <- attr_obj$read()
      safe_h5_close(attr_obj)
    }
  }
  header
}

.read_lna_header <- function(h5, run_id = NULL) {
  header <- list()
  if (!h5$exists("header")) return(header)
  hdr_grp <- h5[["header"]]
  on.exit(safe_h5_close(hdr_grp), add = TRUE)

  if (hdr_grp$exists("runs")) {
    runs_group <- hdr_grp[["runs"]]
    on.exit(safe_h5_close(runs_group), add = TRUE)
    available <- names(runs_group)
    selected <- run_id
    if (is.null(selected)) {
      if (length(available) == 1L) {
        selected <- available[[1L]]
      } else if (length(available) > 1L) {
        abort_lna(
          paste0(
            "run-scoped geometry is ambiguous; specify one of: ",
            paste(available, collapse = ", ")
          ),
          runs = available,
          .subclass = "lna_error_run_geometry",
          location = ".read_lna_header"
        )
      }
    }
    if (!is.null(selected)) {
      if (!is.character(selected) || length(selected) != 1L ||
          is.na(selected) || !nzchar(selected)) {
        abort_lna(
          "run_id must be one non-empty string when reading geometry",
          .subclass = "lna_error_run_geometry",
          location = ".read_lna_header"
        )
      }
      if (!runs_group$exists(selected)) {
        abort_lna(
          paste0("run geometry is missing for ", selected),
          run_id = selected,
          .subclass = "lna_error_run_geometry",
          location = ".read_lna_header"
        )
      }
      group <- runs_group[[selected]]
      on.exit(safe_h5_close(group), add = TRUE)
      header <- .read_lna_header_group(group)
      attr(header, "lna_geometry_source") <- "run-scoped"
      return(header)
    }
  }

  if (!hdr_grp$exists("global")) return(header)
  g <- hdr_grp[["global"]]
  on.exit(safe_h5_close(g), add = TRUE)
  header <- .read_lna_header_group(g)
  attr(header, "lna_geometry_source") <- "legacy-global"
  header
}

#' Read header safely - try open handle first, fall back to file path
#' @keywords internal
.read_lna_header_safe <- function(h5, file_path = NULL, run_id = NULL) {
  # Try the handle if it's valid
  if (!is.null(h5)) {
    valid <- tryCatch(h5$is_valid, error = function(e) FALSE)
    if (isTRUE(valid)) {
      return(.read_lna_header(h5, run_id = run_id))
    }
  }

  # Fall back to file path
  if (!is.null(file_path) && is.character(file_path) && file.exists(file_path)) {
    h5_new <- hdf5r::H5File$new(file_path, mode = "r")
    on.exit(try(h5_new$close_all(), silent = TRUE), add = TRUE)
    return(.read_lna_header(h5_new, run_id = run_id))
  }

  # Try to extract filename from the closed handle
  file_from_h5 <- tryCatch(h5$filename, error = function(e) NULL)
  if (!is.null(file_from_h5) && file.exists(file_from_h5)) {
    h5_new <- hdf5r::H5File$new(file_from_h5, mode = "r")
    on.exit(try(h5_new$close_all(), silent = TRUE), add = TRUE)
    return(.read_lna_header(h5_new, run_id = run_id))
  }

  list()
}
