#' @include as_neurovec.R transform_capabilities.R
#' @importFrom neuroim2 NeuroSpace LogicalNeuroVol IndexLookupVol SparseNeuroVol SparseNeuroVec
#' @importFrom neuroim2 space spacing origin indices
#' @importClassesFrom neuroim2 AbstractSparseNeuroVec SparseNeuroVec NeuroVec ArrayLike4D NeuroSpace LogicalNeuroVol IndexLookupVol
#' @importMethodsFrom neuroim2 linear_access matricized_access series sub_vector lookup
NULL

# LNANeuroVec: lazy subset-dispatch NeuroVec backed by an LNA archive
#
# Inherits from SparseNeuroVec for S4 dispatch compatibility. The @data
# slot is a zero placeholder. All real data access goes through @.state
# which holds an lna_reader and dispatches subsetted reads through the
# inverse pipeline when the pipeline supports it.

#' LNANeuroVec Class
#'
#' A lazy 4D neuroimaging vector backed by an LNA archive. Spatial
#' and temporal accessors dispatch subsetted reads through the inverse
#' pipeline when all transforms in the pipeline support subsetting.
#' For pipelines that do not fully support ROI subsetting, the full
#' masked matrix is materialized on first spatial access (with a
#' warning for large datasets).
#'
#' Time-only subsetting (\code{[[t]]}, \code{sub_vector()}) is physically
#' partial when every transform advertises physical \code{time_idx} support;
#' QPP pipelines read and decode only intersecting dependency tiles.
#'
#' @slot .state Environment with mutable state: \code{$reader},
#'   \code{$full_mat}, \code{$pipeline_caps}, \code{$last_roi}/\code{$last_mat}.
#' @slot .file Character path to the source LNA file.
#' @param x An \code{LNANeuroVec} object.
#' @param i Numeric index (first spatial dimension or linear index).
#' @param j Numeric index (second spatial dimension).
#' @param k Numeric index (third spatial dimension).
#' @param m Numeric index (fourth, time dimension).
#' @param drop Logical; whether singleton dimensions may be dropped.
#' @param ... Additional arguments passed to methods.
#' @param object An \code{LNANeuroVec} object (for \code{show}).
#' @export
#' @rdname LNANeuroVec-class
setClass("LNANeuroVec",
  slots = c(
    .state = "environment",
    .file  = "character"
  ),
  contains = "SparseNeuroVec"
)

#' Construct an LNANeuroVec
#'
#' Creates a lazy NeuroVec from an LNA archive. Construction reads
#' only metadata (header, mask, transform descriptors). No data is
#' decompressed until an accessor is called.
#'
#' @param file Path to an \code{.lna.h5} file, or an existing
#'   \code{lna_reader} object.
#' @param mask Optional \code{LogicalNeuroVol} or logical 3D array.
#' @param run_id Optional run identifier for multi-run archives.
#' @param label Optional character label.
#' @return An \code{LNANeuroVec} object.
#' @section Lifecycle:
#' Stable. This is the primary lazy `neuroim2` bridge for LNA archives.
#' @export
LNANeuroVec <- function(file, mask = NULL, run_id = NULL, label = "") {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    abort_lna("neuroim2 is required for LNANeuroVec",
              .subclass = "lna_error_missing_package",
              location = "LNANeuroVec")
  }

  # 1. Create or accept reader
  if (inherits(file, "lna_reader")) {
    reader <- file
    file_path <- reader$file
  } else {
    stopifnot(is.character(file), length(file) == 1)
    file_path <- file
    reader <- read_lna(file, lazy = TRUE, run_id = run_id)
  }

  # 2. Read header (cheap — just HDF5 attributes)
  header <- .read_lna_header(reader$h5, run_id = reader$current_run_id)
  dims <- header$dims
  sp <- header$spacing %||% rep(1, 3)
  orig <- header$origin %||% rep(0, 3)

  # 3. Discover transform types for capability check
  tf_group <- reader$h5[["transforms"]]
  transforms <- discover_transforms(tf_group)
  safe_h5_close(tf_group)
  transform_types <- transforms$type
  pipeline_caps <- .pipeline_supports_subset(transform_types)

  # 4. Determine mask and dims with a lightweight peek
  if (is.null(mask) || is.null(dims) || length(dims) < 3) {
    peek <- reader$data(time_idx = 1L)
    if (is.null(dims) || length(dims) < 3) {
      raw <- peek$stash$input
      if (is.array(raw) && length(dim(raw)) >= 3) {
        dims <- dim(raw)[1:3]
      } else if (!is.null(peek$mask_info$mask)) {
        dims <- dim(peek$mask_info$mask)
      } else {
        dims <- c(ncol(as_dense_mat(raw)), 1L, 1L)
      }
    }
    if (is.null(mask) && !is.null(peek$mask_info$mask)) {
      mask <- peek$mask_info$mask
    }
  }
  spatial_dims <- as.integer(dims[1:3])

  # 5. Build mask
  if (is.null(mask)) {
    mask_arr <- array(TRUE, dim = spatial_dims)
  } else if (inherits(mask, "LogicalNeuroVol")) {
    mask_arr <- as.array(mask)
  } else {
    mask_arr <- as.logical(mask)
    dim(mask_arr) <- spatial_dims
  }
  mask_indices <- as.integer(which(mask_arr))
  n_mask <- length(mask_indices)

  # 6. Determine n_time without full read
  n_time <- .lna_determine_ntime(reader, dims, header, n_mask)

  # 7. Build spaces and map. Prefer the complete stored affine so rotation,
  # shear, and axis direction are not lost in the lazy bridge. The
  # spacing/origin path remains for legacy archives without an affine.
  if (!is.null(header$transform)) {
    resolved <- .lna_neurospaces_from_header(
      header,
      inferred_dims = spatial_dims,
      n_time = n_time,
      assume_identity = FALSE,
      location = "LNANeuroVec"
    )
    mask_space <- resolved$space3d
    space4d <- resolved$space4d
  } else {
    mask_space <- neuroim2::NeuroSpace(
      spatial_dims,
      spacing = sp,
      origin = orig
    )
    space4d <- neuroim2::NeuroSpace(
      c(spatial_dims, as.integer(n_time)),
      spacing = sp,
      origin = orig
    )
  }
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, mask_space)
  map <- neuroim2::IndexLookupVol(mask_space, mask_indices)

  # 8. State environment
  state <- new.env(parent = emptyenv())
  state$reader <- reader
  state$pipeline_caps <- pipeline_caps
  state$full_mat <- NULL       # lazy full materialization cache
  state$last_roi <- NULL       # last ROI mask used (for LRU)
  state$last_mat <- NULL       # last subsetted result
  state$full_mat_warned <- FALSE

  # 9. Placeholder data slot
  placeholder <- matrix(0, nrow = as.integer(n_time), ncol = n_mask)

  new("LNANeuroVec",
      data = placeholder,
      space = space4d,
      mask = mask_vol,
      map = map,
      label = as.character(label),
      .state = state,
      .file = file_path)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.lna_get_reader <- function(x) {
  reader <- x@.state$reader
  if (is.null(reader)) stop("LNANeuroVec is closed")
  valid <- tryCatch(reader$h5$is_valid, error = function(e) FALSE)
  if (!isTRUE(valid)) stop("LNANeuroVec reader is closed")
  reader
}

#' Dispatch a subsetted read through the inverse pipeline
#' @keywords internal
.lna_subsetted_read <- function(x, roi_mask = NULL, time_idx = NULL) {
  reader <- .lna_get_reader(x)
  handle <- reader$data(roi_mask = roi_mask, time_idx = time_idx)
  raw <- handle$stash$input
  mat <- as_dense_mat(raw)

  # Normalize orientation to time x voxels
  if (!is.null(roi_mask)) {
    n_roi <- sum(as.logical(roi_mask))
    if (nrow(mat) == n_roi && ncol(mat) != n_roi && ncol(mat) > 1) {
      mat <- t(mat)
    }
  }
  mat
}

#' Full materialization with memory warning (fallback path)
#' @keywords internal
.lna_full_mat <- function(x) {
  if (!is.null(x@.state$full_mat)) return(x@.state$full_mat)

  reader <- .lna_get_reader(x)

  # Warn on large materializations
  n_mask <- length(neuroim2::indices(x))
  n_time <- dim(x)[4]
  est_mb <- n_mask * n_time * 8 / 1e6
  if (est_mb > 100 && !isTRUE(x@.state$full_mat_warned)) {
    warning(sprintf(
      "LNANeuroVec: materializing full matrix (%.0f MB) because pipeline does not fully support ROI subsetting. Consider using [[t]] or sub_vector() for time-only access.",
      est_mb
    ), call. = FALSE)
    x@.state$full_mat_warned <- TRUE
  }

  handle <- reader$data()
  raw <- handle$stash$input
  mat <- as_dense_mat(raw)

  n_spatial <- prod(dim(x)[1:3])

  # Normalize orientation
  if (nrow(mat) == n_mask && ncol(mat) != n_mask && ncol(mat) > 1) {
    mat <- t(mat)
  }
  if (ncol(mat) == n_spatial && n_mask < n_spatial) {
    mat <- mat[, neuroim2::indices(x), drop = FALSE]
  }

  x@.state$full_mat <- mat
  mat
}

#' Read columns for specific sparse voxel indices (smart dispatch)
#' @keywords internal
.lna_read_columns <- function(x, sparse_idx, time_idx = NULL) {
  # If full mat already cached, use it
  if (!is.null(x@.state$full_mat)) {
    rows <- time_idx %||% seq_len(nrow(x@.state$full_mat))
    return(x@.state$full_mat[rows, sparse_idx, drop = FALSE])
  }

  # Check pipeline capability for ROI subsetting
  if (isTRUE(x@.state$pipeline_caps$roi_mask)) {
    # Build ROI mask from sparse indices
    all_indices <- neuroim2::indices(x)
    linear_idx <- all_indices[sparse_idx]
    spatial_dims <- dim(x)[1:3]
    roi_mask <- array(FALSE, dim = spatial_dims)
    roi_mask[linear_idx] <- TRUE

    # Check LRU cache
    roi_key <- paste0(
      paste(sparse_idx, collapse = ","), "|t=",
      if (is.null(time_idx)) "all" else paste(time_idx, collapse = ",")
    )
    if (identical(x@.state$last_roi, roi_key) && !is.null(x@.state$last_mat)) {
      return(x@.state$last_mat)
    }

    mat <- tryCatch(
      .lna_subsetted_read(x, roi_mask = roi_mask, time_idx = time_idx),
      error = function(e) {
        msg <- conditionMessage(e)
        # Closed-reader errors are real lifecycle errors, not a signal to
        # fall back to full materialization.
        if (grepl("closed", msg, ignore.case = TRUE)) {
          stop(e)
        }
        # ROI subsetting failed — fall back to full materialization
        warning(sprintf(
          "LNANeuroVec: ROI-subsetted read failed (%s); falling back to full materialization.",
          msg
        ), call. = FALSE)
        x@.state$pipeline_caps$roi_mask <- FALSE  # don't retry
        full <- .lna_full_mat(x)
        rows <- time_idx %||% seq_len(nrow(full))
        full[rows, sparse_idx, drop = FALSE]
      }
    )

    # Cache
    x@.state$last_roi <- roi_key
    x@.state$last_mat <- mat
    return(mat)
  }

  # Fallback: full materialization (with warning)
  mat <- .lna_full_mat(x)
  rows <- time_idx %||% seq_len(nrow(mat))
  mat[rows, sparse_idx, drop = FALSE]
}

.lna_determine_ntime <- function(reader, dims, header, n_mask) {
  if (length(dims) >= 4) return(as.integer(dims[4]))

  h5 <- reader$h5
  run_id <- reader$current_run_id %||% "run-01"
  for (role in c("quantized", "coefficients", "values")) {
    ds_path <- paste0("scans/", run_id, "/", role)
    if (h5$exists(ds_path)) {
      ds <- h5[[ds_path]]
      ds_dims <- ds$dims
      safe_h5_close(ds)
      if (length(ds_dims) == 2) {
        if (ds_dims[2] == n_mask) return(as.integer(ds_dims[1]))
        if (ds_dims[1] == n_mask) return(as.integer(ds_dims[2]))
        return(as.integer(ds_dims[1]))
      }
      if (length(ds_dims) >= 4) return(as.integer(ds_dims[4]))
    }
  }

  if (h5$exists("transforms")) {
    tf_group <- h5[["transforms"]]
    on.exit(safe_h5_close(tf_group), add = TRUE)
    transforms <- discover_transforms(tf_group)
    qpp_names <- transforms$name[transforms$type == "qpp"]
    for (name in qpp_names) {
      desc <- read_json_descriptor(tf_group, name)
      owner <- desc$run_id
      if ((is.null(owner) || identical(as.character(owner), run_id)) &&
          is.numeric(desc$params$time_count) &&
          length(desc$params$time_count) == 1L &&
          is.finite(desc$params$time_count) && desc$params$time_count >= 1) {
        return(as.integer(desc$params$time_count))
      }
    }
  }

  # Fallback: full peek
  handle <- reader$data()
  mat <- as_dense_mat(handle$stash$input)
  if (nrow(mat) == n_mask) as.integer(ncol(mat)) else as.integer(nrow(mat))
}

# ---------------------------------------------------------------------------
# S4 method overrides
# ---------------------------------------------------------------------------

#' @rdname LNANeuroVec-class
#' @export
setMethod("matricized_access", signature(x = "LNANeuroVec", i = "integer"),
  function(x, i) {
    .lna_read_columns(x, i)
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("matricized_access", signature(x = "LNANeuroVec", i = "numeric"),
  function(x, i) {
    .lna_read_columns(x, as.integer(i))
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("linear_access", signature(x = "LNANeuroVec", i = "numeric"),
  function(x, i) {
    dims <- dim(x)
    spatial_nels <- prod(dims[1:3])
    timepoints <- ((i - 1L) %/% spatial_nels) + 1L
    spatial_offsets <- ((i - 1L) %% spatial_nels) + 1L
    mapped <- neuroim2::lookup(x@map, spatial_offsets)

    output <- numeric(length(i))
    nz <- which(mapped > 0L)
    if (length(nz) == 0) return(output)

    unique_sparse <- sort(unique(mapped[nz]))
    unique_time <- sort(unique(timepoints[nz]))
    mat <- .lna_read_columns(x, unique_sparse, time_idx = unique_time)
    col_remap <- match(mapped[nz], unique_sparse)
    idx_mat <- cbind(match(timepoints[nz], unique_time), col_remap)
    output[nz] <- mat[idx_mat]
    output
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod(
  "[",
  signature(x = "LNANeuroVec", i = "numeric", j = "numeric"),
  function(x, i, j, k, m, ..., drop = TRUE) {
    spatial_dims <- dim(x)[1:3]
    if (missing(k)) {
      k <- seq_len(spatial_dims[3])
    }
    if (missing(m)) {
      m <- seq_len(dim(x)[4])
    }

    vmat <- as.matrix(expand.grid(i, j, k))
    ind <- vmat[, 1] +
      (vmat[, 2] - 1L) * spatial_dims[1] +
      (vmat[, 3] - 1L) * spatial_dims[1] * spatial_dims[2]

    mapped <- neuroim2::lookup(x, ind)
    keep <- mapped > 0
    dimout <- c(length(i), length(j), length(k), length(m))

    if (sum(keep) == 0) {
      out <- array(0, dimout)
      return(if (drop) base::drop(out) else out)
    }

    out <- array(0, dimout)
    nz_idx <- which(keep)
    unique_sparse <- sort(unique(mapped[nz_idx]))
    vals <- .lna_read_columns(x, unique_sparse, time_idx = m)
    col_remap <- match(mapped[nz_idx], unique_sparse)
    coords <- vmat[nz_idx, , drop = FALSE]

    for (col in seq_along(nz_idx)) {
      ii <- match(coords[col, 1], i)
      jj <- match(coords[col, 2], j)
      kk <- match(coords[col, 3], k)
      out[ii, jj, kk, ] <- vals[, col_remap[col]]
    }

    if (drop) {
      base::drop(out)
    } else {
      out
    }
  }
)

# ---------------------------------------------------------------------------
# Direct series overrides (parent accesses @data directly, not matricized_access)
# ---------------------------------------------------------------------------

#' @rdname LNANeuroVec-class
#' @export
setMethod("series", signature(x = "LNANeuroVec", i = "integer"),
  function(x, i, j, k, ...) {
    if (!missing(j) && !missing(k)) {
      # series(x, i, j, k) — 3D coordinate form
      spatial_dims <- dim(x)[1:3]
      lin <- i + (j - 1L) * spatial_dims[1] + (k - 1L) * spatial_dims[1] * spatial_dims[2]
      mapped <- neuroim2::lookup(x@map, lin)
      if (all(mapped == 0)) return(rep(0, dim(x)[4]))
      nz <- which(mapped > 0)
      mat <- .lna_read_columns(x, mapped[nz])
      if (length(nz) == 1) return(as.numeric(mat[, 1]))
      return(mat)
    }
    # series(x, i) — linear voxel index form
    mapped <- neuroim2::lookup(x@map, i)
    nz <- which(mapped > 0)
    if (length(nz) == 0) {
      return(matrix(0, nrow = dim(x)[4], ncol = length(i)))
    }
    out <- matrix(0, nrow = dim(x)[4], ncol = length(i))
    mat <- .lna_read_columns(x, mapped[nz])
    out[, nz] <- mat
    if (ncol(out) == 1) return(as.numeric(out[, 1]))
    out
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("series", signature(x = "LNANeuroVec", i = "numeric"),
  function(x, i, j, k, ...) {
    if (!missing(j) && !missing(k)) {
      return(series(x, as.integer(i), as.integer(j), as.integer(k), ...))
    }
    series(x, as.integer(i), ...)
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("series", signature(x = "LNANeuroVec", i = "matrix"),
  function(x, i, ...) {
    # i is an N x 3 matrix of coordinates
    spatial_dims <- dim(x)[1:3]
    lin <- i[, 1] + (i[, 2] - 1L) * spatial_dims[1] +
           (i[, 3] - 1L) * spatial_dims[1] * spatial_dims[2]
    series(x, as.integer(lin), ...)
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("series", signature(x = "LNANeuroVec", i = "LogicalNeuroVol"),
  function(x, i, ...) {
    idx <- which(as.logical(i))
    series(x, as.integer(idx), ...)
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("[[", signature(x = "LNANeuroVec", i = "numeric"),
  function(x, i) {
    stopifnot(length(i) == 1L)
    # Always use true time-subsetted read
    mat <- .lna_subsetted_read(x, time_idx = as.integer(i))
    dat <- if (nrow(mat) == 1) as.numeric(mat[1, ]) else as.numeric(mat)
    xs <- neuroim2::space(x)
    newdim <- dim(xs)[1:3]
    bspace <- neuroim2::NeuroSpace(newdim,
      spacing = neuroim2::spacing(xs),
      origin = neuroim2::origin(xs))
    neuroim2::SparseNeuroVol(dat, bspace, indices = neuroim2::indices(x))
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("sub_vector", signature(x = "LNANeuroVec", i = "numeric"),
  function(x, i) {
    # Always use true time-subsetted read
    mat <- .lna_subsetted_read(x, time_idx = as.integer(i))
    xs <- neuroim2::space(x)
    newdim <- c(dim(x)[1:3], length(i))
    bspace <- neuroim2::NeuroSpace(newdim,
      spacing = neuroim2::spacing(xs),
      origin = neuroim2::origin(xs))
    neuroim2::SparseNeuroVec(mat, bspace, mask = x@mask)
  })

#' @rdname LNANeuroVec-class
#' @export
setMethod("show", signature(object = "LNANeuroVec"),
  function(object) {
    dims <- dim(object)
    n_mask <- length(neuroim2::indices(object))
    has_reader <- !is.null(object@.state$reader)
    has_full <- !is.null(object@.state$full_mat)
    caps <- object@.state$pipeline_caps
    roi_ok <- isTRUE(caps$roi_mask)
    strategy <- caps$roi_strategy %||% if (roi_ok) "direct" else "fallback"
    strategy_label <- switch(
      strategy,
      "coefficient" = "subset (coefficient)",
      "shared-coefficient" = "subset (shared coefficient)",
      "parcel-gather" = "subset (QPP parcel gather)",
      "dense-pass-through" = "subset (dense pass-through)",
      "mixed" = "subset (mixed)",
      if (roi_ok) "subset" else "full-mat fallback"
    )

    status <- if (has_full) "materialized" else if (has_reader) "lazy (open)" else "closed"

    cat("LNANeuroVec\n")
    cat("  File:", object@.file, "\n")
    cat("  Dimensions:", paste(dims[1:3], collapse = " x "), "\n")
    cat("  Time Points:", dims[4], "\n")
    cat("  Masked Voxels:", n_mask, "/", prod(dims[1:3]), "\n")
    cat("  Status:", status, "\n")
    cat("  ROI dispatch:", strategy_label, "\n")
    cat("  Physical time I/O:", if (isTRUE(caps$physical_time_idx)) "partial" else "full", "\n")
    cat("  Physical ROI I/O:", if (isTRUE(caps$physical_roi_mask)) "partial" else "full", "\n")
  })

# ---------------------------------------------------------------------------
# Lifecycle
# ---------------------------------------------------------------------------

#' Close an LNANeuroVec
#'
#' Releases the HDF5 handle and clears caches.
#' @param x An \code{LNANeuroVec} object.
#' @return Invisibly, the object.
#' @export
lna_close <- function(x) {
  UseMethod("lna_close")
}

#' @export
lna_close.LNANeuroVec <- function(x) {
  if (!is.null(x@.state$reader)) {
    x@.state$reader$close()
    x@.state$reader <- NULL
  }
  x@.state$last_roi <- NULL
  x@.state$last_mat <- NULL
  invisible(x)
}

#' Check if data has been fully materialized
#' @param x An \code{LNANeuroVec} object.
#' @return Logical.
#' @export
lna_loaded <- function(x) {
  !is.null(x@.state$full_mat)
}

# ---------------------------------------------------------------------------
# Convenience generic
# ---------------------------------------------------------------------------

#' Coerce to LNANeuroVec
#' @param x A file path or \code{lna_reader}.
#' @param ... Additional arguments passed to \code{LNANeuroVec()}.
#' @return An \code{LNANeuroVec} object.
#' @export
as_lna_neurovec <- function(x, ...) {
  UseMethod("as_lna_neurovec")
}

#' @export
as_lna_neurovec.character <- function(x, mask = NULL, ...) {
  LNANeuroVec(file = x, mask = mask, ...)
}

#' @export
as_lna_neurovec.lna_reader <- function(x, mask = NULL, ...) {
  LNANeuroVec(file = x, mask = mask, ...)
}
