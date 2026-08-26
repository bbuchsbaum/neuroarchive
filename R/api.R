#' Write data to an LNA file
#'
#' Compresses one or more fMRI runs using a sequence of transforms and
#' stores the result in an `.lna.h5` file. Parameter values for each
#' transform are resolved by merging the JSON schema defaults, package
#' options set via `lna_options()`, and any user supplied
#' `transform_params` (later values override earlier ones).
#'
#' @param x Numeric array or `DenseNeuroVec` (or list of those). Each
#'   array must have at least three dimensions (`x`, `y`, `z`, and
#'   optionally `time`). 3D inputs (`DenseNeuroVol` or 3D array) are
#'   expanded to 4D. Lists denote multiple runs. The transform pipeline is
#'   restarted and fitted independently for each list element unless the first
#'   transform is the explicit cross-run `aggregate_runs` operation.
#' @param file Path to the output `.h5` file. If `NULL`, writing occurs
#'   in memory using the HDF5 core driver and no file is created. The
#'   returned result then contains `file = NULL`.
#' @param transforms Character vector naming the transforms to apply in
#'   forward order (e.g., `c("quant", "basis")`).
#' @param transform_params Named list of parameters for the specified
#'   transforms.
#' @param mask Optional: a `LogicalNeuroVol` or 3D logical array used to
#'   subset voxels prior to compression.
#' @param header Optional named list of header attributes to store under
#'   `/header`. When `NULL` and `x` inherits from `NeuroObj` the header is
#'   created from its `NeuroSpace`.
#' @param run_geometry Optional run-scoped spatial descriptors. For multiple
#'   runs, this must be a named list keyed by the resolved `run_id` values.
#'   Each descriptor requires `dims` and a full `transform` (or `affine`).
#' @param geometry_policy `"optional"` preserves the generic low-level API;
#'   `"fmri"` requires geometry for every run and never fabricates identity
#'   geometry.
#' @param plugins Optional named list saved under the `/plugins` group.
#' @param block_table Optional data frame specifying spatial block
#'   coordinates stored at `/spatial/block_table`. Columns must contain
#'   1-based voxel indices in masked space when a mask is provided.
#' @param run_id Optional character vector of run identifiers. When `x`
#'   is a list these override `names(x)`; otherwise a single identifier
#'   is used for the lone run. Values must be unique and the vector length must
#'   match the number of runs.
#' @param checksum Character string specifying checksum mode. One of
#'   `"sha256"` (default) or `"none"`. When `"sha256"` a checksum of the
#'   final file is computed and stored in the `/lna_checksum` attribute.
#'   This requires closing and reopening the file once writing has
#'   finished.
#' @param ... Reserved for method-specific extensions.
#' @return Invisibly returns a list with elements `file`, `plan`, `header`,
#'   `input_bytes`, `output_bytes`, `ratio`, `transforms_applied`, and
#'   `checksum` with class `"lna_write_result"`.
#' @section Lifecycle:
#' Stable. This is the primary LNA archive writer and part of the core
#' compatibility contract.
#' @details For parallel workflows create a unique temporary file and
#'   rename it into place once writing succeeds. The underlying HDF5 file
#'   is opened with mode `"w"`, truncating any existing file at `file`.
#'   Run-local transform descriptors and payloads are stored separately below
#'   `/scans/<run_id>`. A shared mask, when supplied, must match every run.
#'   Explicit geometry is stored under `/header/runs/<run_id>`; a compatible
#'   global spatial header is retained only for a single or identical grid.
#' @seealso read_lna, validate_lna
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
#' result <- write_lna(arr, tmp, transforms = "quant")
#' read_lna(tmp)
#' @importFrom hdf5r H5File
#' @export
write_lna <- function(x, file = NULL, transforms = character(),
                      transform_params = list(), mask = NULL,
                      header = NULL, plugins = NULL, block_table = NULL,
                      run_id = NULL, checksum = c("sha256", "none"),
                      run_geometry = NULL,
                      geometry_policy = c("optional", "fmri"), ...) {
  UseMethod("write_lna")
}

#' @export
write_lna.default <- function(x, file = NULL, transforms = character(),
                              transform_params = list(), mask = NULL,
                              header = NULL, plugins = NULL, block_table = NULL,
                              run_id = NULL, checksum = c("sha256", "none"),
                              run_geometry = NULL,
                              geometry_policy = c("optional", "fmri"), ...) {

  checksum <- match.arg(checksum)
  geometry_policy <- match.arg(geometry_policy)

  # Compute input size for compression summary
  input_bytes <- as.numeric(object.size(x))

  info <- open_output_h5(file)
  h5 <- info$h5
  archive_complete <- FALSE
  on.exit(
    close_output_h5(info, cleanup_temp = !archive_complete),
    add = TRUE
  )

  is_multi_input <- is.list(x) && !methods::is(x, "NeuroObj") && length(x) > 1L
  if (is.null(header) && is.null(run_geometry) && !is_multi_input) {
    header <- derive_header_from_input(x)
  }

  transform_params <- .normalise_transform_pipeline(
    transforms,
    transform_params,
    warn_semantics = TRUE
  )

  result <- core_write(
    x = x,
    transforms = transforms,
    transform_params = transform_params,
    mask = mask,
    header = header,
    run_geometry = run_geometry,
    geometry_policy = geometry_policy,
    plugins = plugins,
    run_id = run_id,
    h5 = h5
  )
 
  validate_block_table(block_table, result$handle$mask_info$active_voxels)
  
  plugins_from_handle <- result$handle$meta$plugins
  if (length(plugins_from_handle) == 0) plugins_from_handle <- NULL
  header_from_handle <- result$global_header %||%
    result$handle$meta$header %||% list()
  run_geometry_from_handle <- result$run_geometry %||% list()

  write_block_table_dataset(h5, block_table)

  materialise_plan(
    h5,
    result$plan,
    checksum = checksum,
    header = header_from_handle,
    run_geometry = run_geometry_from_handle,
    geometry_policy = geometry_policy,
    plugins = plugins_from_handle,
    run_ids = result$handle$run_ids
  )

  # From this point onward the temporary archive is complete. If publication
  # fails, preserve it for recovery rather than treating it as a partial write.
  archive_complete <- TRUE

  # Close the HDF5 file before atomic rename (but don't cleanup temp file yet)
  close_output_h5(info, cleanup_temp = FALSE)

  # Atomic rename: move temp file to target path
  info <- finalise_output_h5(info)

  # Compute output size and compression ratio
  output_file <- info$file
  output_bytes <- if (!is.null(output_file) && file.exists(output_file)) {
    file.size(output_file)
  } else {
    NA_real_
  }

  ratio <- if (!is.na(output_bytes) && output_bytes > 0) {
    input_bytes / output_bytes
  } else {
    NA_real_
  }

  handle_out <- result$handle
  if (!is.null(handle_out) && "h5" %in% names(handle_out)) {
    handle_out$h5 <- NULL
  }

  lnaobj <- list(
    file = if (info$in_memory) NULL else output_file,
    plan = result$plan,
    handle = handle_out,
    transform_params = header_from_handle$transform_params,
    header = header_from_handle,
    run_geometry = run_geometry_from_handle,
    geometry_policy = geometry_policy,
    input_bytes = input_bytes,
    output_bytes = output_bytes,
    ratio = ratio,
    transforms_applied = transforms,
    checksum = checksum
  )
  class(lnaobj) <- c("lna_write_result", "list")
  invisible(lnaobj)
}

#' Read data from an LNA file
#'
#' Loads data from an `.lna.h5` file using `core_read`.  When
#' `lazy = TRUE` the function returns an `lna_reader` object that keeps
#' the HDF5 handle open for on-demand reconstruction of the data. If the
#' selection matches multiple runs, both eager and lazy reads return a named
#' list with one independent result per run.
#'
#' @param file Path to an LNA file on disk.
#' @param run_id Character vector of run identifiers or glob patterns. Passed to
#'   `core_read` for selection of specific runs.
#' @param allow_plugins Character string specifying how to handle
#'   transforms that require optional packages. One of
#'   \code{"installed"} (default), \code{"none"}, or \code{"prompt"}.
#'   Non-interactive sessions treat \code{"prompt"} the same as
#'   \code{"installed"}.  When a required transform implementation is
#'   missing, \code{"installed"} emits a warning and skips that
#'   transform. Interactive use of \code{"prompt"} will ask whether to
#'   continue; declining aborts reading.
#' @param validate Validation mode. One of `"quick"` (default), `TRUE`, or
#'   `FALSE`. `"quick"` performs bounded structural, required-version, QPP
#'   relationship, extent, and index checks. `TRUE` first runs full
#'   [validate_lna()] validation, including checksums and every QPP tile.
#'   `FALSE` skips that preflight. For current QPP archives, the decoder still
#'   verifies selected tile bytes and parcel steps before returning or caching
#'   values in every mode.
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param roi_mask Optional ROI mask used to subset voxels before
#'   applying transforms.
#' @param time_idx Optional vector of time indices for subsetting
#'   volumes prior to transformation.
#' @param lazy Logical. If `TRUE`, the HDF5 file remains open and the
#'   returned `lna_reader` can load data lazily.
#' @return For one selected run, an `lna_reader` when `lazy = TRUE` or a
#'   `DataHandle` otherwise. For multiple selected runs, a named list of
#'   independent `lna_reader` or `DataHandle` objects, respectively.
#' @section Lifecycle:
#' Stable. This is the primary LNA archive reader and part of the core
#' compatibility contract.
#' @seealso write_lna, validate_lna
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(16), dim = c(4, 4, 1, 1))
#' write_lna(arr, tmp, transforms = "quant")
#' read_lna(tmp)
#' @export
read_lna <- function(file, run_id = NULL,
                     allow_plugins = c("installed", "none", "prompt"),
                     validate = "quick",
                     output_dtype = c("float32", "float64", "float16"),
                     roi_mask = NULL, time_idx = NULL,
                     lazy = FALSE) {
  if (!(is.character(file) && length(file) == 1)) {
    abort_lna(
      "file must be a path",
      .subclass = "lna_error_validation",
      location = "read_lna:file"
    )
  }
  output_dtype <- match.arg(output_dtype)
  allow_plugins <- match.arg(allow_plugins)

  # Handle validation parameter - can be "quick", TRUE, FALSE, or logical
  if (is.character(validate)) {
    validate <- match.arg(validate, choices = "quick")
  } else if (!is.logical(validate) || length(validate) != 1) {
    abort_lna(
      "validate must be one of 'quick', TRUE, or FALSE",
      .subclass = "lna_error_validation",
      location = "read_lna:validate"
    )
  }

  # Perform quick validation if requested
  if (identical(validate, "quick")) {
    quick_validate_lna(file)
    # Convert to FALSE for core_read since we've already validated
    validate_for_core <- FALSE
  } else {
    validate_for_core <- validate
  }

  args <- list(
    file = file,
    run_id = run_id,
    allow_plugins = allow_plugins,
    validate = validate_for_core,
    output_dtype = output_dtype,
    roi_mask = roi_mask,
    time_idx = time_idx
  )

  if (lazy) {
    h5 <- open_h5(file, mode = "r")
    on.exit(close_h5_safely(h5), add = TRUE)
    runs <- resolve_run_ids(run_id, discover_run_ids(h5))
    if (length(runs) == 0L) {
      abort_lna(
        "run_id did not match any runs",
        .subclass = "lna_error_run_id",
        location = "read_lna:run_id"
      )
    }
    close_h5_safely(h5)
    h5 <- NULL

    readers <- list()
    on.exit({
      lapply(readers, function(reader) try(reader$close(), silent = TRUE))
    }, add = TRUE)
    for (rid in runs) {
      run_args <- args
      run_args$run_id <- rid
      readers[[rid]] <- lna_reader$new(
        file = file,
        core_read_args = run_args
      )
    }
    on.exit(NULL, add = FALSE)
    if (length(readers) == 1L) readers[[1L]] else readers
  } else {
    args$lazy <- FALSE
    do.call(core_read, args)
  }
}

#' Convenience alias for `write_lna`
#'
#' `compress_fmri()` simply forwards its arguments to `write_lna()` without
#' altering the dimensions of the input.
#'
#' @param ... Arguments forwarded to `write_lna()`.
#' @seealso write_lna
#' @export
compress_fmri <- function(...) write_lna(...)

#' Convenience alias for `read_lna`
#'
#' `open_lna()` simply forwards its arguments to `read_lna()`.
#'
#' @inheritParams read_lna
#' @seealso read_lna
#' @export
open_lna <- read_lna

# -------------------------------------------------------------------------
# Internal helper functions

open_output_h5 <- function(path) {
  if (is.null(path)) {
    tmp <- tempfile(fileext = ".h5")
    # For in-memory operation, we'll create a regular file in temp directory
    # Modern hdf5r doesn't support the core driver parameters in the same way
    h5 <- open_h5(tmp, mode = "w")
    if (is.null(h5) || !h5$is_valid) {
      stop("Failed to create temporary HDF5 file")
    }
    # Keep explicit cleanup flag separate from in_memory semantics.
    list(h5 = h5, file = tmp, target_path = NULL, in_memory = TRUE, cleanup_file = TRUE)
  } else {
    # Write to temp file in same directory for atomic rename
    target_dir <- dirname(path)
    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    tmp <- tempfile(tmpdir = target_dir, fileext = ".lna.tmp")
    h5 <- open_h5(tmp, mode = "w")
    if (is.null(h5) || !h5$is_valid) {
      unlink(tmp, force = TRUE)
      stop(sprintf("Failed to open HDF5 file for writing (target: '%s')", path))
    }
    list(h5 = h5, file = tmp, target_path = path, in_memory = FALSE, cleanup_file = FALSE)
  }
}

close_output_h5 <- function(info, cleanup_temp = TRUE) {
  if (!is.null(info$h5) && info$h5$is_valid) {
    # Best-effort flush/close to avoid transient read-after-write lock races
    tryCatch(info$h5$flush(), error = function(e) invisible(NULL))
    tryCatch(info$h5$close_all(), error = function(e) {
      warning(paste("Error closing output HDF5 handle:", conditionMessage(e)))
    })
  }
  if (isTRUE(info$cleanup_file) && file.exists(info$file)) {
    unlink(info$file, force = TRUE)
  }
  # Clean up temp file if atomic rename hasn't happened yet (error path)
  # Only cleanup if requested (to allow explicit close before finalise)
  if (cleanup_temp && !is.null(info$target_path) && file.exists(info$file)) {
    unlink(info$file, force = TRUE)
  }
  invisible(NULL)
}

finalise_output_h5 <- function(info) {
  if (!is.null(info$target_path)) {
    ok <- file.rename(info$file, info$target_path)
    if (!isTRUE(ok)) {
      # file.rename can fail across filesystems; fall back to copy+delete
      copy_ok <- file.copy(info$file, info$target_path, overwrite = TRUE)
      if (!isTRUE(copy_ok)) {
        abort_lna(
          paste0(
            "Failed to publish completed LNA archive; temporary archive retained at '",
            info$file,
            "'"
          ),
          temporary_file = info$file,
          target_file = info$target_path,
          .subclass = "lna_error_write",
          location = "api::finalise_output_h5"
        )
      }
      unlink(info$file, force = TRUE)
    }
    # Update file path in info so downstream code (like output_bytes) uses the right path
    info$file <- info$target_path
    info$target_path <- NULL
  }
  info
}

derive_header_from_input <- function(x) {
  src <- if (is.list(x)) x[[1]] else x
  if (methods::is(src, "NeuroObj")) {
    # Try global space function first (for testing), then neuroim2 namespace
    spc <- tryCatch({
      if (exists("space", envir = .GlobalEnv, mode = "function")) {
        get("space", envir = .GlobalEnv)(src)
      } else {
        space(src)
      }
    }, error = function(e) NULL)
    if (!is.null(spc)) {
      return(neuroim2_space_to_lna_header(spc))
    }
  }
  NULL
}

validate_block_table <- function(block_table, max_idx) {
  if (is.null(block_table)) return(invisible(NULL))
  if (!is.data.frame(block_table)) {
    abort_lna(
      "block_table must be a data frame",
      .subclass = "lna_error_validation",
      location = "write_lna:block_table"
    )
  }
  if (nrow(block_table) > 0) {
    num_cols <- vapply(block_table, is.numeric, logical(1))
    if (!all(num_cols)) {
      abort_lna(
        "block_table columns must be numeric",
        .subclass = "lna_error_validation",
        location = "write_lna:block_table"
      )
    }
    coords <- unlist(block_table)
    if (any(is.na(coords)) || any(coords < 1, na.rm = TRUE)) {
      abort_lna(
        "block_table coordinates must be non-missing and >= 1",
        .subclass = "lna_error_validation",
        location = "write_lna:block_table"
      )
    }
    if (!is.null(max_idx) && any(coords > max_idx, na.rm = TRUE)) {
      abort_lna(
        "block_table coordinates exceed masked voxel count",
        .subclass = "lna_error_validation",
        location = "write_lna:block_table"
      )
    }
  }
  invisible(NULL)
}

write_block_table_dataset <- function(h5, block_table) {
  if (is.null(block_table)) return(invisible(NULL))
  bt_matrix <- as.matrix(block_table)
  h5_write_dataset(h5[["/"]], "spatial/block_table", bt_matrix)
  invisible(NULL)
}

#' Print method for lna_write_result
#'
#' @param x An object of class lna_write_result
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @method print lna_write_result
#' @export
print.lna_write_result <- function(x, ...) {
  fmt_bytes <- function(b) {
    if (is.na(b)) return("NA")
    if (b >= 1e9) sprintf("%.1f GB", b / 1e9)
    else if (b >= 1e6) sprintf("%.1f MB", b / 1e6)
    else if (b >= 1e3) sprintf("%.1f KB", b / 1e3)
    else sprintf("%d B", as.integer(b))
  }

  transforms_str <- if (length(x$transforms_applied) > 0) {
    paste0("[", paste(x$transforms_applied, collapse = "+"), "]")
  } else {
    "[none]"
  }

  checksum_str <- if (!is.null(x$checksum) && nchar(x$checksum) > 0 && x$checksum != "none") {
    "sha256:ok"
  } else {
    "no checksum"
  }

  ratio_str <- if (!is.na(x$ratio)) {
    sprintf("(%.1fx)", x$ratio)
  } else {
    "(NA)"
  }

  cat(sprintf("LNA written: %s -> %s %s %s %s\n",
    fmt_bytes(x$input_bytes),
    fmt_bytes(x$output_bytes),
    ratio_str,
    transforms_str,
    checksum_str))
  invisible(x)
}
