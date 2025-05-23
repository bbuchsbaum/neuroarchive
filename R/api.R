#' @importFrom hdf5r H5File
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
#'   expanded to 4D. Lists denote multiple runs.
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
#' @param plugins Optional named list saved under the `/plugins` group.
#' @param block_table Optional data frame specifying spatial block
#'   coordinates stored at `/spatial/block_table`. Columns must contain
#'   1-based voxel indices in masked space when a mask is provided.
#' @param run_id Optional character vector of run identifiers. When `x`
#'   is a list these override `names(x)`; otherwise a single identifier
#'   is used for the lone run.
#' @param checksum Character string specifying checksum mode. One of
#'   `"none"` (default) or `"sha256"`. When `"sha256"` a checksum of the
#'   final file is computed and stored in the `/lna_checksum` attribute.
#'   This requires closing and reopening the file once writing has
#'   finished.
#' @return Invisibly returns a list with elements `file`, `plan`, and
#'   `header` and class `"lna_write_result"`.
#' @details For parallel workflows create a unique temporary file and
#'   rename it into place once writing succeeds. The underlying HDF5 file
#'   is opened with mode `"w"`, truncating any existing file at `file`.
#' @seealso read_lna, validate_lna
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
#' write_lna(arr, tmp, transforms = "quant")
#' read_lna(tmp)
#' @export
write_lna <- function(x, ...) {
  UseMethod("write_lna")
}

#' @export
write_lna.default <- function(x, file = NULL, transforms = character(),
                      transform_params = list(), mask = NULL,
                      header = NULL, plugins = NULL, block_table = NULL,
                      run_id = NULL, checksum = c("none", "sha256")) {

  checksum <- match.arg(checksum)

  info <- open_output_h5(file)
  h5 <- info$h5
  on.exit(close_output_h5(info), add = TRUE)

  if (is.null(header)) {
    header <- derive_header_from_input(x)
  }

  result <- core_write(
    x = x,
    transforms = transforms,
    transform_params = transform_params,
    mask = mask,
    header = header,
    plugins = plugins,
    run_id = run_id
  )
 
  validate_block_table(block_table, result$handle$mask_info$active_voxels)
  
  plugins_from_handle <- result$handle$meta$plugins
  if (length(plugins_from_handle) == 0) plugins_from_handle <- NULL
  header_from_handle <- result$handle$meta$header %||% list()

  materialise_plan(
    h5,
    result$plan,
    checksum = checksum,
    header = header_from_handle,
    plugins = plugins_from_handle
  )

  write_block_table_dataset(h5, block_table)

  lnaobj <- list(
    file = if (info$in_memory) NULL else info$file,
    plan = result$plan,
    transform_params = header_from_handle$transform_params,
    header = header_from_handle
  )
  class(lnaobj) <- c("lna_write_result", "list")
  lnaobj
}

#' Read data from an LNA file
#'
#' Loads data from an `.lna.h5` file using `core_read`.  When
#' `lazy = TRUE` the function returns an `lna_reader` object that keeps
#' the HDF5 handle open for on-demand reconstruction of the data.
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
#' @param validate Logical flag for validation; forwarded to `core_read`.
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param roi_mask Optional ROI mask used to subset voxels before
#'   applying transforms.
#' @param time_idx Optional vector of time indices for subsetting
#'   volumes prior to transformation.
#' @param lazy Logical. If `TRUE`, the HDF5 file remains open and the
#'   returned `lna_reader` can load data lazily.
#' @return When `lazy = TRUE`, an `lna_reader` object.  Otherwise the result of
#'   `core_read`: a `DataHandle` for a single run or a list of `DataHandle`
#'   objects when multiple runs are loaded.
#' @seealso write_lna, validate_lna
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(16), dim = c(4, 4, 1, 1))
#' write_lna(arr, tmp, transforms = "quant")
#' read_lna(tmp)
#' @export
read_lna <- function(file, run_id = NULL,
                     allow_plugins = c("installed", "none", "prompt"),
                     validate = FALSE,
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

  args <- list(
    file = file,
    run_id = run_id,
    allow_plugins = allow_plugins,
    validate = validate,
    output_dtype = output_dtype
  )

  if (!is.null(roi_mask)) args$roi_mask <- roi_mask
  if (!is.null(time_idx)) args$time_idx <- time_idx

  if (lazy) {
    lna_reader$new(
      file = file,
      core_read_args = args
    )
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
#' @inheritParams write_lna
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
    h5 <- hdf5r::H5File$new(
      tmp,
      mode = "w",
      driver = "core",
      driver_info = list(backing_store = FALSE)
    )
    if (is.null(h5) || !h5$is_valid) {
      stop("Failed to create in-memory HDF5 file with H5File$new")
    }
    warning(
      sprintf(
        "In-memory HDF5 file (core driver) created via H5File$new using temp name: %s",
        tmp
      ),
      call. = FALSE
    )
    list(h5 = h5, file = tmp, in_memory = TRUE)
  } else {
    h5 <- open_h5(path, mode = "w")
    if (is.null(h5) || !h5$is_valid) {
      stop(sprintf("Failed to open HDF5 file '%s'", path))
    }
    list(h5 = h5, file = path, in_memory = FALSE)
  }
}

close_output_h5 <- function(info) {
  if (!is.null(info$h5) && info$h5$is_valid) {
    info$h5$close_all()
  }
  if (isTRUE(info$in_memory) && file.exists(info$file)) {
    unlink(info$file, force = TRUE)
  }
  invisible(NULL)
}

derive_header_from_input <- function(x) {
  src <- if (is.list(x)) x[[1]] else x
  if (methods::is(src, "NeuroObj")) {
    spc <- tryCatch(space(src), error = function(e) NULL)
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
