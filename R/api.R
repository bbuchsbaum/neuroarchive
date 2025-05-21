#' Write data to an LNA file
#'
#' Compresses and stores a neuroimaging object using the specified
#' transform pipeline.  Parameter values are resolved by merging
#' transform schema defaults, package options set via
#' `lna_options()`, and the user supplied `transform_params`
#' (later values override earlier ones).
#'
#' @param x Input object passed to `core_write`.
#' @param file Path to output `.h5` file. If `NULL`, writing is performed
#'   in memory using the HDF5 core driver and no file is created. The
#'   returned list then has `file = NULL`.
#' @param transforms Character vector of transform types.
#' @param transform_params Named list of transform parameters.
#' @param mask Optional mask passed to `core_write`.
#' @param header Optional named list of header attributes.
#' @param plugins Optional named list saved under `/plugins/`.
#' @param block_table Optional data frame specifying spatial block coordinates
#'   stored at `/spatial/block_table`. Coordinate columns must contain
#'   1-based voxel indices in masked space when a mask is provided.
#' @return Invisibly returns a list with elements `file`, `plan`, and
#'   `header` with class `"lna_write_result"`.
#' @details For parallel workflows use a unique temporary file and
#'   `file.rename()` it to the final path once writing succeeds.
#'   The underlying HDF5 handle is opened with mode `"w"` which
#'   truncates any existing file at `file`.
#' @seealso read_lna, validate_lna
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(16), dim = c(4, 4, 1, 1))
#' write_lna(arr, tmp, transforms = "quant")
#' @export
write_lna <- function(x, file = NULL, transforms = character(),
                      transform_params = list(), mask = NULL,
                      header = NULL, plugins = NULL, block_table = NULL) {

  in_memory <- FALSE
  if (is.null(file)) {
    tmp <- tempfile(fileext = ".h5")
    file <- tmp
    in_memory <- TRUE
  }

  result <- core_write(x = x, transforms = transforms,
                       transform_params = transform_params,
                       mask = mask, header = header, plugins = plugins)

  if (!is.null(block_table)) {
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
      max_idx <- result$handle$mask_info$active_voxels
      if (!is.null(max_idx) && any(coords > max_idx, na.rm = TRUE)) {
        abort_lna(
          "block_table coordinates exceed masked voxel count",
          .subclass = "lna_error_validation",
          location = "write_lna:block_table"
        )
      }
    }
  }

  if (in_memory) {
    h5 <- open_h5(file, mode = "w", driver = "core", backing_store = FALSE)
  } else {
    h5 <- open_h5(file, mode = "w")
  }

  on.exit(close_h5_safely(h5))

  materialise_plan(h5, result$plan,
                   header = result$handle$meta$header,
                   plugins = result$handle$meta$plugins)


  if (!is.null(block_table)) {
    h5_write_dataset(h5[["/"]], "spatial/block_table", as.matrix(block_table))
  }

  close_h5_safely(h5)

  out_file <- if (in_memory) NULL else file
  out <- list(file = out_file, plan = result$plan,
              header = result$handle$meta$header)
  class(out) <- c("lna_write_result", class(out))
  invisible(out)
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
#' @param lazy Logical. If `TRUE`, the HDF5 file remains open and the
#'   returned `lna_reader` can load data lazily.
#' @return The result of `core_read`: a `DataHandle` for a single run or a list
#'   of `DataHandle` objects when multiple runs are loaded.
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

  if (lazy) {
    lna_reader$new(
      file = file,
      core_read_args = list(
        run_id = run_id,
        allow_plugins = allow_plugins,
        validate = validate,
        output_dtype = output_dtype,
        roi_mask = roi_mask,
        time_idx = time_idx
      )
    )
  } else {
    core_read(
      file = file,
      run_id = run_id,
      allow_plugins = allow_plugins,
      validate = validate,
      output_dtype = output_dtype,
      roi_mask = roi_mask,
      time_idx = time_idx,
      lazy = FALSE
    )
  }
}

#' Convenience alias for `write_lna`
#'
#' `compress_fmri()` simply forwards its arguments to `write_lna()`.
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
