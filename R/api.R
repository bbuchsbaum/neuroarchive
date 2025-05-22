#' Write data to an LNA file
#'
#' Compresses one or more fMRI runs using a sequence of transforms and
#' stores the result in an `.lna.h5` file. Parameter values for each
#' transform are resolved by merging the JSON schema defaults, package
#' options set via `lna_options()`, and any user supplied
#' `transform_params` (later values override earlier ones).
#'
#' @param x Numeric array or list of arrays. Each array must have at
#'   least three dimensions (`x`, `y`, `z`, and optionally `time`). If a
#'   list is supplied each element represents a run. When `x` is a single
#'   array it is treated as one run.
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
#'   `/header`.
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
  cat("\n[write_lna] Entry\n")
  cat(paste0("[write_lna] Input file arg: ", ifelse(is.null(file), "NULL", file), "\n"))
  checksum <- match.arg(checksum)

  in_memory <- FALSE
  file_to_use <- file # Will be updated if file is NULL

  if (is.null(file)) {
    cat("[write_lna] file is NULL, preparing for in-memory HDF5.\n")
    # Ensure tmp is uniquely named to avoid clashes if this function is called multiple times with file=NULL
    tmp_for_mem <- tempfile(fileext = ".h5") 
    cat(paste0("[write_lna] Temporary file for in-memory mode: ", tmp_for_mem, "\n"))
    file_to_use <- tmp_for_mem
    in_memory <- TRUE
  } else {
    cat(paste0("[write_lna] file is provided: ", file, ", preparing for disk-based HDF5.\n"))
    file_to_use <- file
  }

  cat("[write_lna] Calling core_write...\n")
  result <- core_write(x = x, transforms = transforms,
                       transform_params = transform_params,
                       mask = mask, header = header, plugins = plugins, run_id = run_id)
  cat("[write_lna] core_write returned. Plan object: ", class(result$plan), " Handle object: ", class(result$handle), "\n")

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

  cat(paste0("[write_lna] Attempting to open HDF5 file: ", file_to_use, " in_memory: ", in_memory, "\n"))
  h5 <- NULL # Initialize h5 to NULL
  tryCatch({
    if (in_memory) {
      cat("[write_lna] Opening HDF5 for in-memory via core driver.\n")
      h5 <- open_h5(file_to_use, mode = "w", driver = "core", backing_store = FALSE)
    } else {
      cat("[write_lna] Opening HDF5 for disk-based write.\n")
      h5 <- open_h5(file_to_use, mode = "w")
    }
    cat(paste0("[write_lna] open_h5 call completed. H5 object valid: ", ifelse(!is.null(h5) && h5$is_valid, "TRUE", "FALSE"), "\n"))
  }, error = function(e) {
    cat(paste0("[write_lna] ERROR during open_h5: ", conditionMessage(e), "\n"))
    # To ensure h5 is NULL if open_h5 failed before assignment or with an invalid object
    h5 <<- NULL 
    stop(e) # Re-throw the error to halt execution as expected
  })
  
  # If h5 is NULL here, it means open_h5 failed and error was re-thrown, 
  # or it failed in a way that didn't assign to h5 before erroring and stop() was called.
  # However, if we want to be super defensive for the on.exit:
  if (is.null(h5) || !inherits(h5, "H5File") || !h5$is_valid) {
     cat("[write_lna] HDF5 handle is NULL or invalid after open_h5 attempt and before materialise_plan. Aborting write_lna.\n")
     # Depending on desired behavior, could return an error or a specific result indicating failure.
     # For now, let it proceed to on.exit if h5 is NULL, close_h5_safely handles NULL.
     # But if it should stop, this is a place.
  }

  cat("[write_lna] Setting up on.exit handler to close HDF5 file.\n")
  on.exit({
    cat(paste0("[write_lna] on.exit: Attempting to close HDF5 file: ", file_to_use, "\n"))
    cat(paste0("[write_lna] on.exit: Is h5 object NULL? ", is.null(h5), "\n"))
    if (!is.null(h5)) {
      cat(paste0("[write_lna] on.exit: Is h5 valid before close? ", h5$is_valid, "\n"))
    }
    neuroarchive:::close_h5_safely(h5)
    cat(paste0("[write_lna] on.exit: close_h5_safely completed for ", file_to_use, "\n"))
    # Check if file exists after attempted close, especially if in_memory was false
    if (!in_memory && !is.null(file_to_use)) {
        cat(paste0("[write_lna] on.exit: Checking file existence for ", file_to_use, " after close: ", file.exists(file_to_use), "\n"))
    }
  }, add = TRUE)

  # plugins_from_handle <- result$handle$meta$plugins # Original line
  # Using a safer access pattern in case result$handle$meta or plugins is NULL
  plugins_from_handle <- list()
  if (!is.null(result) && !is.null(result$handle) && !is.null(result$handle$meta) && !is.null(result$handle$meta$plugins)) {
    plugins_from_handle <- result$handle$meta$plugins
  }
  if (length(plugins_from_handle) == 0) plugins_from_handle <- NULL

  header_from_handle <- list()
  if (!is.null(result) && !is.null(result$handle) && !is.null(result$handle$meta) && !is.null(result$handle$meta$header)) {
    header_from_handle <- result$handle$meta$header
  }

  cat("[write_lna] Calling materialise_plan...\n")
  materialise_plan(h5, result$plan,
                   checksum = checksum, # Pass the checksum argument
                   header = header_from_handle, # use safe version
                   plugins = plugins_from_handle) # use safe version
  cat("[write_lna] materialise_plan returned.\n")

  if (!is.null(block_table)) {
    cat("[write_lna] Writing block_table dataset.\n")
    bt_matrix <- as.matrix(block_table)
    h5_write_dataset(h5[["/"]], "spatial/block_table", bt_matrix)
    cat("[write_lna] Finished writing block_table dataset.\n")
  }

  final_out_file <- if (in_memory) NULL else file_to_use
  cat(paste0("[write_lna] Final out file determination: ", ifelse(is.null(final_out_file), "NULL", final_out_file), "\n"))
  
out <- list(file = final_out_file, plan = result$plan,
              header = header_from_handle)
  class(out) <- c("lna_write_result", class(out))
  cat("[write_lna] Exiting successfully.\n")
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