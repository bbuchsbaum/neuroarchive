#' Core LNA Read Routine
#'
#' @description Opens an LNA HDF5 file, discovers available transform
#'   descriptors and runs the inverse transform chain. This is a
#'   minimal skeleton used during early development.
#'
#' @param file Path to an LNA file on disk.
#' @param allow_plugins Character. How to handle transforms requiring
#'   external packages. One of "installed" (default), "none", or "prompt".
#'   "none" errors on missing implementations. "installed" warns and
#'   proceeds if a transform implementation is unavailable. "prompt"
#'   behaves like "installed" when \code{!rlang::is_interactive()}.
#' @param validate Logical flag indicating if validation should be
#'   performed via `validate_lna()` before reading.
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param lazy Logical. If `TRUE`, the HDF5 file handle remains open
#'   after return (for lazy reading).
#'
#' @return A `DataHandle` object representing the loaded data.
#' @import hdf5r
#' @keywords internal
core_read <- function(file, allow_plugins = c("installed", "none", "prompt"), validate = FALSE,
                      output_dtype = c("float32", "float64", "float16"),
                      lazy = FALSE) {
  allow_plugins <- match.arg(allow_plugins)
  if (identical(allow_plugins, "prompt") && !rlang::is_interactive()) {
    allow_plugins <- "installed"
  }
  output_dtype <- match.arg(output_dtype)
  h5 <- open_h5(file, mode = "r")
  if (!lazy) {
    on.exit(close_h5_safely(h5))
  }

  if (validate) {
    validate_lna(file)
  }

  handle <- DataHandle$new(h5 = h5)
  tf_group <- h5[["transforms"]]

  transforms <- discover_transforms(tf_group)

  missing_methods <- transforms$type[
    vapply(
      transforms$type,
      function(t) is.null(getS3method("invert_step", t, optional = TRUE)),
      logical(1)
    )
  ]
  if (length(missing_methods) > 0) {
    msg <- paste0(
      "Missing invert_step implementation for transform(s): ",
      paste(unique(missing_methods), collapse = ", ")
    )
    if (identical(allow_plugins, "none")) {
      abort_lna(msg, .subclass = "lna_error_no_method")
    } else {
      warning(msg)
    }
  }

  if (nrow(transforms) > 0) {
    progress_enabled <- !progressr::handlers_is_empty()
    loop <- function() {
      p <- if (progress_enabled) progressr::progressor(steps = nrow(transforms)) else NULL
      for (i in rev(seq_len(nrow(transforms)))) {
        if (!is.null(p)) p(message = transforms$type[[i]])
        name <- transforms$name[[i]]
        type <- transforms$type[[i]]
        desc <- read_json_descriptor(tf_group, name)
        handle <<- invert_step(type, desc, handle)
      }
    }
    if (progress_enabled) {
      progressr::with_progress(loop())
    } else {
      loop()
    }
  }

  if (identical(output_dtype, "float16") && !has_float16_support()) {
    abort_lna(
      "float16 output not supported",
      .subclass = "lna_error_float16_unsupported",
      location = sprintf("core_read:%s", file)
    )
  }
  handle$meta$output_dtype <- output_dtype
  handle$meta$allow_plugins <- allow_plugins

  return(handle)
}
