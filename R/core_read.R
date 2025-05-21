#' Core LNA Read Routine
#'
#' @description Opens an LNA HDF5 file, discovers available transform
#'   descriptors and runs the inverse transform chain. This is a
#'   minimal skeleton used during early development.
#'
#' @param file Path to an LNA file on disk.
#' @param run_id Character vector of run identifiers or glob patterns. If
#'   `NULL`, the first available run is used. Glob patterns are matched
#'   against available run groups under `/scans`.
#' @param allow_plugins Character. How to handle transforms requiring
#'   external packages. One of "installed" (default), "none", or "prompt".
#'   "installed" attempts to use a transform if its implementation is
#'   available, issuing a warning when missing. "prompt" in interactive
#'   sessions asks whether to continue when a transform implementation is
#'   missing (non-affirmative answers abort). Non-interactive sessions treat
#'   "prompt" the same as "installed". "none" errors immediately when an
#'   implementation is missing.
#' @param validate Logical flag indicating if validation should be
#'   performed via `validate_lna()` before reading.
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param lazy Logical. If `TRUE`, the HDF5 file handle remains open
#'   after return (for lazy reading).
#'
#' @return If a single run is selected, a `DataHandle` object. When
#'   multiple runs match and `lazy = FALSE`, a named list of `DataHandle`
#'   objects is returned.
#' @import hdf5r
#' @keywords internal
core_read <- function(file, run_id = NULL,
                      allow_plugins = c("installed", "none", "prompt"), validate = FALSE,
                      output_dtype = c("float32", "float64", "float16"),
                      roi_mask = NULL, time_idx = NULL,
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

  available_runs <- discover_run_ids(h5)
  runs <- resolve_run_ids(run_id, available_runs)
  if (length(runs) == 0) {
    abort_lna("run_id did not match any runs", .subclass = "lna_error_run_id")
  }
  if (lazy && length(runs) > 1) {
    warning("Multiple runs matched; using first match in lazy mode")
    runs <- runs[1]
  }

  if (validate) {
    validate_lna(file)
  }


  subset_params <- list()
  if (!is.null(roi_mask)) {
    if (inherits(roi_mask, "LogicalNeuroVol")) {
      roi_mask <- as.array(roi_mask)
    }
    subset_params$roi_mask <- roi_mask
  }
  if (!is.null(time_idx)) {
    subset_params$time_idx <- as.integer(time_idx)
  }

  handle <- DataHandle$new(h5 = h5, subset = subset_params)

  tf_group <- h5[["transforms"]]

  transforms <- discover_transforms(tf_group)

  missing_methods <- transforms$type[
    vapply(
      transforms$type,
      function(t) is.null(getS3method("invert_step", t, optional = TRUE)),
      logical(1)
    )
  ]
  skip_types <- handle_missing_methods(
    missing_methods,
    allow_plugins,
    location = sprintf("core_read:%s", file)
  )
  if (length(skip_types) > 0) {
    transforms <- transforms[!transforms$type %in% skip_types, , drop = FALSE]
  }

  process_run <- function(rid) {
    handle <- DataHandle$new(h5 = h5, run_ids = runs, current_run_id = rid)

    if (nrow(transforms) > 0) {
      progress_enabled <- !progressr::handlers_is_empty()
      loop <- function() {
        p <- if (progress_enabled) progressr::progressor(steps = nrow(transforms)) else NULL
        for (i in rev(seq_len(nrow(transforms)))) {
          if (!is.null(p)) p(message = transforms$type[[i]])
          name <- transforms$name[[i]]
          type <- transforms$type[[i]]
          step_idx <- transforms$index[[i]]
          desc <- read_json_descriptor(tf_group, name)
          handle <<- run_transform_step("invert", type, desc, handle, step_idx)
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
    handle
  }

  results <- lapply(runs, process_run)
  names(results) <- runs

  if (length(results) == 1) {
    results[[1]]
  } else {
    results
  }
}
