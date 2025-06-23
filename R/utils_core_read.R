#' Helper utilities for core_read
#'
#' These functions modularize pieces of the core_read
#' implementation for clarity.
#'
#' @keywords internal
NULL

#' Normalize allow_plugins argument
#'
#' Ensures prompt mode degrades to installed when not interactive.
#'
#' @keywords internal
normalize_allow_plugins <- function(choice) {
  choice <- match.arg(choice, c("installed", "none", "prompt"))
  if (identical(choice, "prompt") && !rlang::is_interactive()) {
    choice <- "installed"
  }
  choice
}

#' Resolve runs for reading
#'
#' Determines which run identifiers to process and handles lazy mode.
#'
#' @keywords internal
resolve_runs_for_read <- function(h5, run_id, lazy) {
  available <- discover_run_ids(h5)
  runs <- resolve_run_ids(run_id, available)
  if (length(runs) == 0) {
    abort_lna("run_id did not match any runs", .subclass = "lna_error_run_id")
  }
  if (lazy && length(runs) > 1) {
    warning("Multiple runs matched; using first match in lazy mode")
    runs <- runs[1]
  }
  runs
}

#' Collect subset parameters for DataHandle
#'
#' Validates `roi_mask` and `time_idx` parameters and returns a list
#' suitable for the DataHandle constructor.
#'
#' @keywords internal
collect_subset_params <- function(roi_mask, time_idx) {
  subset <- list()
  if (!is.null(roi_mask)) {
    if (inherits(roi_mask, "LogicalNeuroVol")) {
      roi_mask <- as.array(roi_mask)
    }
    if (!(is.logical(roi_mask) && (is.vector(roi_mask) ||
                                   (is.array(roi_mask) && length(dim(roi_mask)) == 3)))) {
      abort_lna(
        "roi_mask must be logical vector or 3D logical array",
        .subclass = "lna_error_validation",
        location = "core_read:roi_mask"
      )
    }
    subset$roi_mask <- roi_mask
  }
  if (!is.null(time_idx)) {
    if (!is.numeric(time_idx)) {
      abort_lna(
        "time_idx must be numeric",
        .subclass = "lna_error_validation",
        location = "core_read:time_idx"
      )
    }
    subset$time_idx <- as.integer(time_idx)
  }
  subset
}

#' Prepare transforms for inverse pass
#'
#' Discovers descriptors and filters those without an implementation,
#' respecting the `allow_plugins` policy.
#'
#' @keywords internal
prepare_transforms_for_read <- function(tf_group, allow_plugins, file) {
  transforms <- discover_transforms(tf_group)

  missing <- transforms$type[
    vapply(
      transforms$type,
      function(t) is.null(getS3method("invert_step", t, optional = TRUE)),
      logical(1)
    )
  ]
  skip <- handle_missing_methods(
    missing,
    allow_plugins,
    location = sprintf("core_read:%s", file)
  )
  if (length(skip) > 0) {
    transforms <- transforms[!transforms$type %in% skip, , drop = FALSE]
  }
  transforms
}

#' Apply transforms in reverse order
#'
#' Runs the inverse transform chain on a `DataHandle`.
#'
#' @keywords internal
apply_invert_transforms <- function(handle, transforms, tf_group, validate, h5) {
  step_loop <- function(p) {
    h <- handle
    for (i in rev(seq_len(nrow(transforms)))) {
      if (!is.null(p)) p(message = transforms$type[[i]])
      name <- transforms$name[[i]]
      type <- transforms$type[[i]]
      step_idx <- transforms$index[[i]]
      desc <- read_json_descriptor(tf_group, name)

      h <- run_transform_step("invert", type, desc, h, step_idx)

      if (validate) runtime_validate_step(type, desc, h5)
    }
    h
  }

  with_progress_loop(nrow(transforms), step_loop)
}

#' Finalize handle for return
#'
#' Ensures dtype support and stores metadata.
#'
#' @keywords internal
finalize_handle_for_read <- function(handle, output_dtype, allow_plugins, file) {
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

#' Process a single run for core_read
#'
#' Constructs a `DataHandle`, applies transforms and finalizes the result.
#'
#' @keywords internal
process_run_core_read <- function(rid, h5, runs, subset_params, transforms,
                                 tf_group, validate, output_dtype, allow_plugins,
                                 file) {
  handle <- DataHandle$new(
    h5 = h5,
    subset = subset_params,
    run_ids = runs,
    current_run_id = rid
  )

  if (nrow(transforms) > 0) {
    handle <- apply_invert_transforms(handle, transforms, tf_group, validate, h5)
  } else {
    root <- h5[["/"]]
      on.exit(safe_h5_close(root), add = TRUE)
    path <- file.path("scans", rid, "data", "values")
    data <- h5_read(root, path)
    handle <- handle$with(stash = list(input = data))
  }

  finalize_handle_for_read(handle, output_dtype, allow_plugins, file)
}

