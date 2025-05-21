#' Core Write Pipeline (Skeleton)
#'
#' @description Internal function orchestrating the forward transform pass.
#'   This version implements the bare structure used for early milestones.
#'
#' @param x An input object to be written.
#' @param transforms Character vector of transform types (forward order).
#' @param transform_params Named list of user supplied parameters for transforms.
#' @param mask Optional mask object passed through to transforms.
#' @param header Optional named list of header attributes.
#' @param plugins Optional named list of plugin metadata to store under
#'   `/plugins`.
#'
#' @return A list with `handle` and `plan` objects.
#' @keywords internal
core_write <- function(x, transforms, transform_params = list(),
                       mask = NULL, header = NULL, plugins = NULL) {
  stopifnot(is.character(transforms))
  stopifnot(is.list(transform_params))

  validate_input_data(x)

  mask_info <- validate_mask(mask)
  mask_array <- mask_info$array
  active_voxels <- mask_info$active_voxels

  header_list <- validate_named_list(header, "header")

  plugin_list <- validate_named_list(plugins, "plugins")

  # --- Determine run identifiers ---
  if (is.list(x)) {
    raw_run_ids <- names(x)
    if (is.null(raw_run_ids) || any(raw_run_ids == "")) {
      raw_run_ids <- sprintf("run-%02d", seq_along(x))
    }
    run_ids <- vapply(raw_run_ids, sanitize_run_id, character(1), USE.NAMES = FALSE)
    # TODO: what if x was list() ? run_ids would be character(0). current_run_id would fail.
    # For now, assume x is non-empty if list, or is a single array.
  } else {
    run_ids <- sanitize_run_id("run-01") # Ensure this is sanitized
  }

  current_run_id_for_handle <- if (length(run_ids) > 0) run_ids[1] else sanitize_run_id("run-01")

  if (!is.null(mask_array)) {
    # Strict check: mask must cover all voxels of the input data. This
    # may be relaxed when partial masking is supported.
    check_mask <- function(obj) {
      dims <- dim(obj)
      if (is.null(dims) || length(dims) < 3) {
        abort_lna(
          "input data must have at least 3 dimensions for mask check",
          .subclass = "lna_error_validation",
          location = "core_write:mask_check"
        )
      }
      voxel_count <- prod(dims[1:3])
      if (active_voxels != voxel_count) {
        abort_lna(
          "mask voxel count mismatch",
          .subclass = "lna_error_validation",
          location = "core_write:mask_check"
        )
      }
    }
    if (is.list(x)) {
      lapply(x, check_mask)
    } else {
      check_mask(x)
    }
  }

  # --- Create plan and initial handle ---
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(input = x),
    initial_meta = list(mask = mask_array, header = header_list,
                        plugins = plugin_list),
    plan = plan,
    run_ids = run_ids, # Already sanitized if from list, or is the single sanitized one
    current_run_id = current_run_id_for_handle, # Use sanitized version
    mask_info = if (!is.null(mask_array)) list(mask = mask_array,
                                             active_voxels = active_voxels) else NULL
  )

  # --- Resolve parameters with defaults and package options ---
  merged_params <- resolve_transform_params(transforms, transform_params)

  # --- Loop through transforms calling forward_step ---
  progress_enabled <- is_progress_globally_enabled()
  loop <- function() {
    p <- if (progress_enabled) progressr::progressor(steps = length(transforms)) else NULL
    for (type in transforms) {
      if (!is.null(p)) p(message = type)
      desc <- list(type = type, params = merged_params[[type]])
      step_idx <- handle$plan$next_index
      handle <<- run_transform_step("forward", type, desc, handle, step_idx)
    }
  }
  if (progress_enabled) {
    progressr::with_progress(loop())
  } else {
    loop()
  }

  # If there were no transforms, the initial input data needs to be planned for writing.
  if (length(transforms) == 0) {
    current_input <- handle$stash$input # This could be a list if x was a list
    
    add_initial_data_to_plan <- function(data_to_write, run_label, plan_to_update) {
      payload_key <- paste0(run_label, "_initial_data_payload")
      plan_to_update$add_payload(payload_key, data_to_write, overwrite = TRUE)
      
      # Construct the HDF5 path for the data
      # Assumes sanitize_run_id has been called on run_label if it came from multi-run list names
      # For single run_ids determined by core_write, it's already sanitized.
      hdf5_path <- file.path("/scans", run_label, "data", "values")
      
      plan_to_update$add_dataset_def(
        path = hdf5_path,
        role = "raw_data",
        producer = "core_write_initial_input",
        origin = run_label,
        step_index = 0L, # Or -1L for initial data before transforms
        params_json = "{}",
        payload_key = payload_key,
        write_mode = "eager"
      )
    }

    if (is.list(current_input) && !is.null(names(current_input))) {
      # If x was a named list (multi-run) - names(current_input) are the original, unsanitized names
      # run_ids from handle are the sanitized ones. We need to iterate over sanitized run_ids.
      if (length(handle$run_ids) == length(current_input)) { # sanity check
          for (i in seq_along(handle$run_ids)) {
            # Use the sanitized run_id from handle$run_ids
            # And the corresponding data from current_input (which should align by position)
            add_initial_data_to_plan(current_input[[i]], handle$run_ids[i], handle$plan)
          }
      } else {
        # This case should ideally not happen if run_ids were derived from names(x)
        warning("Mismatch between number of runs in handle$run_ids and current_input list.")
      }
    } else {
      # If x was a single array (or unnamed list treated as single run)
      # Use handle$current_run_id which is guaranteed to be sanitized.
      add_initial_data_to_plan(current_input, handle$current_run_id, handle$plan)
    }
  }

  list(handle = handle, plan = plan)
}

#' Validate and normalise mask argument
#'
#' @param mask Optional mask object.
#' @return List with `array` and `active_voxels` entries.
#' @keywords internal
validate_mask <- function(mask) {
  if (is.null(mask)) {
    return(list(array = NULL, active_voxels = NULL))
  }

  if (inherits(mask, "LogicalNeuroVol")) {
    arr <- as.array(mask)
  } else if (is.array(mask) && length(dim(mask)) == 3 && is.logical(mask)) {
    arr <- mask
  } else {
    abort_lna(
      "mask must be LogicalNeuroVol or 3D logical array",
      .subclass = "lna_error_validation",
      location = "core_write:mask"
    )
  }

  list(array = arr, active_voxels = sum(arr))
}

#' Validate optional named lists
#'
#' Used for the `header` and `plugins` arguments in `core_write`.
#'
#' @param lst List or `NULL`.
#' @param field Field name used in error messages.
#' @return The validated list or an empty list if `NULL` or empty.
#' @keywords internal
validate_named_list <- function(lst, field) {
  if (is.null(lst)) {
    return(list())
  }

  stopifnot(is.list(lst))

  if (length(lst) == 0) {
    return(list())
  }

  if (is.null(names(lst)) || any(names(lst) == "")) {
    abort_lna(
      sprintf("%s must be a named list", field),
      .subclass = "lna_error_validation",
      location = sprintf("core_write:%s", field)
    )
  }
  lst
}

#' Validate input data
#'
#' Ensures that `x` (or each element of a list `x`) has at least three
#' dimensions. This prevents ambiguous handling of 2D matrices.
#'
#' @param x Input object for `core_write`.
#' @keywords internal
validate_input_data <- function(x) {
  check_dims <- function(obj) {
    dims <- dim(obj)
    if (is.null(dims) || length(dims) < 3) {
      abort_lna(
        "input data must have at least 3 dimensions",
        .subclass = "lna_error_validation",
        location = "core_write:input"
      )
    }
    invisible(TRUE)
  }

  if (is.list(x)) {
    lapply(x, check_dims)
  } else {
    check_dims(x)
  }

  invisible(TRUE)
}
