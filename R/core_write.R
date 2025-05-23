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
#' @param run_id Optional run_id parameter to override names(x) logic.
#'
#' @return A list with `handle` and `plan` objects.
#' @keywords internal
core_write <- function(x, transforms, transform_params = list(),
                       mask = NULL, header = NULL, plugins = NULL, run_id = NULL) {
  # cat("\n[core_write] Entry\n")
  # cat(paste0("[core_write] Number of transforms: ", length(transforms), "\n"))
  # cat(paste0("[core_write] Input data 'x' class: ", class(x), " Is list? ", is.list(x), "\n"))

  stopifnot(is.character(transforms))
  stopifnot(is.list(transform_params))

  # Determine required dimensionality from first transform
  first_type <- if (length(transforms) > 0) transforms[[1]] else ""
  req_dims <- transform_min_dims(first_type)
  if (!is.null(mask)) req_dims <- max(req_dims, 3L)

  # cat("[core_write] Validating input data...\n")
  x <- validate_input_data(x, min_dims = req_dims)
  # cat("[core_write] Input data validated.\n")

  # cat("[core_write] Validating mask...\n")
  mask_info <- validate_mask(mask)
  mask_array <- mask_info$array
  active_voxels <- mask_info$active_voxels
  # cat(paste0("[core_write] Mask validated. Is mask_array NULL? ", is.null(mask_array), " Active voxels: ", ifelse(is.null(active_voxels), "N/A", active_voxels), "\n"))

  header_list <- validate_named_list(header, "header")
  plugin_list <- validate_named_list(plugins, "plugins")
  # cat(paste0("[core_write] Header list length: ", length(header_list), " Plugin list length: ", length(plugin_list), "\n"))

  # cat("[core_write] Determining run identifiers...\n")
  if (!is.null(run_id)) {
    # cat(paste0("[core_write] run_id parameter provided: ", paste(run_id, collapse=", "), "\n"))
    if (is.list(x)) {
        sanitized_run_ids <- vapply(run_id, sanitize_run_id, character(1), USE.NAMES = FALSE)
        # cat(paste0("[core_write] Using provided run_id for list input, sanitized: ", paste(sanitized_run_ids, collapse=", "), "\n"))
    } else {
        sanitized_run_ids <- vapply(run_id, sanitize_run_id, character(1), USE.NAMES = FALSE)
        # cat(paste0("[core_write] Using provided run_id for non-list input, sanitized: ", paste(sanitized_run_ids, collapse=", "), "\n"))
    }
    final_run_ids <- sanitized_run_ids 
  } else if (is.list(x)) {
    # cat("[core_write] Input x is a list, deriving run_ids from names(x).\n")
    raw_run_ids <- names(x)
    if (is.null(raw_run_ids) || any(raw_run_ids == "")) {
      # cat("[core_write] names(x) are NULL or empty, generating default run_ids.\n")
      raw_run_ids <- sprintf("run-%02d", seq_along(x))
    }
    final_run_ids <- vapply(raw_run_ids, sanitize_run_id, character(1), USE.NAMES = FALSE)
    # cat(paste0("[core_write] Derived run_ids from list: ", paste(final_run_ids, collapse=", "), "\n"))
  } else {
    # cat("[core_write] Input x is not a list, using default run_id 'run-01'.\n")
    final_run_ids <- sanitize_run_id("run-01")
    # cat(paste0("[core_write] Default run_id: ", final_run_ids, "\n"))
  }
  
  if (length(final_run_ids) == 0 && is.list(x) && length(x) == 0) {
      # cat("[core_write] x is an empty list and no run_ids derived, defaulting to run-01 for handle continuity.\n")
      final_run_ids <- sanitize_run_id("run-01") 
  }

  current_run_id_for_handle <- if (length(final_run_ids) > 0) final_run_ids[1] else sanitize_run_id("run-01")
  # cat(paste0("[core_write] current_run_id_for_handle set to: ", current_run_id_for_handle, "\n"))

  if (!is.null(mask_array)) {
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

  # cat("[core_write] Creating Plan and DataHandle objects...\n")
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(input = x),
    initial_meta = list(mask = mask_array, header = header_list,
                        plugins = plugin_list),
    plan = plan,
    run_ids = final_run_ids, 
    current_run_id = current_run_id_for_handle, 
    mask_info = if (!is.null(mask_array)) list(mask = mask_array,
                                             active_voxels = active_voxels) else NULL
  )
  # cat("[core_write] Plan and DataHandle created.\n")

  # cat("[core_write] Resolving transform parameters...\n")
  merged_params <- resolve_transform_params(transforms, transform_params)
  # cat("[core_write] Transform parameters resolved.\n")

  # cat("[core_write] Starting transform loop...\n")
  progress_enabled <- is_progress_globally_enabled()
  loop <- function() {
    p <- if (progress_enabled) progressr::progressor(steps = length(transforms)) else NULL
    current_input_key <- "input" # Initial data key

    for (i in seq_along(transforms)) {
      type <- transforms[i]
      # cat(paste0("[core_write] Applying transform: ", type, "\n"))
      if (!is.null(p)) p(message = type)
      
      step_idx <- handle$plan$next_index # Get index before it's potentially incremented by add_descriptor
      
      # Define output key for this step
      # Ensure it's unique enough if a transform type is repeated, though full plan should handle this better.
      output_key <- paste0(type, "_s", step_idx, "_out") 

      desc <- list(
        type = type, 
        params = merged_params[[type]],
        inputs = list(current_input_key),
        outputs = list(output_key) # Define what this step will produce
      )
      
      # The run_transform_step will use desc$inputs and should update the stash with desc$outputs
      handle <<- run_transform_step("forward", type, desc, handle, step_idx)
      
      # The output of this step becomes the input for the next
      current_input_key <- output_key
      
      # cat(paste0("[core_write] Finished transform: ", type, "\n"))
    }
  }
  if (progress_enabled) {
    progressr::with_progress(loop())
  } else {
    loop()
  }
  # cat("[core_write] Transform loop finished.\n")

  if (length(transforms) == 0) {
    # cat("[core_write] No transforms specified, adding initial data to plan.\n")
    current_input <- handle$stash$input
    
    add_initial_data_to_plan <- function(data_to_write, run_label, plan_to_update) {
      # cat(paste0("[core_write] add_initial_data_to_plan for run_label: ", run_label, "\n"))
      payload_key <- paste0(run_label, "_initial_data_payload")
      plan_to_update$add_payload(payload_key, data_to_write, overwrite = TRUE)
      hdf5_path <- file.path("/scans", run_label, "data", "values")
      # cat(paste0("[core_write] Adding dataset definition for path: ", hdf5_path, "\n"))
      plan_to_update$add_dataset_def(
        path = hdf5_path,
        role = "raw_data",
        producer = "core_write_initial_input",
        origin = run_label,
        step_index = 0L,
        params_json = "{}",
        payload_key = payload_key,
        write_mode = "eager",
        dtype = NA_character_
      )
    }

    if (is.list(current_input) && !is.null(names(current_input))) {
      # cat("[core_write] Processing initial data for a list input.\n")
      if (length(handle$run_ids) == length(current_input)) {
          for (i in seq_along(handle$run_ids)) {
            add_initial_data_to_plan(current_input[[i]], handle$run_ids[i], handle$plan)
          }
      } else {
        # cat("[core_write] WARNING: Mismatch between number of runs in handle$run_ids and current_input list.\n")
        warning("Mismatch between number of runs in handle$run_ids and current_input list.")
      }
    } else {
      # cat("[core_write] Processing initial data for a single array input.\n")
      add_initial_data_to_plan(current_input, handle$current_run_id, handle$plan)
    }
    # cat("[core_write] Finished adding initial data to plan.\n")
  }

  # cat("[core_write] Exiting successfully.\n")
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
#' Ensures that `x` (or each element of a list `x`) has at least
#' `min_dims` dimensions. This prevents ambiguous handling of 2D
#' matrices when transforms expect 3D input.
#'
#' @param x Input object for `core_write`.
#' @param min_dims Integer, minimum required number of dimensions.
#'
#' @keywords internal
validate_input_data <- function(x, min_dims = 3L) {
  check_dims <- function(obj) {
    dims <- dim(obj)
    if (is.null(dims)) {
      base_dim <- length(obj)
      new_dims <- c(base_dim, rep(1L, max(0L, min_dims - 1L)))
      return(array(obj, dim = new_dims))
    }

    if (length(dims) < min_dims) {
      new_dims <- c(dims, rep(1L, min_dims - length(dims)))
      obj <- array(obj, dim = new_dims)
    }
    obj
  }

  if (is.list(x)) {
    lapply(x, check_dims)
  } else {
    check_dims(x)
  }
}
