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
  stopifnot(is.character(transforms))
  stopifnot(is.list(transform_params))

  req_dims <- determine_required_dims(transforms, mask)
  x <- validate_input_data(x, min_dims = req_dims)

  mask_info <- validate_mask(mask, x)
  header_list <- validate_named_list(header, "header")
  plugin_list <- validate_named_list(plugins, "plugins")

  run_ids <- determine_run_ids(x, run_id)
  current_run_id <- if (length(run_ids) > 0) run_ids[1] else sanitize_run_id("run-01")

  x <- convert_inputs(x)
  if (!is.null(mask_info$array)) {
    check_mask_against_input(mask_info$array, mask_info$active_voxels, x)
  }

  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(input = x),
    initial_meta = list(mask = mask_info$array, header = header_list,
                        plugins = plugin_list),
    plan = plan,
    run_ids = run_ids,
    current_run_id = current_run_id,
    mask_info = if (!is.null(mask_info$array)) list(
      mask = mask_info$array,
      active_voxels = mask_info$active_voxels
    ) else NULL
  )

  merged_params <- resolve_transform_params(transforms, transform_params)
  handle <- run_transform_loop(handle, transforms, merged_params)

  if (length(transforms) == 0) {
    add_initial_data(handle)
  }

  list(handle = handle, plan = plan)
}

# Helper utilities -------------------------------------------------------

determine_required_dims <- function(transforms, mask) {
  first_type <- if (length(transforms) > 0) transforms[[1]] else ""
  req <- transform_min_dims(first_type)
  if (!is.null(mask)) req <- max(req, 3L)
  req
}

determine_run_ids <- function(x, run_id) {
  if (!is.null(run_id)) {
    return(vapply(run_id, sanitize_run_id, character(1), USE.NAMES = FALSE))
  }

  if (is.list(x)) {
    raw <- names(x)
    if (is.null(raw) || any(raw == "")) {
      raw <- sprintf("run-%02d", seq_along(x))
    }
    ids <- vapply(raw, sanitize_run_id, character(1), USE.NAMES = FALSE)
  } else {
    ids <- sanitize_run_id("run-01")
  }

  if (length(ids) == 0 && is.list(x) && length(x) == 0) {
    ids <- sanitize_run_id("run-01")
  }
  ids
}

convert_inputs <- function(x) {
  if (is.list(x)) {
    lapply(x, ensure_lna_array_input)
  } else {
    ensure_lna_array_input(x)
  }
}

check_mask_against_input <- function(mask_array, active_voxels, x) {
  check_one <- function(obj) {
    dims <- dim(obj)
    if (is.null(dims) || length(dims) < 3) {
      abort_lna(
        "input data must have at least 3 dimensions for mask check",
        .subclass = "lna_error_validation",
        location = "core_write:mask_check"
      )
    }
    if (active_voxels != prod(dims[1:3])) {
      abort_lna(
        "mask voxel count mismatch",
        .subclass = "lna_error_validation",
        location = "core_write:mask_check"
      )
    }
  }
  if (is.list(x)) {
    lapply(x, check_one)
  } else {
    check_one(x)
  }
}

run_transform_loop <- function(handle, transforms, params) {
  if (length(transforms) == 0) return(handle)

  progress_enabled <- is_progress_globally_enabled()
  step_loop <- function(h) {
    p <- if (progress_enabled) progressr::progressor(steps = length(transforms)) else NULL
    current_key <- "input"
    for (type in transforms) {
      if (!is.null(p)) p(message = type)
      step_idx <- h$plan$next_index
      output_key <- paste0(type, "_s", step_idx, "_out")
      desc <- list(
        type = type,
        params = params[[type]],
        inputs = list(current_key),
        outputs = list(output_key)
      )
      h <- run_transform_step("forward", type, desc, h, step_idx)
      current_key <- output_key
    }
    h
  }
  if (progress_enabled) {
    progressr::with_progress(step_loop(handle))
  } else {
    step_loop(handle)
  }
}

add_initial_data <- function(handle) {
  current_input <- handle$stash$input
  plan <- handle$plan
  add_payload <- function(data, run_label) {
    payload_key <- paste0(run_label, "_initial_data_payload")
    plan$add_payload(payload_key, data, overwrite = TRUE)
    plan$add_dataset_def(
      path = file.path("/scans", run_label, "data", "values"),
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
    if (length(handle$run_ids) == length(current_input)) {
      for (i in seq_along(handle$run_ids)) {
        add_payload(current_input[[i]], handle$run_ids[i])
      }
    } else {
      warning("Mismatch between number of runs in handle$run_ids and current_input list.")
    }
  } else {
    add_payload(current_input, handle$current_run_id)
  }
  invisible(handle)
}

#' Validate and normalise mask argument
#'
#' @param mask Optional mask object.
#' @param input Optional input object used to compare mask space/orientation.
#' @return List with `array` and `active_voxels` entries.
#' @keywords internal
validate_mask <- function(mask, input = NULL) {
  if (is.null(mask)) {
    return(list(array = NULL, active_voxels = NULL))
  }

  if (inherits(mask, "LogicalNeuroVol")) {
    mask_space <- tryCatch(space(mask), error = function(e) NULL)

    if (!is.null(input)) {
      x_run <- if (is.list(input)) input[[1]] else input
      input_space <- tryCatch(space(x_run), error = function(e) NULL)
      if (!is.null(mask_space) && !is.null(input_space)) {
        mask_dims <- tryCatch(dim(mask_space)[1:3], error = function(e) NULL)
        input_dims <- tryCatch(dim(input_space)[1:3], error = function(e) NULL)
        mask_trans <- tryCatch(trans(mask_space), error = function(e) NULL)
        input_trans <- tryCatch(trans(input_space), error = function(e) NULL)

        if (!is.null(mask_dims) && !is.null(input_dims) &&
            (!identical(mask_dims, input_dims) ||
             (!is.null(mask_trans) && !is.null(input_trans) && !identical(mask_trans, input_trans)))) {
          warn_lna("Mask orientation/space differs from input data; reslice mask to data space for accurate application.")
        }
      }
    }

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

#' Ensure LNA-Compatible Array Input
#'
#' Detects common `neuroim2` objects and converts them to 4D
#' arrays expected by the LNA engine. If the object is a 3D array
#' (or `DenseNeuroVol`), a singleton fourth dimension is appended and
#' an attribute `lna.was_3d` is set to `TRUE`.
#'
#' @param obj Input object.
#' @return A 4D array with possible `lna.was_3d` attribute.
#' @keywords internal
ensure_lna_array_input <- function(obj) {
  if (inherits(obj, "DenseNeuroVec")) {
    arr <- as.array(obj)
    attr(arr, "lna.was_3d") <- FALSE
    return(arr)
  }

  if (inherits(obj, "DenseNeuroVol")) {
    arr <- as.array(obj)
    arr <- array(arr, dim = c(dim(arr), 1L))
    attr(arr, "lna.was_3d") <- TRUE
    return(arr)
  }

  if (is.array(obj)) {
    d <- dim(obj)
    if (length(d) == 4) {
      attr(obj, "lna.was_3d") <- FALSE
      return(obj)
    }
    if (length(d) == 3) {
      arr <- array(obj, dim = c(d, 1L))
      attr(arr, "lna.was_3d") <- TRUE
      return(arr)
    }
  }

  abort_lna(
    "input must be DenseNeuroVec, DenseNeuroVol, or 3D/4D array",
    .subclass = "lna_error_validation",
    location = "ensure_lna_array_input"
  )
}
