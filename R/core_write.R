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
#'
#' @return A list with `handle` and `plan` objects.
#' @keywords internal
core_write <- function(x, transforms, transform_params = list(),
                       mask = NULL, header = NULL) {
  stopifnot(is.character(transforms))
  stopifnot(is.list(transform_params))

  # --- Validate and store mask ---
  mask_array <- NULL
  active_voxels <- NULL
  if (!is.null(mask)) {
    if (inherits(mask, "LogicalNeuroVol")) {
      mask_array <- as.array(mask)
    } else if (is.array(mask) && length(dim(mask)) == 3 && is.logical(mask)) {
      mask_array <- mask
    } else {
      abort_lna("mask must be LogicalNeuroVol or 3D logical array",
                .subclass = "lna_error_validation")
    }
    active_voxels <- sum(mask_array)
  }
  if (!is.null(header)) {
    stopifnot(is.list(header))
    if (is.null(names(header)) || any(names(header) == "")) {
      abort_lna("header must be a named list", .subclass = "lna_error_validation")
    }
    header_list <- header
  } else {
    header_list <- list()
  }

  # --- Determine run identifiers ---
  if (is.list(x)) {
    if (is.null(names(x)) || any(names(x) == "")) {
      names(x) <- sprintf("run-%02d", seq_along(x))
    }
    run_ids <- names(x)
  } else {
    run_ids <- "run-01"
  }

  if (!is.null(mask_array)) {
    check_mask <- function(obj) {
      dims <- dim(obj)
      if (is.null(dims) || length(dims) < 3) {
        abort_lna("input data must have at least 3 dimensions for mask check",
                  .subclass = "lna_error_validation")
      }
      voxel_count <- prod(dims[1:3])
      if (active_voxels != voxel_count) {
        abort_lna("mask voxel count mismatch", .subclass = "lna_error_validation")
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
    initial_meta = list(mask = mask_array, header = header_list),
    plan = plan,
    run_ids = run_ids,
    current_run_id = run_ids[1],
    mask_info = if (!is.null(mask_array)) list(mask = mask_array,
                                             active_voxels = active_voxels) else NULL
  )

  # --- Resolve parameters with defaults and package options ---
  merged_params <- resolve_transform_params(transforms, transform_params)

  # --- Loop through transforms calling forward_step ---
  for (type in transforms) {
    desc <- list(type = type, params = merged_params[[type]])
    cls <- structure(type, class = c(type, "character"))
    handle <- forward_step(cls, desc, handle)
  }

  list(handle = handle, plan = plan)
}
