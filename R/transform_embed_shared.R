# embed.shared_basis transform
#
# Projects data onto a shared basis artifact stored externally. Only
# per-subject coefficients (and optional offset) are stored in the archive.
# The basis itself is resolved at read time via resolve_basis().

#' Forward step for embed.shared_basis transform
#'
#' Projects input data onto a pre-computed shared basis. The basis is loaded
#' from an external `.lna_basis.h5` file identified by `basis_id` +
#' `basis_checksum`. Only the projection coefficients and per-subject offset
#' are stored in the archive.
#'
#' @param type The transform type string.
#' @param desc The transform descriptor list.
#' @param handle The DataHandle object.
#' @return The updated DataHandle with coefficients in the stash.
#' @method forward_step embed.shared_basis
#' @export
#' @keywords internal
forward_step.embed.shared_basis <- function(type, desc, handle) {
  p <- desc$params %||% list()

  basis_id <- p$basis_id
  basis_checksum <- p$basis_checksum
  basis_locator <- p$basis_locator
  do_center <- p$center %||% TRUE

  if (is.null(basis_id) || !nzchar(basis_id)) {
    abort_lna("'basis_id' must be provided",
              .subclass = "lna_error_validation",
              location = "forward_step.embed.shared_basis:basis_id")
  }
  if (is.null(basis_checksum) || !nzchar(basis_checksum)) {
    abort_lna("'basis_checksum' must be provided",
              .subclass = "lna_error_validation",
              location = "forward_step.embed.shared_basis:checksum")
  }

  # Resolve and load shared basis
  basis_file <- resolve_basis(
    basis_id = basis_id,
    checksum = basis_checksum,
    locator = basis_locator,
    archive_path = handle$h5$filename
  )
  basis_data <- read_shared_basis(basis_file, verify_checksum = TRUE)
  .validate_shared_basis_mask(handle, basis_data, "forward_step.embed.shared_basis:mask")
  loadings <- basis_data$loadings  # voxels x atoms

  # Get input data
  inp <- handle$pull_first(c("input_dense_mat", "dense_mat", "input"))
  input_key <- inp$key
  X <- as_dense_mat(inp$value)  # time x voxels

  # Validate dimensions
  if (ncol(X) != nrow(loadings)) {
    abort_lna(
      sprintf("Voxel dimension mismatch: data has %d, basis has %d",
              ncol(X), nrow(loadings)),
      .subclass = "lna_error_validation",
      location = "forward_step.embed.shared_basis:dims"
    )
  }

  # Center if requested
  offset <- numeric(0)
  X_proj <- X
  if (isTRUE(do_center)) {
    offset <- colMeans(X)
    X_proj <- sweep(X, 2L, offset, "-")
  }

  # Project: solve L'L c = L' x for coefficients.
  gram <- crossprod(loadings)
  rhs <- X_proj %*% loadings  # time x atoms
  coeff <- t(.solve_shared_basis_gram(gram, t(rhs)))  # time x atoms

  # Build paths
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  scans_root <- lna_options("paths.scans_root")[[1]]
  coef_path <- paste0(scans_root, run_id, "/", base_name, "/coefficients")
  offset_path <- paste0(scans_root, run_id, "/", base_name, "/offset")
  step_index <- plan$next_index
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  # Build descriptor
  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  output_key <- desc$outputs[[1]] %||% "coefficients"
  desc$outputs <- c(output_key)

  datasets <- list(
    list(path = coef_path, role = "coefficients")
  )
  if (length(offset) > 0) {
    datasets[[length(datasets) + 1L]] <- list(path = offset_path, role = "offset")
  }
  desc$datasets <- datasets

  # Register in plan
  plan$add_descriptor(fname, desc)
  plan$add_payload(coef_path, as.matrix(coeff))
  plan$add_dataset_def(coef_path, "coefficients", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       coef_path, "eager", dtype = NA_character_)
  if (length(offset) > 0) {
    plan$add_payload(offset_path, offset)
    plan$add_dataset_def(offset_path, "offset", as.character(type), run_id,
                         as.integer(step_index), params_json,
                         offset_path, "eager", dtype = NA_character_)
  }

  handle$plan <- plan
  handle <- handle$update_stash(keys = input_key,
                      new_values = setNames(list(as.matrix(coeff)), output_key))
  return(handle)
}

#' Inverse step for embed.shared_basis transform
#'
#' Reconstructs dense data from coefficients using the shared basis resolved
#' from the descriptor's basis_id + checksum.
#'
#' @param type The transform type string.
#' @param desc The transform descriptor list.
#' @param handle The DataHandle object.
#' @return The updated DataHandle with reconstructed data in the stash.
#' @method invert_step embed.shared_basis
#' @export
#' @keywords internal
invert_step.embed.shared_basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_id <- p$basis_id
  basis_checksum <- p$basis_checksum
  basis_locator <- p$basis_locator
  do_center <- p$center %||% TRUE

  if (is.null(basis_id) || !nzchar(basis_id)) {
    abort_lna("'basis_id' missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.embed.shared_basis:basis_id")
  }
  if (is.null(basis_checksum) || !nzchar(basis_checksum)) {
    abort_lna("'basis_checksum' missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.embed.shared_basis:checksum")
  }

  # Resolve and load shared basis
  basis_file <- resolve_basis(
    basis_id = basis_id,
    checksum = basis_checksum,
    locator = basis_locator,
    archive_path = handle$h5$filename
  )
  basis_data <- read_shared_basis(basis_file, verify_checksum = TRUE)
  .validate_shared_basis_mask(handle, basis_data, "invert_step.embed.shared_basis:mask")
  loadings <- basis_data$loadings  # voxels x atoms

  # Get coefficients from stash
  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  input_key <- desc$inputs[[1]] %||% "input"

  if (handle$has_key(coeff_key)) {
    coeff <- handle$get_inputs(coeff_key)[[coeff_key]]  # time x atoms
  } else {
    coeff_entry <- Filter(function(d) identical(d$role, "coefficients"), desc$datasets %||% list())
    if (length(coeff_entry) == 0) {
      return(handle)
    }
    root <- handle$h5[["/"]]
    on.exit(safe_h5_close(root), add = TRUE)
    coeff <- h5_read(root, coeff_entry[[1]]$path)
  }

  # Apply subsetting
  subset <- handle$subset
  applied <- handle$subset_applied %||% list(roi_mask = FALSE, time_idx = FALSE)
  roi_mask <- subset$roi_mask %||% subset$roi
  time_idx <- subset$time_idx %||% subset$time

  if (!isTRUE(applied$time_idx) && !is.null(time_idx)) {
    idx <- as.integer(time_idx)
    idx <- idx[!is.na(idx) & idx >= 1L & idx <= nrow(coeff)]
    coeff <- coeff[idx, , drop = FALSE]
    handle$subset_applied$time_idx <- TRUE
  }

  # Subset loadings for ROI
  vox_idx <- NULL
  if (!is.null(roi_mask)) {
    roi_vec <- as.logical(roi_mask)
    n_vox <- nrow(loadings)

    # Map full-volume ROI through file mask if needed
    if (length(roi_vec) != n_vox && !is.null(handle$mask_info) &&
        !is.null(handle$mask_info$mask)) {
      file_mask <- as.logical(handle$mask_info$mask)
      if (length(roi_vec) == length(file_mask) && n_vox == sum(file_mask)) {
        roi_vec <- roi_vec[file_mask]
      }
    }

    if (length(roi_vec) != n_vox) {
      abort_lna(
        "roi_mask length does not match basis voxel dimension",
        .subclass = "lna_error_validation",
        location = "invert_step.embed.shared_basis:roi_mask"
      )
    }

    vox_idx <- which(roi_vec)
    loadings <- loadings[vox_idx, , drop = FALSE]
    handle$subset_applied$roi_mask <- TRUE
  }

  if (ncol(coeff) != ncol(loadings)) {
    abort_lna(
      sprintf("Coefficient dimension mismatch: coeff has %d columns, basis has %d atoms",
              ncol(coeff), ncol(loadings)),
      .subclass = "lna_error_runtime",
      location = "invert_step.embed.shared_basis:dims"
    )
  }

  # Reconstruct: dense = coeff %*% t(loadings)
  dense <- as.matrix(coeff) %*% t(loadings)

  # Add back offset
  if (isTRUE(do_center)) {
    # Read offset from HDF5 if available
    offset_entries <- Filter(function(d) identical(d$role, "offset"), desc$datasets)
    if (length(offset_entries) > 0) {
      root <- handle$h5[["/"]]
      offset <- h5_read(root, offset_entries[[1]]$path)
      if (!is.null(roi_mask) && !is.null(vox_idx)) {
        offset <- offset[vox_idx]
      }
      dense <- sweep(dense, 2L, offset, "+")
    }
  }

  handle <- handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
  return(handle)
}

#' Default parameters for embed.shared_basis
#' @export
#' @keywords internal
lna_default.embed.shared_basis <- function() {
  list(
    basis_id = "",
    basis_checksum = "",
    basis_kind = "",
    basis_locator = NULL,
    center = TRUE
  )
}

.validate_shared_basis_mask <- function(handle, basis_data, location) {
  aliases <- basis_data$basis_ids %||% character()
  asset <- basis_data$basis_id
  if (is.null(asset) && length(aliases) > 0L) asset <- aliases[[1L]]
  asset <- asset %||% basis_data$checksum %||% "<unknown-shared-basis>"
  fail <- function(message, field, expected = NA_character_,
                   actual = NA_character_) {
    abort_lna(
      message,
      asset = asset,
      field = field,
      expected_identity = expected,
      actual_identity = actual,
      .subclass = "lna_error_spatial_mismatch",
      location = location
    )
  }

  expected_grid <- basis_data$grid_id %||% ""
  expected_support <- basis_data$support_id %||% ""
  support_type <- basis_data$support_type %||% ""
  if (!nzchar(expected_grid)) {
    fail("Shared basis is missing required grid identity", "basis_grid_id")
  }
  if (!nzchar(expected_support)) {
    fail(
      "Shared basis is missing required support identity",
      "basis_support_id",
      expected = expected_grid
    )
  }
  if (!support_type %in% c("masked", "full-grid")) {
    fail(
      "Shared basis is missing a valid support descriptor",
      "basis_support_type",
      expected = expected_support,
      actual = support_type
    )
  }

  header <- handle$meta$header %||% list()
  if ((is.null(header$dims) || is.null(header$transform)) &&
      !is.null(handle$h5)) {
    h5_valid <- tryCatch(isTRUE(handle$h5$is_valid), error = function(e) FALSE)
    if (h5_valid) {
      header <- .read_lna_header(
        handle$h5,
        run_id = handle$current_run_id
      )
    }
  }
  if (is.null(header$dims)) {
    fail(
      "Input is missing dimensions required to validate shared basis grid",
      "input_dims",
      expected = expected_grid
    )
  }
  if (is.null(header$transform)) {
    fail(
      "Input is missing affine required to validate shared basis grid",
      "input_transform",
      expected = expected_grid
    )
  }

  actual_grid <- tryCatch(
    .lna_grid_identity_from_header(
      header,
      units = header$spatial_units %||% "mm",
      voxel_ordering = header$voxel_ordering %||% "first-axis-fastest",
      location = location
    ),
    error = function(e) {
      if (inherits(e, "lna_error_spatial_identity")) {
        fail(
          conditionMessage(e),
          e$field %||% "input_grid",
          expected = expected_grid
        )
      }
      stop(e)
    }
  )
  if (!identical(actual_grid$grid_id, expected_grid)) {
    fail(
      "Shared basis grid identity does not match input grid",
      "grid_id",
      expected = expected_grid,
      actual = actual_grid$grid_id
    )
  }

  current_mask <- if (!is.null(handle$mask_info)) {
    handle$mask_info$mask %||% NULL
  } else {
    NULL
  }
  if (identical(support_type, "masked") && is.null(current_mask)) {
    fail(
      "Masked shared basis requires explicit input support",
      "input_support",
      expected = expected_support
    )
  }
  if (is.null(current_mask)) {
    current_mask <- rep(TRUE, prod(actual_grid$dims))
  } else {
    current_mask <- as.array(current_mask)
    if (length(current_mask) != prod(actual_grid$dims)) {
      fail(
        "Input support length does not match its grid dimensions",
        "input_support",
        expected = expected_support,
        actual = paste0("length:", length(current_mask))
      )
    }
  }

  actual_support <- .lna_support_id(actual_grid$grid_id, current_mask)
  if (!identical(actual_support, expected_support)) {
    fail(
      "Shared basis support identity does not match input support",
      "support_id",
      expected = expected_support,
      actual = actual_support
    )
  }

  invisible(TRUE)
}

.solve_shared_basis_gram <- function(gram, rhs_t) {
  gram_mat <- as.matrix(gram)
  rhs_mat <- as.matrix(rhs_t)
  if (is.null(dim(rhs_mat))) {
    rhs_mat <- matrix(rhs_mat, ncol = 1L)
  }
  tryCatch(
    chol2inv(chol(gram_mat)) %*% rhs_mat,
    error = function(e_chol) {
      warning(
        "Cholesky decomposition failed for Gram matrix, falling back to base::solve: ",
        conditionMessage(e_chol),
        call. = FALSE
      )
      tryCatch(
        solve(gram_mat, rhs_mat),
        error = function(e_solve) {
          abort_lna(
            sprintf("Gram matrix solve failed: %s", conditionMessage(e_solve)),
            .subclass = "lna_error_runtime",
            location = ".solve_shared_basis_gram"
          )
        }
      )
    }
  )
}
