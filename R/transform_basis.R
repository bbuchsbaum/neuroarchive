#' Basis Transform - Inverse Step
#'
#' Reconstructs data from coefficients using a stored basis matrix.
#' The `basis` dataset may be stored either as a component-by-voxel matrix
#' (`storage_order = "component_x_voxel"`) or as a voxel-by-component matrix
#' (`storage_order = "voxel_x_component"`). After optional transposition, it
#' should have dimensions `n_voxel x n_component` for reconstruction.
#' @param type the transform type string
#' @param desc the transform descriptor
#' @param handle the LNA data handle
#' @return the updated LNA data handle
#' @method invert_step basis
#' @export
#' @keywords internal
invert_step.basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  storage_order <- p$storage_order %||% "component_x_voxel"
  allowed_orders <- c("component_x_voxel", "voxel_x_component")
  if (!storage_order %in% allowed_orders) {
    abort_lna(
      sprintf(
        "Invalid storage_order '%s'. Allowed values: %s",
        storage_order,
        paste(allowed_orders, collapse = ", ")
      ),
      .subclass = "lna_error_validation",
      location = "invert_step.basis:storage_order"
    )
  }

  basis_path <- NULL
  if (!is.null(desc$datasets) && length(desc$datasets) > 0) {
    idx <- which(vapply(desc$datasets, function(d) d$role, character(1)) == "basis_matrix")
    if (length(idx) > 0) basis_path <- desc$datasets[[idx[1]]]$path
  }
  if (is.null(basis_path)) {
    abort_lna(
      "basis_matrix path not found in descriptor",
      .subclass = "lna_error_descriptor",
      location = "invert_step.basis"
    )
  }

  # In the default pipeline this step *produces* the basis matrix during the
  # forward pass, while the following step (typically `embed`) produces the
  # coefficient matrix. On inversion, the preceding step (e.g. `invert_step.embed`)
  # may already have reconstructed the dense matrix and stored it under
  # `desc$outputs[[1]]`. We therefore treat `desc$outputs[[1]]` as the *incoming*
  # key and `desc$inputs[[1]]` as the *target* key.
  incoming_key <- desc$outputs[[1]] %||% "coefficients"
  target_key <- desc$inputs[[1]] %||% "input"

  if (!handle$has_key(incoming_key)) {
    # Nothing to reconstruct; return handle unchanged
    return(handle)
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)

  incoming <- handle$get_inputs(incoming_key)[[incoming_key]]

  subset <- handle$subset
  roi_mask <- subset$roi_mask %||% subset$roi
  # Determine K (components) and V (voxels) from stored basis matrix.
  k_components <- if (identical(storage_order, "component_x_voxel")) nrow(basis) else ncol(basis)
  v_voxels <- if (identical(storage_order, "component_x_voxel")) ncol(basis) else nrow(basis)
  v_voxels_effective <- v_voxels
  roi_vec <- NULL
  if (!is.null(roi_mask)) {
    roi_vec <- as.logical(roi_mask)
    if (length(roi_vec) != v_voxels && !is.null(handle$mask_info) && !is.null(handle$mask_info$mask)) {
      file_mask <- as.logical(handle$mask_info$mask)
      if (length(roi_vec) == length(file_mask) && v_voxels == sum(file_mask)) {
        roi_vec <- roi_vec[file_mask]
      }
    }
    if (length(roi_vec) != v_voxels) {
      abort_lna(
        "roi_mask length does not match basis voxel dimension",
        .subclass = "lna_error_validation",
        location = "invert_step.basis:roi_mask"
      )
    }
    v_voxels_effective <- sum(roi_vec)
  }

  dense <- NULL

  applied <- handle$subset_applied %||% list(roi_mask = FALSE, time_idx = FALSE)
  if (is.matrix(incoming) && ncol(incoming) == v_voxels_effective) {
    # Already reconstructed dense matrix (time x voxel); just pass through.
    dense <- incoming
    if (!isTRUE(applied$time_idx)) {
      time_idx <- subset$time_idx %||% subset$time
      if (!is.null(time_idx)) {
        time_idx <- as.integer(time_idx)
        if (length(time_idx) > 0 && max(time_idx) <= nrow(dense) && min(time_idx) > 0) {
          dense <- dense[time_idx, , drop = FALSE]
        }
        handle$subset_applied$time_idx <- TRUE
      }
    }
  } else {
    # Treat incoming as coefficients (time x component).
    coeff <- as.matrix(incoming)

    if (!is.null(roi_vec)) {
      vox_idx <- which(roi_vec)
      if (identical(storage_order, "component_x_voxel")) {
        basis <- basis[, vox_idx, drop = FALSE]
      } else {
        basis <- basis[vox_idx, , drop = FALSE]
      }
      v_voxels_effective <- if (identical(storage_order, "component_x_voxel")) ncol(basis) else nrow(basis)
      handle$subset_applied$roi_mask <- TRUE
    }
    if (!isTRUE(applied$time_idx)) {
      time_idx <- subset$time_idx %||% subset$time
      if (!is.null(time_idx)) {
        time_idx <- as.integer(time_idx)
        if (length(time_idx) > 0 && max(time_idx) <= nrow(coeff) && min(time_idx) > 0) {
          coeff <- coeff[time_idx, , drop = FALSE]
        }
        handle$subset_applied$time_idx <- TRUE
      }
    }

    if (ncol(coeff) != k_components) {
      abort_lna(
        sprintf(
          "basis invert expects %d components but incoming matrix has %d columns",
          as.integer(k_components), as.integer(ncol(coeff))
        ),
        .subclass = "lna_error_runtime",
        location = "invert_step.basis:dim_mismatch"
      )
    }

    if (identical(storage_order, "component_x_voxel")) {
      dense <- coeff %*% basis
    } else {
      dense <- coeff %*% t(basis)
    }
  }

  handle$update_stash(keys = incoming_key,
                      new_values = setNames(list(dense), target_key))
}

#' Basis Transform - Forward Step
#'
#' Computes a PCA basis matrix from the input data and registers the
#' datasets in the write plan. Only the "pca" method is currently
#' implemented; specifying any other `method` results in an error.
#' @method forward_step basis
#' @export
#' @keywords internal
forward_step.basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  method <- p$method %||% "pca"
  if (!identical(method, "pca")) {
    abort_lna(
      sprintf("method '%s' is not supported; only 'pca' is implemented", method),
      .subclass = "lna_error_validation",
      location = "forward_step.basis:method"
    )
  }
  k <- p$k %||% 20
  center <- p$center %||% TRUE
  scale <- p$scale %||% FALSE
  storage_order <- p$storage_order %||% "component_x_voxel"
  allowed_orders <- c("component_x_voxel", "voxel_x_component")
  if (!storage_order %in% allowed_orders) {
    abort_lna(
      sprintf(
        "Invalid storage_order '%s'. Allowed values: %s",
        storage_order,
        paste(allowed_orders, collapse = ", ")
      ),
      .subclass = "lna_error_validation",
      location = "forward_step.basis:storage_order"
    )
  }

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  X <- handle$get_inputs(input_key)[[1]]
  X <- as_dense_mat(X)

  if (!is.numeric(X)) {
    abort_lna("basis transform requires numeric input matrix",
              .subclass = "lna_error_validation",
              location = "forward_step.basis:input")
  }

  original_k <- k
  dims <- dim(X)
  if (any(dims == 0)) {
    abort_lna(
      "Input matrix for PCA has zero dimensions.",
      .subclass = "lna_error_validation",
      location = "forward_step.basis:input"
    )
  }
  min_dim <- min(dims)


  if (requireNamespace("irlba", quietly = TRUE)) {
    k_max_allowed <- max(1, min_dim - 1)
    if (k > k_max_allowed) {
      warning(sprintf(
        "Requested k=%d but irlba::prcomp_irlba can only compute %d components for %dx%d data; truncating k to %d.",
        original_k, k_max_allowed, nrow(X), ncol(X), k_max_allowed
      ), call. = FALSE)
      k <- k_max_allowed
    }
    k <- max(1, min(k, k_max_allowed))
    fit <- irlba::prcomp_irlba(X, n = k, center = center, scale. = scale)
  } else {
    k_max_allowed <- max(1, min_dim)
    if (k > k_max_allowed) {
      warning(sprintf(
        "Requested k=%d but stats::prcomp can only compute %d components for %dx%d data; truncating k to %d.",
        original_k, k_max_allowed, nrow(X), ncol(X), k_max_allowed
      ), call. = FALSE)
      k <- k_max_allowed
    }
    k <- max(1, min(k, k_max_allowed))
    fit <- stats::prcomp(X, rank. = k, center = center, scale. = scale)
  }

  k_effective <- ncol(fit$rotation)
  # This secondary warning might occur if the data is rank deficient 
  # and the PCA method returns fewer components than k_to_use.
  if (k_effective < k && original_k >= k_effective ) { # Check against original_k for the warning message
    warning(sprintf(
      "PCA fit returned %d components when k was set to %d (original request: k=%d); using %d components.",
      k_effective, k, original_k, k_effective
    ), call. = FALSE)
  }
  # p$k should store the number of components actually in 'rotation'
  p$k <- k_effective 
  rotation <- fit$rotation[, seq_len(k_effective), drop = FALSE]
  # The original p$k was updated above to k_effective, 
  # so params_json will store k_effective.
  mean_vec <- if (isTRUE(center)) fit$center else NULL
  scale_vec <- if (isTRUE(scale)) fit$scale else NULL

  basis_mat <- if (identical(storage_order, "component_x_voxel"))
    t(rotation) else rotation

  plan <- handle$plan
  step_index <- plan$next_index
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  matrix_path <- paste0("/basis/", base_name, "/matrix")
  center_path <- paste0("/basis/", base_name, "/center")
  scale_path <- paste0("/basis/", base_name, "/scale")
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))
  desc$params <- p

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  
  output_key_for_stash <- desc$outputs[[1]]

  datasets <- list(list(path = matrix_path, role = "basis_matrix"))
  if (!is.null(mean_vec)) datasets[[length(datasets) + 1]] <- list(path = center_path, role = "center")
  if (!is.null(scale_vec)) datasets[[length(datasets) + 1]] <- list(path = scale_path, role = "scale")
  desc$datasets <- datasets

  plan$add_descriptor(fname, desc)
  plan$add_payload(matrix_path, basis_mat)
  plan$add_dataset_def(matrix_path, "basis_matrix", type,
                       plan$origin_label, as.integer(step_index),
                       params_json, matrix_path, "eager")

  if (!is.null(mean_vec)) {
    plan$add_payload(center_path, mean_vec)
    plan$add_dataset_def(center_path, "center", type,
                         plan$origin_label, as.integer(step_index),
                         params_json, center_path, "eager")
  }
  if (!is.null(scale_vec)) {
    plan$add_payload(scale_path, scale_vec)
    plan$add_dataset_def(scale_path, "scale", type,
                         plan$origin_label, as.integer(step_index),
                         params_json, scale_path, "eager")
  }
  
  handle$plan <- plan
  output_key_for_stash <- desc$outputs[[1]]
  
  # Stash the original input X, not the PCA scores (fit$x).
  # embed transform will take this X and project it using the stored basis.
  handle <- handle$update_stash(keys = input_key, new_values = setNames(list(X), output_key_for_stash))
  return(handle) # Explicitly return the updated handle
}
