#' Basis Transform - Inverse Step
#'
#' Reconstructs data from coefficients using a stored basis matrix.
#'
#' The `basis` dataset may be stored either as a component-by-voxel matrix
#' (`storage_order = "component_x_voxel"`) or as a voxel-by-component matrix
#' (`storage_order = "voxel_x_component"`). After optional transposition, it
#' should have dimensions `n_voxel x n_component` for reconstruction.
#' @keywords internal
invert_step.basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  storage_order <- p$storage_order %||% "component_x_voxel"

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

  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  input_key  <- desc$inputs[[1]] %||% "dense_mat"

  if (!handle$has_key(coeff_key)) {
    # Nothing to reconstruct; return handle unchanged
    return(handle)
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)

  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    vox_idx <- which(as.logical(subset$roi_mask))
    if (identical(storage_order, "component_x_voxel")) {
      basis <- basis[, vox_idx, drop = FALSE]
    } else {
      basis <- basis[vox_idx, , drop = FALSE]
    }
  }
  if (!is.null(subset$time_idx)) {
    coeff <- coeff[subset$time_idx, , drop = FALSE]
  }

  if (identical(storage_order, "component_x_voxel")) {
    basis <- t(basis)
  }

  if (nrow(basis) == ncol(coeff)) {
    dense <- coeff %*% basis
  } else {
    dense <- coeff %*% t(basis)
  }

  handle$update_stash(keys = coeff_key, new_values = setNames(list(dense), input_key))
}

#' Basis Transform - Forward Step
#'
#' Computes a basis matrix (e.g., via PCA) from the input data and
#' registers the datasets in the write plan.
#' @keywords internal
forward_step.basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  method <- p$method %||% "pca"
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
  min_dim <- min(dim(X))

  if (requireNamespace("irlba", quietly = TRUE)) {
    # For irlba, k must be < min_dim. Max k is min_dim - 1.
    k_max_allowed <- max(1, min_dim - 1) 
    if (k > k_max_allowed) {
      warning(sprintf(
        "Requested k=%d but irlba::prcomp_irlba can only compute %d components for %dx%d data; truncating k to %d.",
        original_k, k_max_allowed, nrow(X), ncol(X), k_max_allowed
      ), call. = FALSE)
      k <- k_max_allowed
    }
    # Ensure k is not zero if min_dim was 1 (k_max_allowed would be 0)
    k <- max(1, k) 
    if (k >= min_dim) { # Final safety check if min_dim is 1 or 2
        # This case should ideally not be hit if k_max_allowed is correct
        # but prcomp_irlba errors if k is not strictly less than min_dim
        k <- max(1, min_dim -1)
         if (k == 0 && min_dim > 0) k <- 1 # if min_dim = 1, k becomes 0, fix to 1
         if (k==0 && min_dim ==0) { # edge case of 0-dim input.
            abort_lna("Input matrix for PCA has zero dimensions.", 
                        .subclass="lna_error_validation", 
                        location="forward_step.basis:input")
         }
    }

    fit <- irlba::prcomp_irlba(X, n = k, center = center, scale. = scale)
  } else {
    # For stats::prcomp, rank. can be <= min_dim. Max k is min_dim.
    k_max_allowed <- max(1, min_dim)
    if (k > k_max_allowed) {
      warning(sprintf(
        "Requested k=%d but stats::prcomp can only compute %d components for %dx%d data; truncating k to %d.",
        original_k, k_max_allowed, nrow(X), ncol(X), k_max_allowed
      ), call. = FALSE)
      k <- k_max_allowed
    }
    # Ensure k is not zero
    k <- max(1, k)
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
  desc$outputs <- character()
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
  # keep input in the stash for subsequent transforms (e.g., 'embed')
  handle
}
