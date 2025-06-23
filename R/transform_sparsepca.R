#' Sparse PCA Transform - Forward Step
#'
#' Performs a sparse PCA on the input matrix. If the optional `sparsepca`
#' package is available, the transform uses `sparsepca::spca()` to compute
#' sparse loadings. Otherwise it falls back to a truncated SVD via `irlba`
#' (or base `svd`). Columns may be optionally whitened prior to fitting.
#' The chosen backend and singular values are recorded for later use.
#' This example demonstrates how an external plugin transform
#' might integrate with the LNA pipeline.
#'
#' @param storage_order Character string specifying the orientation of the
#'   basis matrix. Either "component_x_voxel" (default) or
#'   "voxel_x_component".
#' @keywords internal
forward_step.sparsepca <- function(type, desc, handle) {
  p <- desc$params %||% list()
  k <- p$k %||% 2 # Default k from original mock
  alpha <- p$alpha %||% 0.001
  whiten <- p$whiten %||% FALSE
  seed <- p$seed

  input_key <- if (!is.null(desc$inputs) && length(desc$inputs) > 0) {
    desc$inputs[[1]]
  } else {
    inp <- handle$pull_first(c("aggregated_matrix", "dense_mat", "input"))
    inp$key
  }
  
  X <- if (input_key %in% names(handle$stash)) {
    handle$stash[[input_key]]
  } else {
    abort_lna(
      sprintf("Required input key '%s' not found in stash", input_key),
      .subclass = "lna_error_contract",
      location = "forward_step.sparsepca:input"
    )
  }
  
  X <- as_dense_mat(X)

  orig_dims <- dim(X)
  num_dims <- length(orig_dims)

  if (num_dims == 2) {
    X_for_pca <- X
  } else if (num_dims == 3 && orig_dims[num_dims] == 1) {
    X_for_pca <- matrix(X, nrow = orig_dims[1], ncol = orig_dims[2])
  } else if (num_dims == 4) {
    n_time <- orig_dims[4]
    n_voxels <- prod(orig_dims[1:3])
    X_for_pca <- matrix(aperm(X, c(4, 1, 2, 3)), nrow = n_time, ncol = n_voxels)
  } else {
    err_msg <- paste0("Input data X for sparsepca has unexpected dimensions: ", paste(orig_dims, collapse="x"),
                      ". Expected a 2D (Time x Voxels), 3D (Time x Voxels x 1), or 4D (Spatial x Spatial x Spatial x Time) array.")
    abort_lna(
      err_msg,
      .subclass = "lna_error_validation",
      location = "forward_step.sparsepca:input_reshape"
    )
  }

  if (isTRUE(whiten)) {
    X_for_pca <- scale(X_for_pca, center = TRUE, scale = TRUE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  k_requested <- as.integer(p$n_components %||% p$k %||% 10L) 
  k_eff <- 0 

  if (k_requested <= 0) {
    basis_V <- matrix(0, nrow = ncol(X_for_pca), ncol = 0)
    embed <- matrix(0, nrow = nrow(X_for_pca), ncol = 0)
    d <- numeric(0)
    k_eff <- 0
    backend <- "k_zero_case"
  } else if (requireNamespace("sparsepca", quietly = TRUE)) {
    backend <- "sparsepca"
    fit <- sparsepca::spca(X_for_pca, k = k_requested, alpha = alpha, verbose = FALSE, center = FALSE, scale = FALSE)
    basis_V <- fit$loadings 
    embed <- X_for_pca %*% fit$loadings 
    d <- fit$d 
    k_eff <- if (is.matrix(basis_V)) ncol(basis_V) else 0
  } else if (requireNamespace("irlba", quietly = TRUE)) {
    backend <- "irlba"
    actual_k_for_irlba <- min(k_requested, min(dim(X_for_pca)) -1)
    if (actual_k_for_irlba < 1) { 
        backend <- "svd_fallback_small_data"
        sv <- svd(X_for_pca, nu = min(nrow(X_for_pca), k_requested) , nv = min(ncol(X_for_pca), k_requested))
        k_eff <- min(length(sv$d), k_requested)
        basis_V <- sv$v[, 1:k_eff, drop = FALSE] 
        embed <- sv$u[, 1:k_eff, drop = FALSE]   
        d <- sv$d[1:k_eff]                       
    } else {
        fit <- irlba::irlba(X_for_pca, nv = actual_k_for_irlba, nu = actual_k_for_irlba) 
        basis_V <- fit$v 
        embed <- fit$u   
        d <- fit$d       
        k_eff <- length(d)
    }
  } else {
    backend <- "svd"
    sv <- svd(X_for_pca, nu = min(nrow(X_for_pca), k_requested) , nv = min(ncol(X_for_pca), k_requested))
    k_eff <- min(length(sv$d), k_requested)
    basis_V <- sv$v[, 1:k_eff, drop = FALSE] 
    embed <- sv$u[, 1:k_eff, drop = FALSE]   
    d <- sv$d[1:k_eff]                       
  }

  if (k_eff < k_requested && k_requested > 0) {
    if (is.matrix(basis_V)) {
      basis_V_padded <- matrix(0, nrow = nrow(basis_V), ncol = k_requested)
      if (k_eff > 0) basis_V_padded[, 1:k_eff] <- basis_V
      basis_V <- basis_V_padded
    } else { 
      basis_V <- matrix(0, nrow = ncol(X_for_pca), ncol = k_requested)
    }

    if (is.matrix(embed)) {
      embed_padded <- matrix(0, nrow = nrow(embed), ncol = k_requested)
      if (k_eff > 0) embed_padded[, 1:k_eff] <- embed
      embed <- embed_padded
    } else {
      embed <- matrix(0, nrow = nrow(X_for_pca), ncol = k_requested)
    }
    
    d_padded <- numeric(k_requested)
    if (k_eff > 0 && length(d) > 0) d_padded[1:k_eff] <- d
    d <- d_padded
  }
  
  if (backend %in% c("irlba", "svd", "svd_fallback_small_data")) {
    if (is.matrix(embed) && ncol(embed) > 0 && length(d) >= ncol(embed)) {
         diag_d_matrix <- diag(d[1:ncol(embed)], nrow = ncol(embed), ncol = ncol(embed))
         embed <- embed %*% diag_d_matrix 
    } else if (ncol(embed) == 0 && k_requested > 0) {
    } else if (ncol(embed) > 0 && length(d) < ncol(embed)){
        if (k_eff > 0 && length(d) >= k_eff) {
            diag_d_matrix <- diag(d[1:k_eff], nrow = k_eff, ncol = k_eff)
            embed <- embed %*% diag_d_matrix
        } else {
        }
    }
  }

  basis_to_store <- Matrix::t(basis_V) 
  
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  
  current_step_idx_for_path <- plan$next_index - 1 
  run_id_for_path <- handle$current_run_id %||% "run-01"
  run_id_for_path <- sanitize_run_id(run_id_for_path) 

  basis_h5_path <- file.path("/basis", sprintf("%s_s%02d_%s_basis", run_id_for_path, current_step_idx_for_path, desc$type))
  embed_h5_path <- file.path("/scans", run_id_for_path, sprintf("s%02d_%s_embedding", current_step_idx_for_path, desc$type))
  d_h5_path     <- file.path("/basis", sprintf("%s_s%02d_%s_d", run_id_for_path, current_step_idx_for_path, desc$type))
  
  if (is.null(desc$params)) desc$params <- list()
  desc$params$k_eff <- k_eff
  desc_params_for_json <- desc$params 
  desc_params_for_json$original_V_dim <- ncol(X_for_pca) 
  desc_params_for_json$backend <- backend
  params_json <- as.character(jsonlite::toJSON(desc_params_for_json, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  
  main_output_key <- desc$outputs[[1]] 
  
  desc$datasets <- list(
    list(path = basis_h5_path, role = "basis_matrix"),
    list(path = embed_h5_path, role = "coefficients"), 
    list(path = d_h5_path, role = "singular_values")
  )
  
  plan$add_descriptor(fname, desc) 
  
  plan$add_payload(basis_h5_path, basis_to_store)
  plan$add_dataset_def(basis_h5_path, "basis_matrix", as.character(type),
                       run_id_for_path, as.integer(current_step_idx_for_path),
                       params_json, basis_h5_path, "eager", dtype = NA_character_)
                       
  plan$add_payload(embed_h5_path, embed) 
  plan$add_dataset_def(embed_h5_path, "coefficients", as.character(type),
                       run_id_for_path, as.integer(current_step_idx_for_path),
                       params_json, embed_h5_path, "eager", dtype = NA_character_)

  plan$add_payload(d_h5_path, d)
  plan$add_dataset_def(d_h5_path, "singular_values", as.character(type),
                       run_id_for_path, as.integer(current_step_idx_for_path),
                       params_json, d_h5_path, "eager", dtype = NA_character_)

  handle$plan <- plan
  
  stash_updates <- setNames(list(embed), main_output_key)
  stash_updates$sparsepca_basis <- basis_to_store 
  stash_updates$sparsepca_d <- d
  
  handle <- handle$update_stash(keys = c(input_key), new_values = stash_updates)
  handle
}

#' Sparse PCA Transform - Inverse Step
#'
#' Reconstructs data from the sparse PCA coefficients and basis matrix.
#'
#' @param storage_order Character string specifying the orientation of the
#'   basis matrix. Either "component_x_voxel" (default) or
#'   "voxel_x_component".
#' @keywords internal
invert_step.sparsepca <- function(type, desc, handle) {
  ds <- desc$datasets
  basis_path <- ds[[which(vapply(ds, function(d) d$role, character(1)) == "basis_matrix")]]$path
  embed_dataset_info <- Filter(function(d) d$role == "coefficients", desc$datasets)
  if (length(embed_dataset_info) == 0) {
      stop("Path for 'coefficients' (embedding) dataset not found in sparsepca descriptor.")
  }
  embed_path <- embed_dataset_info[[1]]$path
  
  d_idx <- which(vapply(ds, function(d) d$role, character(1)) == "singular_values")
  d_path <- if (length(d_idx) == 1) ds[[d_idx]]$path else NULL

  root <- handle$h5[["/"]]
  basis_raw <- h5_read(root, basis_path) 
  embed_raw <- h5_read(root, embed_path) 
  d <- if (!is.null(d_path)) h5_read(root, d_path) else NULL

  k_components <- desc$params$k %||% lna_default.sparsepca()$k
  if (is.null(k_components) || !is.numeric(k_components) || k_components < 0) {
    k_components <- 0 
  }
  
  k_eff_forward <- desc$params$k_eff
  original_V_dim <- desc$params$original_V_dim

  embed <- embed_raw
  if (!is.matrix(embed_raw) && is.vector(embed_raw) && !is.list(embed_raw)) {
    if (k_components > 0 && length(embed_raw) > 0 && length(embed_raw) %% k_components == 0) {
      n_time_embed <- length(embed_raw) / k_components
      embed <- matrix(embed_raw, nrow = n_time_embed, ncol = k_components)
    } else if (k_components == 0 && length(embed_raw) > 0) { 
      embed <- matrix(numeric(0), nrow = length(embed_raw), ncol = 0)
    } else if (k_components > 0 && length(embed_raw) > 0 && length(embed_raw) %% k_components != 0) {
        if (length(embed_raw) == k_components && NROW(embed_raw) == 1) { # T=1 case, ensure NROW is used for vectors
             embed <- matrix(embed_raw, nrow = 1, ncol = k_components)
        } else {
        }
    } else if (length(embed_raw) == 0 && k_components >= 0) { 
        N_time_assumed <- if (k_components > 0) 0 else 0 
        embed <- matrix(numeric(0), nrow=N_time_assumed, ncol=k_components)
    }
  }

  basis <- basis_raw 
  if (!is.matrix(basis_raw) && is.vector(basis_raw) && !is.list(basis_raw)) {
    v_dim_for_reshape <- original_V_dim %||% (if(k_components > 0 && length(basis_raw) > 0 && length(basis_raw) %% k_components == 0) length(basis_raw)/k_components else 0)
    if (k_components > 0 && length(basis_raw) > 0 && v_dim_for_reshape > 0 && length(basis_raw) == k_components * v_dim_for_reshape) {
      basis <- matrix(basis_raw, nrow = k_components, ncol = v_dim_for_reshape) 
    } else if (k_components == 0 && length(basis_raw) > 0) { 
      basis <- matrix(numeric(0), nrow = 0, ncol = length(basis_raw))
    } else if (k_components > 0 && length(basis_raw) > 0 && (v_dim_for_reshape == 0 || length(basis_raw) != k_components * v_dim_for_reshape )) {
        if (length(basis_raw) == k_components && NCOL(basis_raw) == 1) { # V=1 case, ensure NCOL is used for vectors
            basis <- matrix(basis_raw, nrow = k_components, ncol = 1)
        } else {
        }
    } else if (length(basis_raw) == 0 && k_components >= 0) {
        N_vox_assumed <- original_V_dim %||% (if (k_components > 0) 0 else 0)
        basis <- matrix(numeric(0), nrow=k_components, ncol=N_vox_assumed)
    }
  }
  
  if (!is.matrix(embed) || !is.matrix(basis)) {
    stop("Embed or Basis could not be resolved to a matrix in invert_step.sparsepca.")
  }
  
  if (ncol(embed) != nrow(basis)) {
    stop(sprintf(
      "Matrix dimensions incompatible for reconstruction: embed (%s) %d_cols vs basis (%s) %d_rows. Expected K_embed == K_basis.",
      paste(dim(embed), collapse = "x"), ncol(embed),
      paste(dim(basis), collapse = "x"), nrow(basis)
    ))
  }

  Xhat <- embed %*% basis

  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    vox_idx <- which(as.logical(subset$roi_mask))
    Xhat <- Xhat[, vox_idx, drop = FALSE]
  }
  if (!is.null(subset$time_idx)) {
    Xhat <- Xhat[subset$time_idx, , drop = FALSE]
  }

  output_key <- desc$inputs[[1]] %||% "input" 
  handle$update_stash(keys = names(handle$stash), 
                      new_values = setNames(list(Xhat), output_key))
}

#' Default parameters for myorg.sparsepca
#' @export
#' @keywords internal
lna_default.sparsepca <- function() {
  list(k = 50L, alpha = 1e-3, whiten = FALSE, seed = 42L, n_components=50L)
}
