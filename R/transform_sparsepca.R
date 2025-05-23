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
forward_step.myorg.sparsepca <- function(type, desc, handle) {
  p <- desc$params %||% list()
  k <- p$k %||% 2
  alpha <- p$alpha %||% 0.001
  whiten <- p$whiten %||% FALSE
  seed <- p$seed

  inp <- handle$pull_first(c("aggregated_matrix", "dense_mat", "input"))
  input_key <- inp$key
  X <- inp$value

  orig_dims <- dim(X)
  num_dims <- length(orig_dims)

  if (getOption("neuroarchive.debug", FALSE)) {
    cat(sprintf("DEBUG forward_step.myorg.sparsepca: Input X original dims: %s\n", paste(orig_dims, collapse="x")))
  }

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
      location = "forward_step.myorg.sparsepca:input_reshape"
    )
  }

  inp <- handle$pull_first(c("aggregated_matrix", "dense_mat", "input"))
  input_key <- inp$key
  X <- as_dense_mat(inp$value)
  if (isTRUE(whiten)) {
    X_for_pca <- scale(X_for_pca, center = TRUE, scale = TRUE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (requireNamespace("sparsepca", quietly = TRUE)) {
    backend <- "sparsepca"
    fit <- sparsepca::spca(X_for_pca, k = as.integer(k), alpha = alpha,
                           verbose = FALSE)
    basis_V <- fit$loadings # VxK
    embed <- fit$scores # Raw scores

    # Ensure 'd' (from sparsepca fit) is correctly handled and has length 'k'
    if (!is.null(fit$d) && length(fit$d) > 0) {
      d_temp <- fit$d
      if (length(d_temp) > k) {
        d <- d_temp[1:k]
      } else if (length(d_temp) < k) {
        d <- c(d_temp, rep(0, k - length(d_temp))) # Pad with 0 if shorter
      } else {
        d <- d_temp
      }
    } else {
      # If fit$d is NULL or empty, create a dummy 'd' of length k for consistency in storage
      d <- rep(0.0, k) # Using 0.0 as neutral, indicating no scaling factor from sparsepca itself for this path
      if (!is.null(handle$logger)) {
        handle$logger$warn(sprintf("sparsepca backend: fit$d was NULL or empty for k=%d. Using d=0 vector of length k.", k))
      }
    }
  } else if (requireNamespace("irlba", quietly = TRUE)) {
    backend <- "irlba"
    fit <- irlba::irlba(X_for_pca, nv = as.integer(k))
    basis_V <- fit$v # VxK
    embed <- fit$u # U matrix
    d <- fit$d     # Singular values
  } else {
    backend <- "svd"
    sv <- svd(X_for_pca, nu = as.integer(k), nv = as.integer(k))
    basis_V <- sv$v # VxK
    embed <- sv$u # U matrix
    d <- sv$d[seq_len(k)] # Singular values
  }

  # Basis V is (Voxels x Components). We want to store/use KxV for (TxK) %*% (KxV) reconstruction.
  basis_to_store <- Matrix::t(basis_V) # KxV

  if (!is.null(handle$logger)) {
    handle$logger$info(sprintf("sparsepca backend: %s", backend))
    cat(sprintf("DEBUG forward_step: backend=%s, embed dims before user scaling: %s, d length: %d\\n", 
                backend, paste(dim(embed), collapse="x"), length(d)))
    cat(sprintf("DEBUG forward_step: embed dims after all scaling logic: %s\\n", paste(dim(embed), collapse="x")))
    cat(sprintf("DEBUG forward_step: basis_V (loadings) dims: %s; basis_to_store (t(loadings)) dims: %s\\n", 
                paste(dim(basis_V), collapse="x"), paste(dim(basis_to_store), collapse="x")))
  }

  # Scale embeddings for SVD-based methods (irlba, svd) to get coefficients.
  # sparsepca 'embed' (fit$scores) are typically used directly as coefficients.
  if (backend %in% c("irlba", "svd")) {
    if (!is.null(embed) && is.matrix(embed) && !is.null(d) && length(d) > 0 && ncol(embed) == length(d)) {
      embed <- embed %*% diag(d, nrow = length(d), ncol = length(d))
      if (getOption("neuroarchive.debug", FALSE)) {
        cat(sprintf("DEBUG forward_step: Applied SVD-like scaling to embed. New embed dims: %s\\n", paste(dim(embed), collapse="x")))
      }
    } else {
      if (!is.null(handle$logger)) {
        handle$logger$warn(sprintf(
          "Backend %s: Dimension mismatch, NULL d/embed, or d is empty for SVD-like scaling. embed ncol=%s, d length=%s. Skipping embedding scaling.",
          backend, ncol(embed) %||% "NULL", length(d) %||% "NULL"))
      }
    }
  }

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  basis_path <- paste0("/basis/", base, "/basis")
  d_path <- paste0("/basis/", base, "/singular_values")
  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  embed_path <- paste0("/scans/", run_id, "/", base, "/embedding")
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("sparsepca_basis", "sparsepca_embedding")
  desc$datasets <- list(list(path = basis_path, role = "basis_matrix"),
                        list(path = embed_path, role = "coefficients"),
                        list(path = d_path, role = "singular_values"))
  desc$params <- p
  plan$add_descriptor(fname, desc)

  plan$add_payload(basis_path, basis_to_store)
  plan$add_dataset_def(basis_path, "basis_matrix", as.character(type),
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, basis_path, "eager", dtype = NA_character_)
  plan$add_payload(embed_path, embed)
  plan$add_dataset_def(embed_path, "coefficients", as.character(type),
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, embed_path, "eager", dtype = NA_character_)
  plan$add_payload(d_path, d)
  plan$add_dataset_def(d_path, "singular_values", as.character(type),
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, d_path, "eager", dtype = NA_character_)

  handle$plan <- plan
  out_stash <- list()
  out_stash[[desc$outputs[1]]] <- basis_to_store
  out_stash[[desc$outputs[2]]] <- embed
  handle <- handle$update_stash(keys = c(input_key), new_values = out_stash)
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
invert_step.myorg.sparsepca <- function(type, desc, handle) {
  ds <- desc$datasets
  basis_path <- ds[[which(vapply(ds, function(d) d$role, character(1)) == "basis_matrix")]]$path
  embed_path <- ds[[which(vapply(ds, function(d) d$role, character(1)) == "coefficients")]]$path
  d_idx <- which(vapply(ds, function(d) d$role, character(1)) == "singular_values")
  d_path <- if (length(d_idx) == 1) ds[[d_idx]]$path else NULL

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path) # Expected to be KxV ("components x voxels")
  embed <- h5_read(root, embed_path) # Expected to be TxK ("time x components")
  d <- if (!is.null(d_path)) h5_read(root, d_path) else NULL
  
  if (getOption("neuroarchive.debug", FALSE)) {
    cat(sprintf("DEBUG invert_step: embed T_x_K dims = %s, basis K_x_V dims = %s, d length = %s\\n", 
                paste(dim(embed), collapse="x"), 
                paste(dim(basis), collapse="x"),
                length(d)))
  }
  
  # p <- desc$params %||% list() # Params like storage_order removed
  # storage_order related logic removed.

  # Basis is KxV, Embed is TxK. Xhat = Embed %*% Basis => (TxK %*% KxV) = TxV
  if (!is.null(embed) && !is.null(basis) && is.matrix(embed) && is.matrix(basis)) {
    if (ncol(embed) != nrow(basis)) {
      stop(sprintf(
        "Matrix dimensions incompatible for reconstruction: embed (%s) %d_cols vs basis (%s) %d_rows. Expected K_embed == K_basis.",
        paste(dim(embed), collapse = "x"), ncol(embed),
        paste(dim(basis), collapse = "x"), nrow(basis)
      ))
    }
  } else {
    stop("Embed or Basis is NULL or not a matrix in invert_step.")
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
  handle$update_stash(keys = character(),
                      new_values = setNames(list(Xhat), output_key))
}

#' Default parameters for myorg.sparsepca
#' @export
#' @keywords internal
lna_default.myorg.sparsepca <- function() {
  list(k = 50L, alpha = 1e-3, whiten = FALSE, seed = 42L)
}
