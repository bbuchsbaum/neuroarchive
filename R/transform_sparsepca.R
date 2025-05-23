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

  # Use pipeline-provided input key or fallback to common keys
  input_key <- if (!is.null(desc$inputs) && length(desc$inputs) > 0) {
    desc$inputs[[1]]
  } else {
    # Fallback for standalone usage
    inp <- handle$pull_first(c("aggregated_matrix", "dense_mat", "input"))
    inp$key
  }
  
  # Get the data using the determined input key
  X <- if (input_key %in% names(handle$stash)) {
    handle$stash[[input_key]]
  } else {
    abort_lna(
      sprintf("Required input key '%s' not found in stash", input_key),
      .subclass = "lna_error_contract",
      location = "forward_step.myorg.sparsepca:input"
    )
  }
  
  X <- as_dense_mat(X)

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
    # Also, ensure basis_V and embed have k columns, padding if necessary
    # This is mainly for sparsepca backend when k_effective < k_requested.
    k_requested <- as.integer(k)

    if (!is.null(fit$d) && length(fit$d) > 0) {
      d_temp <- fit$d
      if (length(d_temp) > k_requested) {
        d <- d_temp[1:k_requested]
      } else if (length(d_temp) < k_requested) {
        d <- c(d_temp, rep(0, k_requested - length(d_temp))) # Pad with 0 if shorter
      } else {
        d <- d_temp
      }
    } else {
      d <- rep(0.0, k_requested)
      if (!is.null(handle$logger) && backend == "sparsepca") { # Warning specific to sparsepca behavior
          handle$logger$warn(sprintf("sparsepca backend: fit$d was NULL or empty for k=%d. Using d=0 vector of length k.", k_requested))
      }
    }
    
    # For sparsepca, if effective k (ncol of loadings/scores) is less than requested k, pad them.
    if (backend == "sparsepca") {
      k_effective_basis <- if (is.matrix(basis_V)) ncol(basis_V) else 0
      if (k_effective_basis < k_requested) {
        if (getOption("neuroarchive.debug", TRUE)) {
          cat(sprintf("DEBUG sparsepca_forward: Padding basis_V for sparsepca. k_eff=%d, k_req=%d\n", k_effective_basis, k_requested))
        }
        padding_cols_basis <- k_requested - k_effective_basis
        if (k_effective_basis > 0) { # If basis_V is not empty
          basis_V <- cbind(basis_V, matrix(0, nrow = nrow(basis_V), ncol = padding_cols_basis))
        } else { # basis_V might be Vx0 or NULL if input V=0 or PCA failed badly
          num_voxels <- if (is.matrix(X_for_pca)) ncol(X_for_pca) else 0
          basis_V <- matrix(0, nrow = num_voxels, ncol = k_requested)
        }
      }
      
      k_effective_embed <- if (is.matrix(embed)) ncol(embed) else 0
      if (k_effective_embed < k_requested) {
        if (getOption("neuroarchive.debug", TRUE)) {
          cat(sprintf("DEBUG sparsepca_forward: Padding embed for sparsepca. k_eff=%d, k_req=%d\n", k_effective_embed, k_requested))
        }
        padding_cols_embed <- k_requested - k_effective_embed
        if (k_effective_embed > 0) { # If embed is not empty
          embed <- cbind(embed, matrix(0, nrow = nrow(embed), ncol = padding_cols_embed))
        } else { # embed might be Tx0 or NULL
          num_time <- if (is.matrix(X_for_pca)) nrow(X_for_pca) else 0
          embed <- matrix(0, nrow = num_time, ncol = k_requested)
        }
      }
    }

    # For irlba and svd, k_eff is already min(k_requested, rank). d, basis_V, embed are already consistent.
    # However, d might be shorter than k_requested if rank < k_requested. Pad d.
    if (backend %in% c("irlba", "svd")) {
      if (length(d) < k_requested) {
          if (getOption("neuroarchive.debug", TRUE)) {
              cat(sprintf("DEBUG sparsepca_forward: Padding d for %s. len_d=%d, k_req=%d\n", backend, length(d), k_requested))
          }
          d <- c(d, rep(0, k_requested - length(d)))
      }
      # Ensure basis_V and embed also have k_requested columns if svd/irlba returned fewer due to rank < k_requested
      # basis_V is V x k_eff, embed is T x k_eff
      k_effective_svd <- length(fit$d) # True number of components found by svd/irlba
      if (k_effective_svd < k_requested) {
          if (is.matrix(basis_V) && ncol(basis_V) == k_effective_svd) {
               if (getOption("neuroarchive.debug", TRUE)) {
                  cat(sprintf("DEBUG sparsepca_forward: Padding basis_V for %s. k_eff=%d, k_req=%d\n", backend, k_effective_svd, k_requested))
              }
              basis_V <- cbind(basis_V, matrix(0, nrow = nrow(basis_V), ncol = k_requested - k_effective_svd))
          }
          if (is.matrix(embed) && ncol(embed) == k_effective_svd) {
              if (getOption("neuroarchive.debug", TRUE)) {
                  cat(sprintf("DEBUG sparsepca_forward: Padding embed for %s. k_eff=%d, k_req=%d\n", backend, k_effective_svd, k_requested))
              }
              embed <- cbind(embed, matrix(0, nrow = nrow(embed), ncol = k_requested - k_effective_svd))
          }
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
  basis_to_store <- Matrix::t(basis_V) # KxV (should be k_requested x V)

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
  
  # Use provided outputs or default names
  if (is.null(desc$outputs) || length(desc$outputs) == 0) {
    desc$outputs <- c("sparsepca_basis", "sparsepca_embedding")
    basis_key <- "sparsepca_basis"
    embed_key <- "sparsepca_embedding"
  } else {
    # Use the provided output key(s) from the pipeline
    if (length(desc$outputs) == 1) {
      # If only one output key is provided, use it for the embedding (main output)
      # Keep basis with a secondary name for internal use
      embed_key <- desc$outputs[[1]]
      basis_key <- "sparsepca_basis"
      desc$outputs <- c(basis_key, embed_key)
    } else {
      # If two or more outputs are provided, use them as-is
      basis_key <- desc$outputs[[1]]
      embed_key <- desc$outputs[[2]]
    }
  }
  
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
  out_stash[[basis_key]] <- basis_to_store
  out_stash[[embed_key]] <- embed
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
  basis_raw <- h5_read(root, basis_path)
  embed_raw <- h5_read(root, embed_path)
  d <- if (!is.null(d_path)) h5_read(root, d_path) else NULL

  # Get K (number of components) from descriptor to help reshape if needed
  # Fallback for k is based on the forward step's default if not in params.
  # The forward step ensures p$k is populated.
  k_components <- desc$params$k %||% lna_default.myorg.sparsepca()$k
  if (is.null(k_components) || !is.numeric(k_components) || k_components < 0) {
    # This should ideally not happen if descriptor is well-formed
    k_components <- 0 # A safe default if k is missing or invalid
    warning("k_components not found or invalid in sparsepca descriptor, defaulting to 0.")
  }

  # Attempt to reshape embed if it's a vector
  embed <- embed_raw
  if (!is.matrix(embed_raw) && is.vector(embed_raw) && !is.list(embed_raw)) {
    if (k_components > 0 && length(embed_raw) %% k_components == 0) {
      n_time_embed <- length(embed_raw) / k_components
      embed <- matrix(embed_raw, nrow = n_time_embed, ncol = k_components)
      message(sprintf("[DEBUG sparsepca_invert] Reshaped embed from vector to %dx%d matrix", n_time_embed, k_components))
    } else if (k_components == 0 && length(embed_raw) > 0) { # T x 0 case, length is n_time
      embed <- matrix(numeric(0), nrow = length(embed_raw), ncol = 0)
      message(sprintf("[DEBUG sparsepca_invert] Reshaped embed to %dx0 matrix from vector of length %d", length(embed_raw), length(embed_raw)))
    } else if (k_components > 0 && length(embed_raw) > 0 && length(embed_raw) %% k_components !=0) {
        # Potentially 1xK case if T=1, or Kx1 if V=1 for basis
        # For embed (TxK), if T=1, length is K. Reshape to 1xK.
        if (length(embed_raw) == k_components) {
            embed <- matrix(embed_raw, nrow = 1, ncol = k_components)
            message(sprintf("[DEBUG sparsepca_invert] Reshaped embed from vector to 1x%d matrix", k_components))
        }
    }
  }

  # Attempt to reshape basis if it's a vector
  basis <- basis_raw
  if (!is.matrix(basis_raw) && is.vector(basis_raw) && !is.list(basis_raw)) {
    if (k_components > 0 && length(basis_raw) %% k_components == 0) {
      n_voxels_basis <- length(basis_raw) / k_components
      basis <- matrix(basis_raw, nrow = k_components, ncol = n_voxels_basis) # Basis is KxV
      message(sprintf("[DEBUG sparsepca_invert] Reshaped basis from vector to %dx%d matrix", k_components, n_voxels_basis))
    } else if (k_components == 0 && length(basis_raw) > 0) { # 0 x V case, length is n_voxels
      basis <- matrix(numeric(0), nrow = 0, ncol = length(basis_raw))
       message(sprintf("[DEBUG sparsepca_invert] Reshaped basis to 0x%d matrix from vector of length %d", length(basis_raw), length(basis_raw)))
    } else if (k_components > 0 && length(basis_raw) > 0 && length(basis_raw) %% k_components != 0) {
        # For basis (KxV), if V=1, length is K. Reshape to Kx1.
        if (length(basis_raw) == k_components) {
            basis <- matrix(basis_raw, nrow = k_components, ncol = 1)
            message(sprintf("[DEBUG sparsepca_invert] Reshaped basis from vector to %dx1 matrix", k_components))
        }
    }
  }
  
  message(sprintf("[DEBUG sparsepca_invert] AFTER RESHAPE ATTEMPT:"))
  message(sprintf("[DEBUG sparsepca_invert] embed - class: %s, is.null: %s, is.matrix: %s, length: %s, dim: %s", 
                  class(embed)[1], is.null(embed), is.matrix(embed), length(embed), paste(dim(embed), collapse="x")))
  message(sprintf("[DEBUG sparsepca_invert] basis - class: %s, is.null: %s, is.matrix: %s, length: %s, dim: %s", 
                  class(basis)[1], is.null(basis), is.matrix(basis), length(basis), paste(dim(basis), collapse="x")))

  if (getOption("neuroarchive.debug", FALSE)) {
    cat(sprintf("DEBUG invert_step: embed T_x_K dims = %s, basis K_x_V dims = %s, d length = %s\\n", 
                paste(dim(embed), collapse="x"), 
                paste(dim(basis), collapse="x"),
                length(d)))
  }

  # Basis is KxV, Embed is TxK. Xhat = Embed %*% Basis => (TxK %*% KxV) = TxV
  if (!is.null(embed) && !is.null(basis) && is.matrix(embed) && is.matrix(basis)) {
    if (ncol(embed) != nrow(basis)) {
      stop(sprintf(
        "Matrix dimensions incompatible for reconstruction: embed (%s) %d_cols vs basis (%s) %d_rows. Expected K_embed == K_basis.",
        paste(dim(embed), collapse = "x"), ncol(embed),
        paste(dim(basis), collapse = "x"), nrow(basis)
      ))
    }
  } else if (!is.matrix(embed) || !is.matrix(basis)) {
    # Handle case where embed/basis are stored dimension vectors from empty arrays
    if (length(embed) == 2 && length(basis) == 2 && 
        all(embed >= 0) && all(basis >= 0) && 
        all(embed == as.integer(embed)) && all(basis == as.integer(basis))) {
      # These look like stored dimensions - reconstruct the original empty matrices
      embed_dims <- as.integer(embed)
      basis_dims <- as.integer(basis)
      
      # Reconstruct the matrices
      embed <- array(numeric(0), dim = embed_dims)
      basis <- array(numeric(0), dim = basis_dims)
      
      # Check if we can do matrix multiplication
      if (embed_dims[2] != basis_dims[1]) {
        # Dimensions don't match for multiplication - create empty result
        Xhat <- matrix(numeric(0), nrow = embed_dims[1], ncol = basis_dims[2])
      } else {
        # Dimensions match - do the multiplication (which will result in empty matrix)
        Xhat <- embed %*% basis
      }
    } else {
      stop("Embed or Basis is NULL or not a matrix in invert_step.")
    }
  } else {
    stop("Embed or Basis is NULL or not a matrix in invert_step.")
  }

  if (!exists("Xhat")) {
    Xhat <- embed %*% basis
  }

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
