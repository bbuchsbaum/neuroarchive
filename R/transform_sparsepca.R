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
  storage_order <- p$storage_order %||% "component_x_voxel"
  whiten <- p$whiten %||% FALSE
  seed <- p$seed

  allowed_orders <- c("component_x_voxel", "voxel_x_component")
  if (!storage_order %in% allowed_orders) {
    abort_lna(
      sprintf(
        "Invalid storage_order '%s'. Allowed values: %s",
        storage_order,
        paste(allowed_orders, collapse = ", ")
      ),
      .subclass = "lna_error_validation",
      location = "forward_step.myorg.sparsepca:storage_order"
    )
  }

  inp <- handle$pull_first(c("aggregated_matrix", "dense_mat", "input"))
  input_key <- inp$key
  X <- as_dense_mat(inp$value)
  if (isTRUE(whiten)) {
    X <- scale(X, center = TRUE, scale = TRUE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (requireNamespace("sparsepca", quietly = TRUE)) {
    backend <- "sparsepca"
    fit <- sparsepca::spca(X, k = as.integer(k), alpha = alpha,
                           verbose = FALSE)
    basis <- fit$loadings
    embed <- fit$scores
    d <- svd(X, nu = 0, nv = 0)$d[seq_len(k)]
  } else if (requireNamespace("irlba", quietly = TRUE)) {
    backend <- "irlba"
    fit <- irlba::irlba(X, nv = as.integer(k))
    basis <- fit$v
    embed <- fit$u
    d <- fit$d
  } else {
    backend <- "svd"
    sv <- svd(X, nu = as.integer(k), nv = as.integer(k))
    basis <- sv$v
    embed <- sv$u
    d <- sv$d[seq_len(k)]
  }
  if (!is.null(handle$logger)) {
    handle$logger$info(sprintf("sparsepca backend: %s", backend))
  }
  embed <- embed %*% diag(d)

  if (identical(storage_order, "voxel_x_component")) {
    basis <- Matrix::t(basis)
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

  plan$add_payload(basis_path, basis)
  plan$add_dataset_def(basis_path, "basis_matrix", as.character(type),
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, basis_path, "eager")
  plan$add_payload(embed_path, embed)
  plan$add_dataset_def(embed_path, "coefficients", as.character(type),
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, embed_path, "eager")
  plan$add_payload(d_path, d)
  plan$add_dataset_def(d_path, "singular_values", as.character(type),
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, d_path, "eager")

  handle$plan <- plan
  # Stash the actual basis and embedding, not just TRUE/FALSE flags
  out_stash <- list()
  out_stash[[desc$outputs[1]]] <- basis # e.g., sparsepca_basis
  out_stash[[desc$outputs[2]]] <- embed # e.g., sparsepca_embedding
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
  basis <- h5_read(root, basis_path)
  embed <- h5_read(root, embed_path)
  d <- if (!is.null(d_path)) h5_read(root, d_path) else NULL
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
      location = "invert_step.myorg.sparsepca:storage_order"
    )
  }
  if (identical(storage_order, "voxel_x_component")) {
    basis <- Matrix::t(basis)
  }

  if (!is.null(d) && ncol(embed) == length(d)) {
    embed <- embed %*% diag(d)
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
  list(k = 50L, alpha = 1e-3, whiten = FALSE, seed = 42L,
       storage_order = "component_x_voxel")
}
