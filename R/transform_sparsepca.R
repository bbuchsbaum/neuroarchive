#' Sparse PCA Transform - Forward Step
#'
#' Performs a sparse PCA on the input matrix. If the optional `sparsepca`
#' package is available, the transform uses `sparsepca::spca()` to compute
#' sparse loadings. Otherwise it falls back to a truncated SVD via `irlba`
#' (or base `svd`). This example demonstrates how an external plugin transform
#' might integrate with the LNA pipeline.
#' @keywords internal
forward_step.myorg.sparsepca <- function(type, desc, handle) {
  p <- desc$params %||% list()
  k <- p$k %||% 2
  alpha <- p$alpha %||% 0.001
  storage_order <- p$storage_order %||% "component_x_voxel"

  input_key <- if (handle$has_key("aggregated_matrix")) "aggregated_matrix" else "dense_mat"
  X <- handle$get_inputs(input_key)[[1]]
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (requireNamespace("sparsepca", quietly = TRUE)) {
    fit <- sparsepca::spca(X, k = as.integer(k), alpha = alpha,
                           verbose = FALSE)
    basis <- fit$loadings
    embed <- fit$scores
  } else if (requireNamespace("irlba", quietly = TRUE)) {
    fit <- irlba::irlba(X, nv = as.integer(k))
    basis <- fit$v
    embed <- fit$u %*% diag(fit$d)
  } else {
    sv <- svd(X, nu = as.integer(k), nv = as.integer(k))
    basis <- sv$v
    embed <- sv$u %*% diag(sv$d[seq_len(k)])
  }

  if (identical(storage_order, "voxel_x_component")) {
    basis <- Matrix::t(basis)
  }

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  basis_path <- paste0("/basis/", base, "/basis")
  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  embed_path <- paste0("/scans/", run_id, "/", base, "/embedding")
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("sparsepca_basis", "sparsepca_embedding")
  desc$datasets <- list(list(path = basis_path, role = "basis_matrix"),
                        list(path = embed_path, role = "coefficients"))
  plan$add_descriptor(fname, desc)

  plan$add_payload(basis_path, basis)
  plan$add_dataset_def(basis_path, "basis_matrix", type,
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, basis_path, "eager")
  plan$add_payload(embed_path, embed)
  plan$add_dataset_def(embed_path, "coefficients", type,
                       handle$plan$origin_label, as.integer(plan$next_index - 1),
                       params_json, embed_path, "eager")

  handle$plan <- plan
  out <- list(sparsepca_basis = TRUE, sparsepca_embedding = TRUE)
  handle <- handle$update_stash(c(input_key), out)
  handle
}

#' Sparse PCA Transform - Inverse Step
#'
#' Reconstructs data from the sparse PCA coefficients and basis matrix.
#' @keywords internal
invert_step.myorg.sparsepca <- function(type, desc, handle) {
  ds <- desc$datasets
  basis_path <- ds[[which(vapply(ds, function(d) d$role, character(1)) == "basis_matrix")]]$path
  embed_path <- ds[[which(vapply(ds, function(d) d$role, character(1)) == "coefficients")]]$path

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  embed <- h5_read(root, embed_path)
  p <- desc$params %||% list()
  storage_order <- p$storage_order %||% "component_x_voxel"
  if (identical(storage_order, "voxel_x_component")) {
    basis <- Matrix::t(basis)
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

  output_key <- desc$inputs[[1]] %||% "dense_mat"
  handle$update_stash(keys = character(),
                      new_values = setNames(list(Xhat), output_key))
}

#' Default parameters for myorg.sparsepca
#' @export
#' @keywords internal
lna_default.myorg.sparsepca <- function() {
  list(k = 50L, alpha = 1e-3, whiten = FALSE, storage_order = "component_x_voxel")
}
