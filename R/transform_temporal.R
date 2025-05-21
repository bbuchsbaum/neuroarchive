#' Temporal Transform - Forward Step
#'
#' Projects data onto a temporal basis. Currently supports
#' Discrete Cosine Transform (DCT) and B-spline bases.
#' @keywords internal
forward_step.temporal <- function(type, desc, handle) {
  p <- desc$params %||% list()
  kind <- p$kind %||% "dct"
  n_basis <- p$n_basis %||% 5
  order <- p$order %||% 4
  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  X <- handle$get_inputs(input_key)[[1]]

  if (is.array(X) && length(dim(X)) > 2) {
    d <- dim(X)
    time_dim <- d[length(d)]
    vox_dim <- prod(d[-length(d)])
    mat <- matrix(as.numeric(aperm(X, c(length(d), seq_len(length(d) - 1)))),
                  nrow = time_dim, ncol = vox_dim)
  } else {
    mat <- as.matrix(X)
    time_dim <- nrow(mat)
    vox_dim <- ncol(mat)
  }

  if (kind == "bspline") {
    deg <- max(1L, order - 1L)
    basis <- splines::bs(seq_len(time_dim), df = n_basis, degree = deg,
                         intercept = TRUE)
    knots <- attr(basis, "knots")
  } else if (kind == "dct") {
    basis <- dct_basis(time_dim, n_basis)
    knots <- numeric()
  } else {
    abort_lna(sprintf("Unsupported temporal kind '%s'", kind),
              .subclass = "lna_error_validation",
              location = "forward_step.temporal:kind")
  }

  coeff <- crossprod(basis, mat)
  run_id <- handle$current_run_id %||% "run-01"
  plan <- handle$plan
  step_index <- plan$next_index
  fname <- plan$get_next_filename("temporal")
  base_name <- tools::file_path_sans_ext(fname)
  basis_path <- paste0("/temporal_basis/", base_name, "/matrix")
  coef_path <- paste0("/scans/", run_id,
                      "/temporal_coeffs/", base_name, "/coefficients")
  knot_path <- paste0("/temporal_basis/", base_name, "/knots")
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("temporal_coefficients")
  desc$datasets <- list(
    list(path = basis_path, role = "temporal_basis"),
    list(path = knot_path, role = "knots"),
    list(path = coef_path, role = "coefficients")
  )

  plan$add_descriptor(fname, desc)
  plan$add_payload(basis_path, basis)
  plan$add_dataset_def(basis_path, "temporal_basis", type,
                       plan$origin_label, as.integer(step_index),
                       params_json, basis_path, "eager")
  plan$add_payload(knot_path, as.numeric(knots))
  plan$add_dataset_def(knot_path, "knots", type,
                       plan$origin_label, as.integer(step_index),
                       params_json, knot_path, "eager")
  plan$add_payload(coef_path, coeff)
  plan$add_dataset_def(coef_path, "coefficients", type, run_id,
                       as.integer(step_index), params_json,
                       coef_path, "eager")

  handle$plan <- plan
  handle$update_stash(keys = input_key,
                      new_values = list(input = coeff,
                                         temporal_coefficients = coeff))
}

# Helper to generate a simple DCT-II basis matrix
#' @keywords internal
dct_basis <- function(n_time, k) {
  idx <- seq_len(n_time) - 1/2
  mat <- outer(idx, seq_len(k), function(t, j) cos(pi * (t) * (j - 1) / n_time))
  scale <- sqrt(2 / n_time)
  mat * scale
}

#' Temporal Transform - Inverse Step
#'
#' Reconstructs data from temporal coefficients and stored basis matrix.
#' @keywords internal
invert_step.temporal <- function(type, desc, handle) {
  basis_path <- NULL
  coef_role <- NULL
  knot_path <- NULL
  if (!is.null(desc$datasets) && length(desc$datasets) > 0) {
    for (d in desc$datasets) {
      if (identical(d$role, "temporal_basis")) basis_path <- d$path
      if (identical(d$role, "coefficients")) coef_role <- d$path
      if (identical(d$role, "knots")) knot_path <- d$path
    }
  }
  if (is.null(basis_path) || is.null(coef_role)) {
    abort_lna("descriptor missing paths for temporal transform",
              .subclass = "lna_error_descriptor",
              location = "invert_step.temporal")
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  coeff <- h5_read(root, coef_role)

  subset <- handle$subset
  if (!is.null(subset$time_idx)) {
    basis <- basis[subset$time_idx, , drop = FALSE]
  }
  if (!is.null(subset$roi_mask)) {
    vox_idx <- which(as.logical(subset$roi_mask))
    coeff <- coeff[, vox_idx, drop = FALSE]
  }

  recon <- basis %*% coeff

  input_key <- desc$inputs[[1]] %||% "input"
  coeff_key <- desc$outputs[[1]] %||% "temporal_coefficients"
  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(recon), input_key))
}

