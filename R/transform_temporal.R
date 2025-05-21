#' Temporal Transform - Forward Step
#'
#' Projects data onto a temporal basis (DCT or B-spline).
#' @keywords internal
forward_step.temporal <- function(type, desc, handle) {
  p <- desc$params %||% list()
  kind <- p$kind %||% "dct"
  n_basis <- p$n_basis
  order <- p$order %||% 3
  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"

  X <- handle$get_inputs(input_key)[[1]]
  if (is.array(X) && length(dim(X)) > 2) {
    d <- dim(X)
    tdim <- d[length(d)]
    fdim <- prod(d[-length(d)])
    X <- matrix(as.numeric(aperm(X, c(length(d), seq_len(length(d) - 1)))),
                nrow = tdim, ncol = fdim)
  } else {
    X <- as.matrix(X)
  }

  n_time <- nrow(X)
  if (is.null(n_basis)) n_basis <- n_time
  n_basis <- min(n_basis, n_time)

  if (identical(kind, "dct")) {
    basis <- .dct_basis(n_time, n_basis)
  } else if (identical(kind, "bspline")) {
    basis <- .bspline_basis(n_time, n_basis, order)
  } else {
    abort_lna(
      sprintf("Unsupported temporal kind '%s'", kind),
      .subclass = "lna_error_validation",
      location = "forward_step.temporal:kind"
    )
  }

  coeff <- crossprod(basis, X)

  run_id <- handle$current_run_id %||% "run-01"
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  basis_path <- paste0("/temporal/", base_name, "/basis")
  coef_path <- paste0("/scans/", run_id, "/", base_name, "/coefficients")
  knots_path <- paste0("/temporal/", base_name, "/knots")
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("temporal_coefficients")
  datasets <- list(
    list(path = basis_path, role = "temporal_basis"),
    list(path = coef_path, role = "temporal_coefficients")
  )
  if (identical(kind, "bspline")) {
    datasets[[length(datasets) + 1]] <- list(path = knots_path, role = "knots")
  }
  desc$datasets <- datasets

  plan$add_descriptor(fname, desc)
  plan$add_payload(basis_path, basis)
  plan$add_dataset_def(basis_path, "temporal_basis", type, run_id,
                       as.integer(plan$next_index), params_json,
                       basis_path, "eager")
  if (identical(kind, "bspline")) {
    plan$add_payload(knots_path, attr(basis, "knots"))
    plan$add_dataset_def(knots_path, "knots", type, run_id,
                         as.integer(plan$next_index), params_json,
                         knots_path, "eager")
  }
  plan$add_payload(coef_path, coeff)
  plan$add_dataset_def(coef_path, "temporal_coefficients", type, run_id,
                       as.integer(plan$next_index), params_json,
                       coef_path, "eager")
  handle$plan <- plan

  handle$update_stash(keys = input_key,
                      new_values = list(input = coeff,
                                        temporal_coefficients = coeff))
}

#' Temporal Transform - Inverse Step
#'
#' Reconstructs data from stored temporal basis coefficients.
#' @keywords internal
invert_step.temporal <- function(type, desc, handle) {
  basis_path <- NULL
  if (!is.null(desc$datasets)) {
    idx <- which(vapply(desc$datasets, function(d) d$role, character(1)) ==
                   "temporal_basis")
    if (length(idx) > 0) basis_path <- desc$datasets[[idx[1]]]$path
  }
  if (is.null(basis_path)) {
    abort_lna(
      "temporal_basis path not found in descriptor",
      .subclass = "lna_error_descriptor",
      location = "invert_step.temporal"
    )
  }

  coeff_key <- desc$outputs[[1]] %||% "temporal_coefficients"
  input_key  <- desc$inputs[[1]] %||% "input"
  if (!handle$exists(coeff_key)) {
    return(handle)
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  dense <- basis %*% coeff

  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    roi <- as.logical(subset$roi_mask)
    if (length(roi) == ncol(dense)) {
      dense <- dense[, roi, drop = FALSE]
    }
  }
  if (!is.null(subset$time_idx)) {
    idx <- as.integer(subset$time_idx)
    dense <- dense[idx, , drop = FALSE]
  }

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}

#' Generate an orthonormal DCT basis matrix
#' @keywords internal
.dct_basis <- function(n_time, n_basis) {
  t <- seq_len(n_time) - 0.5
  k <- seq_len(n_basis) - 1
  B <- sqrt(2 / n_time) * cos(outer(t, k, function(ti, ki) pi * ti * ki / n_time))
  B[,1] <- B[,1] / sqrt(2)
  B
}

#' Generate a B-spline basis matrix
#' @keywords internal
.bspline_basis <- function(n_time, n_basis, order) {
  x <- seq_len(n_time)
  splines::bs(x, df = n_basis, degree = order, intercept = TRUE)
}
