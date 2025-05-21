#' Embed Transform - Forward Step
#'
#' Projects data onto a pre-computed basis matrix stored in HDF5.
#' @keywords internal
forward_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path %||% ""
  if (!nzchar(basis_path)) {
    abort_lna("'basis_path' is required",
              .subclass = "lna_error_validation",
              location = "forward_step.embed:basis_path")
  }

  center_path <- p$center_data_with
  scale_path  <- p$scale_data_with

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "dense_mat"
  X <- handle$get_inputs(input_key)[[input_key]]

  if (is.array(X) && length(dim(X)) > 2) {
    d <- dim(X)
    time_dim <- d[length(d)]
    vox_dim  <- prod(d[-length(d)])
    X <- matrix(as.numeric(aperm(X, c(length(d), seq_len(length(d) - 1)))),
                nrow = time_dim, ncol = vox_dim)
  } else {
    X <- as.matrix(X)
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  if (!is.null(center_path)) {
    mvec <- as.numeric(h5_read(root, center_path))
    X <- sweep(X, 2, mvec, FUN = "-")
  }
  if (!is.null(scale_path)) {
    svec <- as.numeric(h5_read(root, scale_path))
    X <- sweep(X, 2, svec, FUN = "/")
  }

  if (nrow(basis) == ncol(X)) {
    coeff <- X %*% t(basis)
  } else if (ncol(basis) == ncol(X)) {
    coeff <- X %*% basis
  } else {
    abort_lna("basis matrix dimensions incompatible with input data",
              .subclass = "lna_error_validation",
              location = "forward_step.embed")
  }

  plan <- handle$plan
  run_id <- handle$current_run_id %||% "run-01"
  step_index <- plan$next_index
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  coef_path <- paste0("/scans/", run_id, "/", base_name, "/coefficients")
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("coefficients")
  desc$datasets <- list(list(path = coef_path, role = "coefficients"))

  plan$add_descriptor(fname, desc)
  plan$add_payload(coef_path, coeff)
  plan$add_dataset_def(coef_path, "coefficients", type, run_id,
                       as.integer(step_index), params_json,
                       coef_path, "eager")
  handle$plan <- plan

  handle <- handle$update_stash(keys = input_key,
                                new_values = setNames(list(coeff), "coefficients"))
  handle
}

#' Embed Transform - Inverse Step
#'
#' Reconstructs data from embedding coefficients using a stored basis matrix.
#' Applies ROI/time subsetting from `handle$subset`.
#' @keywords internal
invert_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path %||% ""
  if (!nzchar(basis_path)) {
    abort_lna("'basis_path' missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.embed")
  }

  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  input_key  <- desc$inputs[[1]] %||% "dense_mat"
  if (!handle$exists(coeff_key)) {
    return(handle)
  }

  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]
  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)

  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    vox_idx <- which(as.logical(subset$roi_mask))
    if (nrow(basis) == length(vox_idx)) {
      basis <- basis[vox_idx, , drop = FALSE]
    } else {
      basis <- basis[, vox_idx, drop = FALSE]
    }
  }
  if (!is.null(subset$time_idx)) {
    coeff <- coeff[subset$time_idx, , drop = FALSE]
  }

  if (nrow(basis) == ncol(coeff)) {
    dense <- tcrossprod(coeff, basis)
  } else if (ncol(basis) == ncol(coeff)) {
    dense <- coeff %*% t(basis)
  } else {
    abort_lna("basis matrix dimensions incompatible with coefficients",
              .subclass = "lna_error_validation",
              location = "invert_step.embed")
  }

  if (!is.null(p$scale_data_with)) {
    svec <- as.numeric(h5_read(root, p$scale_data_with))
    dense <- sweep(dense, 2, svec, FUN = "*")
  }
  if (!is.null(p$center_data_with)) {
    mvec <- as.numeric(h5_read(root, p$center_data_with))
    dense <- sweep(dense, 2, mvec, FUN = "+")
  }

  handle <- handle$update_stash(keys = coeff_key,
                                new_values = setNames(list(dense), input_key))
  handle
}
