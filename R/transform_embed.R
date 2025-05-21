#' Embed Transform - Forward Step
#'
#' Projects data onto a stored basis matrix.
#' @keywords internal
forward_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path %||% ""
  center_path <- p$center_data_with
  scale_path <- p$scale_data_with

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "dense_mat"
  X <- handle$get_inputs(input_key)[[1]]

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
    abort_lna(
      "basis matrix dimensions incompatible with input data",
      .subclass = "lna_error_validation",
      location = "forward_step.embed"
    )
  }

  run_id <- handle$current_run_id %||% "run-01"

  plan <- handle$plan
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
  plan$add_dataset_def(coef_path, "coefficients", type,
                       plan$origin_label, as.integer(step_index),
                       params_json, coef_path, "eager")

  desc$params <- p
  desc$datasets <- list(list(path = coeff_path, role = "coefficients"))

  plan$add_descriptor(fname, desc)
  plan$add_payload(coeff_path, coeff)
  plan$add_dataset_def(coeff_path, "coefficients", type,
                       run_id, as.integer(step_index),
                       params_json, coeff_path, "eager")


  handle$plan <- plan
  handle$update_stash(keys = input_key,
                      new_values = list(coefficients = coeff))
}


#' Embed Transform - Inverse Step
#'
#' Reconstructs data from embedding coefficients using a stored basis matrix.
#' @keywords internal
invert_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path %||% ""
  center_path <- p$center_data_with
  scale_path <- p$scale_data_with

  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  input_key  <- desc$inputs[[1]] %||% "dense_mat"

  if (!handle$exists(coeff_key)) {
    return(handle)
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  if (nrow(basis) == ncol(coeff)) {
    dense <- tcrossprod(coeff, basis)
  } else if (ncol(basis) == ncol(coeff)) {
    dense <- coeff %*% t(basis)
  } else {
    abort_lna(
      "basis matrix dimensions incompatible with coefficients",
      .subclass = "lna_error_validation",
      location = "invert_step.embed"
    )
  }

  if (!is.null(scale_path)) {
    svec <- as.numeric(h5_read(root, scale_path))
    dense <- sweep(dense, 2, svec, FUN = "*")
  }
  if (!is.null(center_path)) {
    mvec <- as.numeric(h5_read(root, center_path))
    dense <- sweep(dense, 2, mvec, FUN = "+")
  }

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}

