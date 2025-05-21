#' Embed Transform - Forward Step
#'
#' Projects data onto a pre-computed basis matrix.
#' @keywords internal
forward_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path %||% ""
  if (!nzchar(basis_path)) {
    abort_lna("'basis_path' is required", .subclass = "lna_error_validation",
              location = "forward_step.embed:basis_path")
  }
  plan <- handle$plan
  basis <- plan$payloads[[basis_path]]
  if (is.null(basis)) {
    abort_lna("basis matrix not found in plan payloads",
              .subclass = "lna_error_contract",
              location = "forward_step.embed:basis")
  }
  mean_vec <- if (!is.null(p$center_data_with)) plan$payloads[[p$center_data_with]] else NULL
  scale_vec <- if (!is.null(p$scale_data_with)) plan$payloads[[p$scale_data_with]] else NULL

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  X <- handle$get_inputs(input_key)[[1]]
  if (is.array(X) && length(dim(X)) > 2) {
    d <- dim(X)
    tdim <- d[length(d)]
    vdim <- prod(d[-length(d)])
    X <- matrix(as.numeric(aperm(X, c(length(d), seq_len(length(d) - 1)))),
                nrow = tdim, ncol = vdim)
  } else {
    X <- as.matrix(X)
  }
  if (!is.null(mean_vec)) X <- sweep(X, 2, mean_vec, "-")
  if (!is.null(scale_vec)) X <- sweep(X, 2, scale_vec, "/")

  if (ncol(basis) == ncol(X)) {
    coeff <- X %*% basis
  } else {
    coeff <- X %*% t(basis)
  }

  run_id <- handle$current_run_id %||% "run-01"
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  coef_path <- paste0("/scans/", run_id, "/", base_name, "/coefficients")
  step_index <- plan$next_index
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
  handle$update_stash(keys = input_key, new_values = list(input = coeff,
                                                         coefficients = coeff))
}

#' Embed Transform - Inverse Step
#'
#' Reconstruct data from coefficients using a stored basis.
#' @keywords internal
invert_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path %||% ""
  if (!nzchar(basis_path)) {
    abort_lna("'basis_path' missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.embed")
  }
  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  mean_vec <- if (!is.null(p$center_data_with)) h5_read(root, p$center_data_with) else NULL
  scale_vec <- if (!is.null(p$scale_data_with)) h5_read(root, p$scale_data_with) else NULL

  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  output_key <- desc$inputs[[1]] %||% "input"

  if (!handle$exists(coeff_key)) {
    return(handle)
  }
  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  if (ncol(basis) == nrow(coeff)) {
    X <- coeff %*% basis
  } else {
    X <- coeff %*% t(basis)
  }
  if (!is.null(scale_vec)) X <- sweep(X, 2, scale_vec, "*")
  if (!is.null(mean_vec)) X <- sweep(X, 2, mean_vec, "+")

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(X), output_key))
}
