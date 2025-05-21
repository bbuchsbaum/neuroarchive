#' Embed Transform - Forward Step
#'
#' Projects data onto a pre-computed basis matrix.
#' @keywords internal
forward_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  plan <- handle$plan

  basis_path <- p$basis_path %||% ""
  if (!nzchar(basis_path)) {
    idx <- tail(which(plan$datasets$role == "basis_matrix"), 1)
    if (length(idx) == 1) {
      basis_path <- plan$datasets$path[[idx]]
      p$basis_path <- basis_path
      same_step <- plan$datasets$step_index[[idx]]
      c_idx <- which(plan$datasets$role == "center" & plan$datasets$step_index == same_step)
      if (length(c_idx) == 1 && is.null(p$center_data_with)) {
        p$center_data_with <- plan$datasets$path[[c_idx]]
      }
      s_idx <- which(plan$datasets$role == "scale" & plan$datasets$step_index == same_step)
      if (length(s_idx) == 1 && is.null(p$scale_data_with)) {
        p$scale_data_with <- plan$datasets$path[[s_idx]]
      }
    } else {
      abort_lna(
        "basis_path not specified and no prior basis_matrix found",
        .subclass = "lna_error_validation",
        location = "forward_step.embed:basis_path"
      )
    }
  }

  fetch <- function(path) {
    if (path %in% names(plan$payloads)) {
      plan$payloads[[path]]
    } else if (!is.null(handle$h5)) {
      root <- handle$h5[["/"]]
      h5_read(root, path)
    } else {
      abort_lna(
        sprintf("Dataset '%s' not available", path),
        .subclass = "lna_error_missing_path",
        location = sprintf("forward_step.embed:%s", path)
      )
    }
  }

  basis_mat <- fetch(basis_path)
  center_vec <- NULL
  scale_vec <- NULL
  if (!is.null(p$center_data_with) && nzchar(p$center_data_with)) {
    center_vec <- as.numeric(fetch(p$center_data_with))
  }
  if (!is.null(p$scale_data_with) && nzchar(p$scale_data_with)) {
    scale_vec <- as.numeric(fetch(p$scale_data_with))
  }

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "dense_mat"
  X <- handle$get_inputs(input_key)[[1]]
  if (is.array(X) && length(dim(X)) > 2) {
    d <- dim(X)
    X <- matrix(as.numeric(aperm(X, c(length(d), seq_len(length(d) - 1)))),
                nrow = d[length(d)], ncol = prod(d[-length(d)]))
  } else {
    X <- as.matrix(X)
  }

  if (!is.null(center_vec)) {
    X <- sweep(X, 2, center_vec, "-")
  }
  if (!is.null(scale_vec)) {
    X <- sweep(X, 2, scale_vec, "/")
  }

  if (nrow(basis_mat) == ncol(X)) {
    coeff <- X %*% basis_mat
  } else if (ncol(basis_mat) == ncol(X)) {
    coeff <- X %*% t(basis_mat)
  } else {
    abort_lna(
      "basis matrix dimensions incompatible with input",
      .subclass = "lna_error_validation",
      location = "forward_step.embed:dim"
    )
  }

  run_id <- handle$current_run_id %||% "run-01"
  step_index <- plan$next_index
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  coeff_path <- paste0("/scans/", run_id, "/", base_name, "/coefficients")
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("coefficients")
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
