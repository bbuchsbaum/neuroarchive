#' Embed Transform - Forward Step
#'

#' Projects data onto a pre-computed basis matrix.
#' @keywords internal
forward_step.embed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  basis_path <- p$basis_path
  if (is.null(basis_path) || !nzchar(basis_path)) {
    abort_lna(
      "'basis_path' must be provided",
      .subclass = "lna_error_validation",
      location = "forward_step.embed:basis_path"
    )
  }
  plan <- handle$plan
  basis <- plan$payloads[[basis_path]]
  if (is.null(basis)) {
    abort_lna("basis matrix not found in plan payloads",
              .subclass = "lna_error_contract",
              location = "forward_step.embed:basis")
  }
  if (!is.numeric(basis)) {
    abort_lna(
      "basis matrix must be numeric",
      .subclass = "lna_error_validation",
      location = "forward_step.embed:basis"
    )
  }
  mean_vec <- if (!is.null(p$center_data_with)) plan$payloads[[p$center_data_with]] else NULL
  scale_vec <- if (!is.null(p$scale_data_with)) plan$payloads[[p$scale_data_with]] else NULL
  if (!is.null(mean_vec) && !is.numeric(mean_vec)) {
    abort_lna(
      "centering vector must be numeric",
      .subclass = "lna_error_validation",
      location = "forward_step.embed:center"
    )
  }
  if (!is.null(scale_vec) && !is.numeric(scale_vec)) {
    abort_lna(
      "scaling vector must be numeric",
      .subclass = "lna_error_validation",
      location = "forward_step.embed:scale"
    )
  }

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  X <- handle$get_inputs(input_key)[[1]]
  X <- as_dense_mat(X)


  if (!is.numeric(X)) {
    abort_lna("embed transform requires numeric input matrix",
              .subclass = "lna_error_validation",
              location = "forward_step.embed:input")
  }
  if (!is.null(mean_vec)) X <- sweep(X, 2, mean_vec, "-")
  if (!is.null(scale_vec)) X <- sweep(X, 2, scale_vec, "/")

  if (nrow(basis) == ncol(X)) {
    coeff <- X %*% basis
  } else if (ncol(basis) == ncol(X)) {
    coeff <- X %*% t(basis)
  } else {
    abort_lna(
      "basis matrix dimensions incompatible with input",
      .subclass = "lna_error_validation",
      location = "forward_step.embed"
    )
  }

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  coef_path <- paste0("/scans/", run_id, "/", base_name, "/coefficients")
  step_index <- plan$next_index
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))
  desc$params <- p

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("coefficients")

  desc$datasets <- list(list(path = coef_path, role = "coefficients"))

  plan$add_descriptor(fname, desc)
  plan$add_payload(coef_path, coeff)

  plan$add_dataset_def(coef_path, "coefficients", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       coef_path, "eager", dtype = NA_character_)

  handle$plan <- plan
  handle$update_stash(keys = input_key, new_values = list(input = coeff,
                                                         coefficients = coeff))
}



#' Embed Transform - Inverse Step
#'
#' Reconstructs data from embedding coefficients using a stored basis matrix.
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
  if (!handle$has_key(coeff_key)) {
    return(handle)
  }

  root <- handle$h5[["/"]]

  if (!path_exists_safely(root, basis_path)) {
    abort_lna(sprintf("dataset '%s' not found", basis_path),
              .subclass = "lna_error_contract",
              location = "invert_step.embed:basis")
  }
  basis <- h5_read(root, basis_path)
  if (!is.numeric(basis)) {
    abort_lna(
      "basis matrix must be numeric",
      .subclass = "lna_error_validation",
      location = "invert_step.embed:basis"
    )
  }
  mean_vec <- if (!is.null(p$center_data_with)) {
    if (!path_exists_safely(root, p$center_data_with)) {
      abort_lna(sprintf("dataset '%s' not found", p$center_data_with),
                .subclass = "lna_error_contract",
                location = "invert_step.embed:center")
    }
    h5_read(root, p$center_data_with)
  } else NULL
  scale_vec <- if (!is.null(p$scale_data_with)) {
    if (!path_exists_safely(root, p$scale_data_with)) {
      abort_lna(sprintf("dataset '%s' not found", p$scale_data_with),
                .subclass = "lna_error_contract",
                location = "invert_step.embed:scale")
    }
    h5_read(root, p$scale_data_with)
  } else NULL
  if (!is.null(mean_vec) && !is.numeric(mean_vec)) {
    abort_lna(
      "centering vector must be numeric",
      .subclass = "lna_error_validation",
      location = "invert_step.embed:center"
    )
  }
  if (!is.null(scale_vec) && !is.numeric(scale_vec)) {
    abort_lna(
      "scaling vector must be numeric",
      .subclass = "lna_error_validation",
      location = "invert_step.embed:scale"
    )
  }


  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  subset <- handle$subset
  roi_mask <- subset$roi_mask %||% subset$roi
  if (!is.null(roi_mask)) {
    vox_idx <- which(as.logical(roi_mask))
    if (nrow(basis) == ncol(coeff)) {
      basis <- basis[, vox_idx, drop = FALSE]
    } else if (ncol(basis) == ncol(coeff)) {
      basis <- basis[vox_idx, , drop = FALSE]
    }
    if (!is.null(mean_vec))  mean_vec <- mean_vec[vox_idx]
    if (!is.null(scale_vec)) scale_vec <- scale_vec[vox_idx]
  }
  time_idx <- subset$time_idx %||% subset$time
  if (!is.null(time_idx)) {
    coeff <- coeff[time_idx, , drop = FALSE]
  }

  if (nrow(basis) == ncol(coeff)) {
    dense <- coeff %*% basis
  } else if (ncol(basis) == ncol(coeff)) {
    dense <- coeff %*% t(basis)
  } else {
    abort_lna(
      "basis matrix dimensions incompatible with coefficients",
      .subclass = "lna_error_validation",
      location = "invert_step.embed"
    )
  }

  if (!is.null(scale_vec)) dense <- sweep(dense, 2, scale_vec, FUN = "*")
  if (!is.null(mean_vec))  dense <- sweep(dense, 2, mean_vec, FUN = "+")

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}


