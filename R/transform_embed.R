#' Embed Transform - Forward Step
#'
#' Projects data onto a pre-computed basis matrix.
#'
#' @param type Character string identifying the transform ("embed").
#' @param desc Descriptor list for this step. `desc$params` must include
#'   `basis_path` and may optionally specify `center_data_with` and
#'   `scale_data_with` dataset paths.
#' @param handle `DataHandle` providing access to the HDF5 file, plan and
#'   runtime stash. The input matrix is retrieved via this handle.
#'
#' @return Invisibly returns the updated `DataHandle` with the computed
#'   coefficients registered in the plan and stash.
#' @keywords internal
forward_step.embed <- function(type, desc, handle) {
  # Capture the designated output key for stashing before desc is modified
  output_key_for_chain <- desc$outputs[[1]]

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
    coeff <- tcrossprod(X, basis)
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
  desc$outputs <- c(output_key_for_chain)

  desc$datasets <- list(list(path = coef_path, role = "coefficients"))

  plan$add_descriptor(fname, desc)
  plan$add_payload(coef_path, coeff)

  plan$add_dataset_def(coef_path, "coefficients", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       coef_path, "eager", dtype = NA_character_)

  handle$plan <- plan
  
  # Use the original output_key_for_chain for stashing, not the modified desc$outputs[[1]]
  handle <- handle$update_stash(keys = input_key, 
                                new_values = setNames(list(coeff), 
                                                      output_key_for_chain))
  return(handle)
}



#' Embed Transform - Inverse Step
#'
#' Reconstructs data from embedding coefficients using a stored basis matrix.
#'
#' @param type Character string identifying the transform ("embed").
#' @param desc Descriptor list describing the inverse step. `desc$params` should
#'   contain `basis_path` along with optional `center_data_with` and
#'   `scale_data_with` dataset paths used for reconstruction.
#' @param handle `DataHandle` with access to the HDF5 file and stash containing
#'   the coefficient matrix.
#'
#' @return The updated `DataHandle` with the reconstructed dense matrix placed in
#'   the stash under `desc$inputs[[1]]`.
#' @keywords internal
invert_step.embed <- function(type, desc, handle) {
  # message("[[DEBUG]] ENTERING invert_step.embed")
  p <- desc$params %||% list()
  basis_path <- p$basis_path
  center_path <- p$center_data_with
  scale_path <- p$scale_data_with
  storage_order <- p$storage_order %||% "component_x_voxel"

  # Key to GET input from (data from invert_step.quant):
  # This should be the key that invert_step.quant used for its output.
  # invert_step.quant uses desc_quant$inputs[[1]] for its output, which is "embed_s1_out".
  # The current desc is desc_embed (from JSON). desc_embed$outputs[[1]] is now "embed_s1_out" due to previous fix.
  coeff_key <- desc$outputs[[1]] %||% "coefficients" 
  # message(sprintf("[[DEBUG invert_step.embed]] Attempting to get input from coeff_key: %s", coeff_key))
  # message(sprintf("[[DEBUG invert_step.embed]] Available stash keys in handle: %s", paste(names(handle$stash), collapse=", ")))

  # Key to PUT output to (data for invert_step.basis):
  # This is desc_embed$inputs[[1]], which is "basis_s0_out".
  input_key <- desc$inputs[[1]] %||% "data" 
  # message(sprintf("[[DEBUG invert_step.embed]] Will stash output to input_key: %s", input_key))

  if (!handle$has_key(coeff_key)) {
    # message(sprintf("[[DEBUG invert_step.embed]] coeff_key '%s' NOT FOUND in handle. Returning handle.", coeff_key))
    return(handle) # Should not happen if chaining is correct
  }

  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]
  # message(sprintf("[[DEBUG invert_step.embed]] Got input. dim(coeff): %s", paste(dim(coeff), collapse="x")))
  
  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  mean_vec <- if (!is.null(center_path)) h5_read(root, center_path) else NULL
  scale_vec <- if (!is.null(scale_path)) h5_read(root, scale_path) else NULL
  
  # Subsetting logic (simplified, assumes full data for now in this context)
  subset <- handle$subset
  time_idx <- subset$time_idx %||% subset$time
  if (!is.null(time_idx) && (is.matrix(coeff) || is.array(coeff))) {
    coeff <- coeff[time_idx, , drop = FALSE]
  }
  # ROI subsetting would apply to basis or reconstructed data, assume not active or handled by basis itself

  # Actual reconstruction
  if (identical(storage_order, "component_x_voxel")) {
    # basis is Components x Voxels (e.g. 4x8), coeff is Samples x Components (e.g. 5x4)
    # dense = coeff %*% basis -> (5x4 %*% 4x8) -> 5x8
    if (ncol(coeff) == nrow(basis)) {
      dense <- coeff %*% basis
    } else {
      abort_lna("Dimension mismatch for coeff %*% basis in invert_step.embed (comp_x_vox)", .subclass="lna_error_runtime")
    }
  } else { # voxel_x_component
    # basis is Voxels x Components (e.g. 8x4), coeff is Samples x Components (e.g. 5x4)
    # dense = coeff %*% t(basis) -> (5x4 %*% t(8x4)) -> (5x4 %*% 4x8) -> 5x8.
    if (ncol(coeff) == ncol(basis)) {
      dense <- coeff %*% t(basis)
    } else {
       abort_lna("Dimension mismatch for coeff %*% t(basis) in invert_step.embed (vox_x_comp)", .subclass="lna_error_runtime")
    }
  }
  # message(sprintf("[[DEBUG invert_step.embed]] dim(dense) after projection: %s", paste(dim(dense), collapse="x")))

  if (!is.null(scale_vec)) dense <- sweep(dense, 2, scale_vec, FUN = "*")
  if (!is.null(mean_vec))  dense <- sweep(dense, 2, mean_vec, FUN = "+")
  # message(sprintf("[[DEBUG invert_step.embed]] dim(dense) after unscale/uncenter: %s", paste(dim(dense), collapse="x")))

  handle <- handle$update_stash(keys = coeff_key,
                                new_values = setNames(list(dense), input_key))
  # message(sprintf("[[DEBUG invert_step.embed]] Stashed dense under key '%s'. Returning handle.", input_key))
  return(handle)
}


