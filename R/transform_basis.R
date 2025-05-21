#' Basis Transform - Inverse Step
#'
#' Reconstructs data from coefficients using a stored basis matrix.
#' @keywords internal
invert_step.basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  storage_order <- p$storage_order %||% "component_x_voxel"

  basis_path <- NULL
  if (!is.null(desc$datasets) && length(desc$datasets) > 0) {
    idx <- which(vapply(desc$datasets, function(d) d$role, character(1)) == "basis_matrix")
    if (length(idx) > 0) basis_path <- desc$datasets[[idx[1]]]$path
  }
  if (is.null(basis_path)) {
    abort_lna(
      "basis_matrix path not found in descriptor",
      .subclass = "lna_error_descriptor",
      location = "invert_step.basis"
    )
  }

  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  input_key  <- desc$inputs[[1]] %||% "dense_mat"

  if (!handle$exists(coeff_key)) {
    # Nothing to reconstruct; return handle unchanged
    return(handle)
  }

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)

  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  if (identical(storage_order, "voxel_x_component")) {
    basis <- t(basis)
  }

  dense <- tcrossprod(coeff, basis)

  handle$update_stash(keys = coeff_key, new_values = setNames(list(dense), input_key))
}
