#' Inverse step for 'basis.empirical_hrbf_compressed'
#'
#' Reconstructs an empirical basis matrix from quantized HRBF
#' dictionary codes and the stored SVD \eqn{V^T} matrix. The HRBF
#' dictionary is regenerated from the descriptor referenced by
#' `hrbf_dictionary_descriptor_path`.
#' @keywords internal
invert_step.basis.empirical_hrbf_compressed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  dict_path <- p$hrbf_dictionary_descriptor_path
  if (is.null(dict_path)) {
    abort_lna("hrbf_dictionary_descriptor_path missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.basis.empirical_hrbf_compressed")
  }

  vt_path <- NULL
  if (!is.null(desc$datasets)) {
    roles <- vapply(desc$datasets, function(d) d$role, character(1))
    idx <- which(roles == "svd_vt")
    if (length(idx) > 0) vt_path <- desc$datasets[[idx[1]]]$path
  }
  if (is.null(vt_path)) {
    abort_lna("svd_vt path not found in descriptor datasets",
              .subclass = "lna_error_descriptor",
              location = "invert_step.basis.empirical_hrbf_compressed:vt_path")
  }

  codes_key <- desc$outputs[[1]] %||% "hrbf_codes"
  input_key <- desc$inputs[[1]] %||% "basis_matrix"
  if (!handle$has_key(codes_key)) {
    return(handle)
  }
  codes <- handle$get_inputs(codes_key)[[codes_key]]

  root <- handle$h5[["/"]]
  Vt <- h5_read(root, vt_path)

  # Load HRBF dictionary descriptor
  path_parts <- strsplit(dict_path, "/")[[1]]
  dname <- tail(path_parts, 1)
  gpath <- paste(head(path_parts, -1), collapse = "/")
  if (gpath == "") gpath <- "/"
  tf_group <- root[[gpath]]
  dict_desc <- read_json_descriptor(tf_group, dname)
  if (inherits(tf_group, "H5Group")) tf_group$close()

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing",
              .subclass = "lna_error_validation",
              location = "invert_step.basis.empirical_hrbf_compressed:mask")
  }
  B_dict <- hrbf_basis_from_params(dict_desc$params, mask_neurovol)

  bits <- p$omp_quant_bits %||% 5
  if (inherits(codes, "integer")) {
    codes_num <- as.numeric(codes)
  } else {
    codes_num <- codes
  }
  codes_num <- codes_num / (2^bits - 1)

  U_sigma <- codes_num %*% B_dict
  basis_reco <- t(Vt) %*% U_sigma

  handle$update_stash(keys = codes_key,
                      new_values = setNames(list(basis_reco), input_key))
}

#' Default parameters for 'basis.empirical_hrbf_compressed'
#' @export
#' @keywords internal
lna_default.basis.empirical_hrbf_compressed <- function() {
  default_params("basis.empirical_hrbf_compressed")
}

