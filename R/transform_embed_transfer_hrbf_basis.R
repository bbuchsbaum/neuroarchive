#' Transfer HRBF Basis - Forward Step
#'
#' Applies an empirical HRBF basis learned in a source LNA file to new
#' target data. The basis is reconstructed on-the-fly from the
#' compressed representation stored in the source file.
#' @keywords internal
forward_step.embed.transfer_hrbf_basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  src_file <- p$source_lna_file_path
  src_desc <- p$source_transform_descriptor_name
  if (is.null(src_file) || !nzchar(src_file)) {
    abort_lna("'source_lna_file_path' must be provided",
              .subclass = "lna_error_validation",
              location = "forward_step.embed.transfer_hrbf_basis:source")
  }
  if (is.null(src_desc) || !nzchar(src_desc)) {
    abort_lna("'source_transform_descriptor_name' must be provided",
              .subclass = "lna_error_validation",
              location = "forward_step.embed.transfer_hrbf_basis:descriptor")
  }
  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing",
              .subclass = "lna_error_validation",
              location = "forward_step.embed.transfer_hrbf_basis:mask")
  }

  B_emp <- .load_empirical_hrbf_basis(src_file, src_desc, mask_neurovol)

  inp <- handle$pull_first(c("input_dense_mat", "dense_mat", "input"))
  input_key <- inp$key
  X <- as_dense_mat(inp$value)

  coeff <- tcrossprod(X, B_emp)

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  scans_root <- lna_options("paths.scans_root")[[1]]
  coef_path <- paste0(scans_root, run_id, "/", base_name, "/coefficients")
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
  handle$update_stash(keys = character(),
                      new_values = list(coefficients = coeff))
}

#' Transfer HRBF Basis - Inverse Step
#'
#' Reconstructs dense data from coefficients using the HRBF basis
#' stored in a source LNA file.
#' @keywords internal
invert_step.embed.transfer_hrbf_basis <- function(type, desc, handle) {
  p <- desc$params %||% list()
  src_file <- p$source_lna_file_path
  src_desc <- p$source_transform_descriptor_name
  if (is.null(src_file) || !nzchar(src_file)) {
    abort_lna("'source_lna_file_path' missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.embed.transfer_hrbf_basis:source")
  }
  if (is.null(src_desc) || !nzchar(src_desc)) {
    abort_lna("'source_transform_descriptor_name' missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.embed.transfer_hrbf_basis:descriptor")
  }
  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing",
              .subclass = "lna_error_validation",
              location = "invert_step.embed.transfer_hrbf_basis:mask")
  }

  B_emp <- .load_empirical_hrbf_basis(src_file, src_desc, mask_neurovol)

  coeff_key <- desc$outputs[[1]] %||% "coefficients"
  input_key <- desc$inputs[[1]] %||% "dense_mat"
  if (!handle$has_key(coeff_key)) {
    return(handle)
  }
  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  subset <- handle$subset
  roi_mask <- subset$roi_mask %||% subset$roi
  if (!is.null(roi_mask)) {
    vox_idx <- which(as.logical(roi_mask))
    B_emp <- B_emp[, vox_idx, drop = FALSE]
  }
  time_idx <- subset$time_idx %||% subset$time
  if (!is.null(time_idx)) {
    coeff <- coeff[time_idx, , drop = FALSE]
  }

  dense <- coeff %*% B_emp
  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}

# ----------------------------------------------------------------------
# Helper
# ----------------------------------------------------------------------

.load_empirical_hrbf_basis <- function(file_path, desc_name, mask_neurovol) {
  h5 <- H5File$new(file_path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  tf_group <- h5[["/transforms"]]
  desc <- read_json_descriptor(tf_group, desc_name)
  safe_h5_close(tf_group)

  if (!identical(desc$type, "basis.empirical_hrbf_compressed")) {
    abort_lna("Source descriptor must be of type 'basis.empirical_hrbf_compressed'",
              .subclass = "lna_error_descriptor",
              location = ".load_empirical_hrbf_basis:type")
  }

  roles <- vapply(desc$datasets, function(d) d$role, character(1))
  vt_path <- desc$datasets[[which(roles == "svd_vt")]]$path
  codes_path <- desc$datasets[[which(roles == "hrbf_codes")]]$path
  Vt <- h5_read(h5[["/"]], vt_path)
  codes <- h5_read(h5[["/"]], codes_path)

  dict_path <- desc$params$hrbf_dictionary_descriptor_path
  dict_desc <- read_json_descriptor(h5[["/"]], dict_path)
  B_dict <- hrbf_basis_from_params(dict_desc$params, mask_neurovol)

  bits <- desc$params$omp_quant_bits %||% 5
  # codes_mat written to HDF5 already had q*scale applied during forward_step
  # So, 'codes' read here are the scaled coefficients.
  codes_num <- if (inherits(codes, "integer")) as.numeric(codes) else codes

  U_sigma <- codes_num %*% B_dict
  res <- t(Vt) %*% U_sigma
  as.matrix(res) # Ensure dense matrix output
}
