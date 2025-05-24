#' HRBF Project Transform - Forward Step
#'
#' Projects input data onto an analytic HRBF basis without storing the
#' basis or centres to HDF5. Only the resulting coefficient matrix is
#' placed in the stash for downstream transforms.
#' @keywords internal
forward_step.spat.hrbf_project <- function(type, desc, handle) {
  p <- desc$params %||% list()
  sigma0 <- p$sigma0 %||% 6
  levels <- p$levels %||% 3L
  radius_factor <- p$radius_factor %||% 2.5
  kernel_type <- p$kernel_type %||% "gaussian"
  seed <- p$seed
  centres_path <- p$centres_path
  sigma_vec_path <- p$sigma_vec_path

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "forward_step.spat.hrbf_project:mask")
  }

  if (!is.null(seed)) {
    p$centres_stored <- FALSE
  } else if (!is.null(centres_path) && !is.null(sigma_vec_path)) {
    p$centres_stored <- TRUE
  } else {
    abort_lna("Either seed or centres_path/sigma_vec_path must be provided",
              .subclass = "lna_error_validation",
              location = "forward_step.spat.hrbf_project:params")
  }

  B_final <- hrbf_basis_from_params(p, mask_neurovol,
                                 if (!is.null(handle$h5)) handle$h5[["/"]] else NULL)
  p$k_actual <- nrow(B_final)
  mask_hash_val <- digest::digest(as.array(mask_neurovol), algo = "sha256")
  p$mask_hash <- paste0("sha256:", mask_hash_val)

  inp <- handle$pull_first(c("input_dense_mat", "dense_mat", "input"))
  input_key <- inp$key
  X <- as_dense_mat(inp$value)
  #B_final_dense <- as.matrix(B_final)  # Convert sparse matrix to dense
  coeff <- Matrix::tcrossprod(X, B_final)

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("hrbf_coefficients")
  desc$datasets <- list()
  plan$add_descriptor(fname, desc)
  handle$plan <- plan

handle$update_stash(keys = character(),
                      new_values = list(hrbf_coefficients = coeff))
}

#' Inverse step for the 'spat.hrbf_project' transform
#'
#' Reconstructs dense data from HRBF coefficients using an analytically
#' generated HRBF basis. Centres are regenerated from the descriptor
#' parameters or loaded from HDF5 if provided.
#' @keywords internal
invert_step.spat.hrbf_project <- function(type, desc, handle) {
  p <- desc$params %||% list()
  sigma0 <- p$sigma0 %||% 6
  levels <- p$levels %||% 3L
  radius_factor <- p$radius_factor %||% 2.5
  kernel_type <- p$kernel_type %||% "gaussian"
  seed <- p$seed
  centres_path <- p$centres_path
  sigma_vec_path <- p$sigma_vec_path
  centres_stored <- isTRUE(p$centres_stored)

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "invert_step.spat.hrbf_project:mask")
  }

  voxel_to_world <- function(vox_mat) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1,1,1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0,0,0))
    sweep(vox_mat - 1, 2, spacing_vec, `*`) +
      matrix(origin_vec, nrow(vox_mat), 3, byrow = TRUE)
  }

  if (centres_stored && !is.null(centres_path) && !is.null(sigma_vec_path)) {
    p_local <- modifyList(p, list(centres_stored = TRUE))
    root <- handle$h5[["/"]]
    B_final <- hrbf_basis_from_params(p_local, mask_neurovol, root)
  } else {
    if (is.null(seed)) {
      abort_lna("seed missing for analytic regeneration",
                .subclass = "lna_error_descriptor",
                location = "invert_step.spat.hrbf_project:seed")
    }
    p_local <- modifyList(p, list(centres_stored = FALSE))
    B_final <- hrbf_basis_from_params(p_local, mask_neurovol,
                                   if (!is.null(handle$h5)) handle$h5[["/"]] else NULL)
  }

  coeff_key <- desc$outputs[[1]] %||% "hrbf_coefficients"
  input_key <- desc$inputs[[1]] %||% "input"
  if (!handle$has_key(coeff_key)) {
    return(handle)
  }
  coeff <- handle$get_inputs(coeff_key)[[coeff_key]]

  subset <- handle$subset
  roi_mask <- subset$roi_mask %||% subset$roi
  if (!is.null(roi_mask)) {
    vox_idx <- which(as.logical(roi_mask))
    B_final <- B_final[, vox_idx, drop = FALSE]
  }
  time_idx <- subset$time_idx %||% subset$time
  if (!is.null(time_idx)) {
    coeff <- coeff[time_idx, , drop = FALSE]
  }

  dense <- coeff %*% B_final  # Let Matrix package handle sparse operations efficiently

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}

#' Default parameters for the 'spat.hrbf_project' transform
#' @export
#' @keywords internal
lna_default.spat.hrbf_project <- function() {
  default_params("spat.hrbf_project")
}
