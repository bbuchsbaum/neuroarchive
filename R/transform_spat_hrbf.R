#' Analytic HRBF Transform - Forward Step (centre generation)
#'
#' Implements the centre generation and parameter bookkeeping portion of
#' the `spat.hrbf` transform. This step either analytically generates
#' RBF centres via Poisson-disk sampling or loads them from an HDF5
#' location. The generated centres and their corresponding sigma values
#' are stored in the handle's stash for use by later steps. Set
#' `num_extra_fine_levels` to generate additional dyadic scales beyond
#' `levels`.
#' @method forward_step spat.hrbf
#' @export
#' @keywords internal
forward_step.spat.hrbf <- function(type, desc, handle) {
  p <- desc$params %||% list()

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "forward_step.spat.hrbf:mask")
  }

  if (is.array(mask_neurovol)) {
    mask_arr <- mask_neurovol
    storage.mode(mask_arr) <- "logical"
  } else {
    mask_arr <- as.array(mask_neurovol)
    # Handle case where as.array doesn't dispatch and returns a list
    if (is.list(mask_arr) && !is.null(mask_arr$arr)) {
      mask_arr <- mask_arr$arr
    }
    if (!is.array(mask_arr)) mask_arr <- as.logical(mask_arr)
    else storage.mode(mask_arr) <- "logical"
  }
  cached_mask_world_coords <- tryCatch({
    neuroim2::grid_to_coord(mask_neurovol, which(mask_arr, arr.ind = TRUE))
  }, error = function(e) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1, 1, 1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0, 0, 0))
    vox_coords <- which(mask_arr, arr.ind = TRUE)
    sweep(vox_coords - 1, 2, spacing_vec, `*`) +
      matrix(origin_vec, nrow(vox_coords), 3, byrow = TRUE)
  })
  mask_linear_indices <- as.integer(which(mask_arr))

  ## Advanced parameter validation / stubs ----
  if (isTRUE(p$use_anisotropic_atoms)) {
    if (is.null(p$anisotropy_source_path)) {
      abort_lna(
        "use_anisotropic_atoms=TRUE requires anisotropy_source_path",
        .subclass = "lna_error_validation",
        location = "forward_step.spat.hrbf:anisotropy_source"
      )
    } else if (!is.null(handle$h5)) {
      root <- handle$h5[["/"]]
      if (!path_exists_safely(root, p$anisotropy_source_path)) {
        abort_lna(
          sprintf("HDF5 path '%s' not found", p$anisotropy_source_path),
          .subclass = "lna_error_missing_path",
          location = "forward_step.spat.hrbf:anisotropy_source"
        )
      }
    }
    warn_lna(
      "Anisotropic atoms not fully implemented in this version; using isotropic."
    )
    p$use_anisotropic_atoms <- FALSE
  }

  if (identical(p$include_gaussian_derivatives, "first_order")) {
    warn_lna(
      "Derivative-of-Gaussian atoms not fully implemented; using standard Gaussians."
    )
    p$include_gaussian_derivatives <- "none"
  }

  cs_map <- p$centre_steering$map_path %||% NULL
  if (!is.null(cs_map)) {
    if (!is.null(handle$h5)) {
      root <- handle$h5[["/"]]
      if (!path_exists_safely(root, cs_map)) {
        abort_lna(
          sprintf("HDF5 path '%s' not found", cs_map),
          .subclass = "lna_error_missing_path",
          location = "forward_step.spat.hrbf:centre_steering_map"
        )
      }
    }
    warn_lna("Centre steering not fully implemented; ignoring map.")
  }

  if (isTRUE(p$use_differential_encoding)) {
    warn_lna(
      "Differential encoding not fully implemented; proceeding without it."
    )
    p$use_differential_encoding <- FALSE
  }

  centre_spec <- na_hrbf_centres_from_params(
    params = p,
    mask_neurovol = mask_neurovol,
    h5_root = if (!is.null(handle$h5)) handle$h5[["/"]] else NULL
  )
  C_total <- centre_spec$centres
  sigma_vec <- centre_spec$sigmas
  p$centres_stored <- isTRUE(centre_spec$centres_stored)

  p$k_actual <- nrow(C_total)
  mask_hash_val <- digest::digest(as.array(mask_neurovol), algo = "sha256")
  p$mask_hash <- paste0("sha256:", mask_hash_val)


  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- desc$inputs %||% character()
  # Preserve pipeline-assigned output key, or default to "coefficients_hrbf"
  pipeline_output_key <- if (!is.null(desc$outputs) && length(desc$outputs) > 0) desc$outputs[[1]] else "coefficients_hrbf"
  desc$outputs <- character()
  datasets <- list()
  desc$datasets <- datasets

  plan <- handle$plan
  step_index <- plan$next_index
  fname <- plan$get_next_filename(type)




  B_final <- hrbf_basis_from_params(
    p, mask_neurovol,
    if (!is.null(handle$h5)) handle$h5[["/"]] else NULL,
    mask_world_coords = cached_mask_world_coords,
    mask_arr = mask_arr,
    mask_linear_indices = mask_linear_indices
  )

  matrix_path <- "/basis/hrbf/analytic/matrix"
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))
  if (isTRUE(p$store_dense_matrix)) {
    plan$add_payload(matrix_path, B_final)
    plan$add_dataset_def(matrix_path, "basis_matrix", as.character(type),
                         plan$origin_label, as.integer(step_index),
                         params_json, matrix_path, "eager", dtype = NA_character_)
    desc$datasets[[length(desc$datasets) + 1L]] <-
      list(path = matrix_path, role = "basis_matrix")
  }

  inp <- handle$pull_first(c("input_dense_mat", "dense_mat", "input"))
  input_key <- inp$key
  X <- as_dense_mat(inp$value)
  coeff <- as.matrix(tcrossprod(X, B_final))

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  coef_path <- file.path("/scans", run_id, "embedding", "coefficients_hrbf")
  desc$inputs <- c(desc$inputs, input_key)
  desc$outputs <- c(desc$outputs, pipeline_output_key)
  desc$datasets[[length(desc$datasets) + 1L]] <-
    list(path = coef_path, role = "coefficients_hrbf")

  plan$add_descriptor(fname, desc)
  plan$add_payload(coef_path, coeff)
  plan$add_dataset_def(coef_path, "coefficients_hrbf", as.character(type), run_id,
                       as.integer(step_index), params_json, coef_path,
                       "eager", dtype = NA_character_)

  handle$plan <- plan

  stash_vals <- list(hrbf_centres = C_total,
                     hrbf_sigmas = sigma_vec,
                     hrbf_basis = B_final,
                     coefficients_hrbf = coeff)
  stash_vals[[pipeline_output_key]] <- coeff
  handle$update_stash(keys = character(),
                      new_values = stash_vals)
}
#' Inverse step for the 'spat.hrbf' transform
#'
#' Reconstructs dense data from HRBF coefficients. The basis is regenerated
#' analytically unless stored in the file. A mask hash mismatch triggers
#' a warning or an error depending on
#' `lna_options("read.strict_mask_hash_validation")`.
#' @method invert_step spat.hrbf
#' @export
#' @keywords internal
invert_step.spat.hrbf <- function(type, desc, handle) {
  p <- desc$params %||% list()
  store_dense <- isTRUE(p$store_dense_matrix)

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "invert_step.spat.hrbf:mask")
  }

  mask_hash_val <- digest::digest(as.array(mask_neurovol), algo = "sha256")
  current_hash <- paste0("sha256:", mask_hash_val)
  if (!is.null(p$mask_hash) && !identical(current_hash, p$mask_hash)) {
    strict <- lna_options("read.strict_mask_hash_validation")$read.strict_mask_hash_validation %||% FALSE
    msg <- sprintf("Mask hash mismatch (descriptor %s vs current %s)", p$mask_hash, current_hash)
    if (isTRUE(strict)) {
      abort_lna(msg, .subclass = "lna_error_validation",
                location = "invert_step.spat.hrbf:mask_hash")
    } else {
      warn_lna(msg, .subclass = "lna_warning_mask_hash",
               location = "invert_step.spat.hrbf:mask_hash")
    }
  }

  basis_path <- NULL
  if (store_dense && !is.null(desc$datasets)) {
    idx <- which(vapply(desc$datasets, function(d) d$role, character(1)) == "basis_matrix")
    if (length(idx) > 0) basis_path <- desc$datasets[[idx[1]]]$path
  }

  if (!is.null(basis_path)) {
    root <- handle$h5[["/"]]
    B_final <- h5_read(root, basis_path)
  } else {
    B_final <- hrbf_basis_from_params(p, mask_neurovol,
                                   if (!is.null(handle$h5)) handle$h5[["/"]] else NULL)
  }

  coeff_key <- desc$outputs[[1]] %||% "coefficients_hrbf"
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

  dense <- coeff %*% B_final

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}

#' Default parameters for the 'spat.hrbf' transform
#' @export
#' @keywords internal
lna_default.spat.hrbf <- function() {
  default_params("spat.hrbf")
}
