#' Analytic HRBF Transform - Forward Step (centre generation)
#'
#' Implements the centre generation and parameter bookkeeping portion of
#' the `spat.hrbf` transform. This step either analytically generates
#' RBF centres via Poisson-disk sampling or loads them from an HDF5
#' location. The generated centres and their corresponding sigma values
#' are stored in the handle's stash for use by later steps. Set
#' `num_extra_fine_levels` to generate additional dyadic scales beyond
#' `levels`.
#' @keywords internal
forward_step.spat.hrbf <- function(type, desc, handle) {
  p <- desc$params %||% list()
  sigma0 <- p$sigma0 %||% 6
  levels <- p$levels %||% 3L
  radius_factor <- p$radius_factor %||% 2.5
  num_extra_fine_levels <- p$num_extra_fine_levels %||% 0L
  kernel_type <- p$kernel_type %||% "gaussian"
  seed <- p$seed
  centres_path <- p$centres_path
  sigma_vec_path <- p$sigma_vec_path

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "forward_step.spat.hrbf:mask")
  }

  mask_arr <- as.array(mask_neurovol)
  cached_mask_world_coords <- neuroim2::grid_to_coord(
    mask_neurovol, which(mask_arr, arr.ind = TRUE)
  )
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

  # helper to convert voxel coordinates to world (mm)
  voxel_to_world <- function(vox_mat) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1,1,1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0,0,0))
    sweep(vox_mat - 1, 2, spacing_vec, `*`) + matrix(origin_vec, nrow(vox_mat), 3, byrow = TRUE)
  }

  C_total <- NULL
  sigma_vec <- NULL

  if (!is.null(seed)) {
    centres_list <- list()
    sigs <- numeric()
    for (j in seq_len(levels + 1L) - 1L) {
      sigma_j <- sigma0 / (2^j)
      r_j <- radius_factor * sigma_j
      vox_centres <- poisson_disk_sample_neuroim2(mask_neurovol, r_j, seed + j)
      if (nrow(vox_centres) > 0) {
        centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
        sigs <- c(sigs, rep(sigma_j, nrow(vox_centres)))
      }
    }
    if (num_extra_fine_levels > 0L) {
      for (j_extra in seq_len(num_extra_fine_levels)) {
        sigma_new <- sigma0 / (2^(levels + j_extra))
        r_new <- radius_factor * sigma_new
        vox_centres <- poisson_disk_sample_neuroim2(
          mask_neurovol, r_new, seed + levels + j_extra
        )
        if (nrow(vox_centres) > 0) {
          centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
          sigs <- c(sigs, rep(sigma_new, nrow(vox_centres)))
        }
      }
    }
    if (length(centres_list) > 0) {
      C_total <- do.call(rbind, centres_list)
    } else {
      C_total <- matrix(numeric(0), ncol = 3)
    }
    sigma_vec <- sigs
    p$centres_stored <- FALSE
  } else if (!is.null(centres_path) && !is.null(sigma_vec_path)) {
    root <- handle$h5[["/"]]
    C_total <- h5_read(root, centres_path)
    sigma_vec <- as.numeric(h5_read(root, sigma_vec_path))
    p$centres_stored <- TRUE
  } else {
    abort_lna("Either seed or centres_path/sigma_vec_path must be provided",
              .subclass = "lna_error_validation",
              location = "forward_step.spat.hrbf:params")
  }

  p$k_actual <- nrow(C_total)
  mask_hash_val <- digest::digest(as.array(mask_neurovol), algo = "sha256")
  p$mask_hash <- paste0("sha256:", mask_hash_val)


  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- desc$inputs %||% character()
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
  coeff <- tcrossprod(X, B_final)

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  coef_path <- file.path("/scans", run_id, "embedding", "coefficients_hrbf")
  desc$inputs <- c(desc$inputs, input_key)
  desc$outputs <- c(desc$outputs, "coefficients_hrbf")
  desc$datasets[[length(desc$datasets) + 1L]] <-
    list(path = coef_path, role = "coefficients_hrbf")

  plan$add_descriptor(fname, desc)
  plan$add_payload(coef_path, coeff)
  plan$add_dataset_def(coef_path, "coefficients_hrbf", as.character(type), run_id,
                       as.integer(step_index), params_json, coef_path,
                       "eager", dtype = NA_character_)

  handle$plan <- plan

  handle$update_stash(keys = character(),
                      new_values = list(hrbf_centres = C_total,
                                         hrbf_sigmas = sigma_vec,
                                         hrbf_basis = B_final,
                                         coefficients_hrbf = coeff))
}
#' Inverse step for the 'spat.hrbf' transform
#'
#' Reconstructs dense data from HRBF coefficients. The basis is regenerated
#' analytically unless stored in the file. A mask hash mismatch triggers
#' a warning or an error depending on
#' `lna_options("read.strict_mask_hash_validation")`.
#' @keywords internal
invert_step.spat.hrbf <- function(type, desc, handle) {
  p <- desc$params %||% list()
  sigma0 <- p$sigma0 %||% 6
  levels <- p$levels %||% 3L
  radius_factor <- p$radius_factor %||% 2.5
  kernel_type <- p$kernel_type %||% "gaussian"
  num_extra_fine_levels <- p$num_extra_fine_levels %||% 0L
  seed <- p$seed
  centres_path <- p$centres_path
  sigma_vec_path <- p$sigma_vec_path
  store_dense <- isTRUE(p$store_dense_matrix)
  centres_stored <- isTRUE(p$centres_stored)

  mask_neurovol <- handle$mask_info$mask
  if (is.null(mask_neurovol)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "invert_step.spat.hrbf:mask")
  }

  voxel_to_world <- function(vox_mat) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1,1,1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0,0,0))
    sweep(vox_mat - 1, 2, spacing_vec, `*`) +
      matrix(origin_vec, nrow(vox_mat), 3, byrow = TRUE)
  }

  if (centres_stored && !is.null(centres_path) && !is.null(sigma_vec_path)) {
    root <- handle$h5[["/"]]
    C_total <- h5_read(root, centres_path)
    sigma_vec <- as.numeric(h5_read(root, sigma_vec_path))
  } else {
    if (is.null(seed)) {
      abort_lna("seed missing for analytic regeneration",
                .subclass = "lna_error_descriptor",
                location = "invert_step.spat.hrbf:seed")
    }
    centres_list <- list()
    sigs <- numeric()
    for (j in seq_len(levels + 1L) - 1L) {
      sigma_j <- sigma0 / (2^j)
      r_j <- radius_factor * sigma_j
      vox_centres <- poisson_disk_sample_neuroim2(mask_neurovol, r_j, seed + j)
      if (nrow(vox_centres) > 0) {
        centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
        sigs <- c(sigs, rep(sigma_j, nrow(vox_centres)))
      }
    }
    if (num_extra_fine_levels > 0L) {
      for (j_extra in seq_len(num_extra_fine_levels)) {
        sigma_new <- sigma0 / (2^(levels + j_extra))
        r_new <- radius_factor * sigma_new
        vox_centres <- poisson_disk_sample_neuroim2(
          mask_neurovol, r_new, seed + levels + j_extra
        )
        if (nrow(vox_centres) > 0) {
          centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
          sigs <- c(sigs, rep(sigma_new, nrow(vox_centres)))
        }
      }
    }
    C_total <- if (length(centres_list) > 0) do.call(rbind, centres_list)
               else matrix(numeric(0), ncol = 3)
    sigma_vec <- sigs
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
