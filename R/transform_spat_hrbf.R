#' Analytic HRBF Transform - Forward Step (centre generation)
#'
#' Implements the centre generation and parameter bookkeeping portion of
#' the `spat.hrbf` transform. This step either analytically generates
#' RBF centres via Poisson-disk sampling or loads them from an HDF5
#' location. The generated centres and their corresponding sigma values
#' are stored in the handle's stash for use by later steps.
#' @keywords internal
forward_step.spat.hrbf <- function(type, desc, handle) {
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
              location = "forward_step.spat.hrbf:mask")
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
  mask_hash_val <- digest::digest(as.array(mask_neurovol), algo = "sha256", serialize = FALSE)
  p$mask_hash <- paste0("sha256:", mask_hash_val)

  # -- Basis Matrix Construction (HRBF-S1-5) ---------------------------------
  mask_arr <- as.array(mask_neurovol)
  mask_idx <- which(mask_arr)
  vox_coords <- which(mask_arr, arr.ind = TRUE)
  mask_coords_world <- voxel_to_world(vox_coords)

  k <- nrow(C_total)
  n_vox <- length(mask_idx)
  B_final <- matrix(0, nrow = k, ncol = n_vox)

  if (k > 0) {
    for (i in seq_len(k)) {
      atom <- generate_hrbf_atom(mask_coords_world, mask_idx,
                                 C_total[i, ], sigma_vec[i],
                                 kernel_type = kernel_type,
                                 normalize_over_mask = TRUE)
      B_final[i, ] <- atom$values
    }
  }

  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- desc$inputs %||% character()
  desc$outputs <- character()
  datasets <- list()

  plan <- handle$plan
  step_index <- plan$next_index
  fname <- plan$get_next_filename(type)

  base_name <- tools::file_path_sans_ext(fname)

  if (isTRUE(p$store_dense_matrix)) {
    matrix_path <- paste0("/basis/", base_name, "/matrix")
    datasets[[1]] <- list(path = matrix_path, role = "basis_matrix")
    params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))
    plan$add_payload(matrix_path, B_final)
    plan$add_dataset_def(matrix_path, "basis_matrix", type,
                         plan$origin_label, as.integer(step_index),
                         params_json, matrix_path, "eager", dtype = NA_character_)
  }

  desc$datasets <- datasets


  mask_arr <- as.array(mask_neurovol)
  mask_coords_vox <- which(mask_arr, arr.ind = TRUE)
  mask_coords_world <- voxel_to_world(mask_coords_vox)
  mask_linear_indices <- as.integer(which(mask_arr))
  n_total_vox <- length(mask_arr)
  k_actual <- nrow(C_total)

  if (k_actual > 0) {
    i_idx <- integer()
    j_idx <- integer()
    x_val <- numeric()
    for (kk in seq_len(k_actual)) {
      atom <- generate_hrbf_atom(mask_coords_world,
                                 mask_linear_indices,
                                 C_total[kk, ],
                                 sigma_vec[kk],
                                 kernel_type)
      i_idx <- c(i_idx, rep.int(kk, length(atom$indices)))
      j_idx <- c(j_idx, atom$indices)
      x_val <- c(x_val, atom$values)
    }
    B_final <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                                    dims = c(k_actual, n_total_vox))
  } else {
    B_final <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                    dims = c(0, n_total_vox))
  }

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

#' Default parameters for the 'spat.hrbf' transform
#' @export
#' @keywords internal
lna_default.spat.hrbf <- function() {
  default_params("spat.hrbf")
}
