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

  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- desc$inputs %||% character()
  desc$outputs <- character()
  desc$datasets <- list()

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  plan$add_descriptor(fname, desc)
  handle$plan <- plan

  handle$update_stash(keys = character(),
                      new_values = list(hrbf_centres = C_total,
                                         hrbf_sigmas = sigma_vec))
}

#' Default parameters for the 'spat.hrbf' transform
#' @export
#' @keywords internal
lna_default.spat.hrbf <- function() {
  default_params("spat.hrbf")
}
