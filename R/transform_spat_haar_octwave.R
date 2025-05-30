#' Forward step for the 'spat.haar_octwave' transform
#'
#' Implements the writer-side lifting analysis using the helpers
#' defined in `haar_octwave_helpers.R`. Coefficients for each level
#' are written to the plan and a concatenated matrix is stashed for
#' downstream transforms.
#' @keywords internal
forward_step.spat.haar_octwave <- function(type, desc, handle) {
  p_in <- desc$params %||% list()
  defaults <- lna_default.spat.haar_octwave()
  p <- modifyList(defaults, p_in)

  levels <- p$levels
  if (is.null(levels) || !is.numeric(levels) || levels < 1) {
    abort_lna("levels parameter required and must be >=1",
              .subclass = "lna_error_validation",
              location = "forward_step.spat.haar_octwave:levels")
  }
  levels <- as.integer(levels)
  z_seed <- p$z_order_seed %||% 42L

  mask_arr <- handle$mask_info$mask
  if (is.null(mask_arr)) {
    abort_lna("mask_info$mask missing",
              .subclass = "lna_error_validation",
              location = "forward_step.spat.haar_octwave:mask")
  }

  inp <- handle$pull_first(c("input_dense_mat", "dense_mat", "input"))
  input_key <- inp$key
  X_masked_vox_time <- convert_to_masked_vox_time_matrix(inp$value, mask_arr)

  coeffs <- perform_haar_lift_analysis(t(X_masked_vox_time), mask_arr,
                                       levels, z_seed)

  # Metadata -------------------------------------------------------------
  vox_coords <- which(as.logical(mask_arr), arr.ind = TRUE)
  bbox <- c(min(vox_coords[, 1]) - 1L, max(vox_coords[, 1]) - 1L,
            min(vox_coords[, 2]) - 1L, max(vox_coords[, 2]) - 1L,
            min(vox_coords[, 3]) - 1L, max(vox_coords[, 3]) - 1L)
  morton_idx <- get_morton_ordered_indices(mask_arr, z_seed)
  morton_hash <- morton_indices_to_hash(morton_idx)
  scalings <- precompute_haar_scalings(mask_arr, levels)
  lowpass_counts <- integer(levels + 1L)
  detail_counts <- integer(levels)
  for (lv in seq_len(levels)) {
    counts <- as.integer(round(scalings[[lv]]$sqrt_nvalid^2))
    lowpass_counts[lv] <- length(counts)
    detail_counts[lv] <- sum(counts)
  }
  lowpass_counts[levels + 1L] <- 1L

  p$num_voxels_in_mask <- length(morton_idx)
  p$octree_bounding_box_mask_space <- bbox
  p$morton_hash_mask_indices <- morton_hash
  p$num_coeffs_per_level <- list(lowpass = lowpass_counts,
                                 detail = detail_counts)

  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("wavelet_coefficients")

  datasets <- list(list(path = "/wavelet/level_ROOT/coefficients",
                        role = "wavelet_coefficients"))
  for (lv in seq_len(levels)) {
    datasets[[length(datasets) + 1L]] <-
      list(path = sprintf("/wavelet/level_%d/detail_coefficients", lv - 1L),
           role = "wavelet_coefficients")
  }
  desc$datasets <- datasets

  plan <- handle$plan
  step_index <- plan$next_index
  fname <- plan$get_next_filename(type)
  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  plan$add_descriptor(fname, desc)

  plan$add_payload(datasets[[1]]$path, coeffs$root)
  plan$add_dataset_def(datasets[[1]]$path, "wavelet_coefficients",
                       as.character(type), run_id, as.integer(step_index),
                       params_json, datasets[[1]]$path, "eager",
                       dtype = NA_character_)

  for (lv in seq_len(levels)) {
    detail_mat <- coeffs$detail[[lv]]
    dpath <- datasets[[lv + 1L]]$path
    plan$add_payload(dpath, detail_mat)
    plan$add_dataset_def(dpath, "wavelet_coefficients", as.character(type),
                         run_id, as.integer(step_index), params_json,
                         dpath, "eager", dtype = NA_character_)
  }

  handle$plan <- plan

  C_total <- cbind(coeffs$root, do.call(cbind, coeffs$detail))

  handle$update_stash(keys = character(),
                      new_values = list(wavelet_coefficients = C_total))
}

#' Inverse step for the 'spat.haar_octwave' transform
#'
#' Loads the coefficient datasets written by the forward step and
#' performs full Haar lifting synthesis to reconstruct the dense
#' voxel-time matrix.  The mask used for reconstruction is verified
#' against the descriptor's stored Morton hash.  Optional `time_idx`
#' and `roi_mask` subsetting is honoured before updating the stash.
#' @keywords internal
invert_step.spat.haar_octwave <- function(type, desc, handle) {
  p <- desc$params %||% list()
  z_seed <- p$z_order_seed %||% 42L

  mask_arr <- handle$mask_info$mask
  if (is.null(mask_arr)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "invert_step.spat.haar_octwave:mask")
  }

  stored_hash <- p$morton_hash_mask_indices
  if (!is.null(stored_hash)) {
    current_hash <- morton_indices_to_hash(
      get_morton_ordered_indices(mask_arr, z_seed)
    )
    if (!identical(stored_hash, current_hash)) {
      strict <- lna_options("read.strict_mask_hash_validation")$read.strict_mask_hash_validation %||% FALSE
      msg <- sprintf(
        "Mask Morton hash mismatch (descriptor %s vs current %s)",
        stored_hash, current_hash
      )
      if (isTRUE(strict)) {
        abort_lna(msg, .subclass = "lna_error_validation",
                  location = "invert_step.spat.haar_octwave:mask_hash")
      } else {
        warn_lna(msg, .subclass = "lna_warning_mask_hash",
                 location = "invert_step.spat.haar_octwave:mask_hash")
      }
    }
  }

  datasets <- desc$datasets %||% list()
  root_path <- NULL
  detail_paths <- list()
  for (d in datasets) {
    p <- d$path
    if (grepl("/wavelet/level_ROOT/coefficients$", p)) {
      root_path <- p
    } else if (grepl("/wavelet/level_[0-9]+/detail_coefficients$", p)) {
      lvl <- as.integer(sub(".*/level_([0-9]+)/detail_coefficients$", "\\1", p))
      detail_paths[[lvl + 1L]] <- p
    }
  }
  if (is.null(root_path)) {
    abort_lna("root coefficient dataset not found",
              .subclass = "lna_error_descriptor",
              location = "invert_step.spat.haar_octwave:root_path")
  }

  h5root <- handle$h5[["/"]]
  root_coeffs <- h5_read(h5root, root_path)
  detail_coeffs <- lapply(detail_paths, function(pp) h5_read(h5root, pp))

  subset <- handle$subset
  time_idx <- subset$time_idx %||% subset$time
  if (!is.null(time_idx)) {
    root_coeffs <- root_coeffs[time_idx, , drop = FALSE]
    detail_coeffs <- lapply(detail_coeffs, function(m) m[time_idx, , drop = FALSE])
  }

  levels <- p$levels %||% length(detail_coeffs)
  reco <- perform_haar_lift_synthesis(
    list(root = root_coeffs, detail = detail_coeffs),
    mask_arr, levels, z_seed
  )

  roi_mask <- subset$roi_mask %||% subset$roi
  if (!is.null(roi_mask)) {
    vox_idx <- which(as.logical(roi_mask))
    reco <- reco[, vox_idx, drop = FALSE]
  }

  coeff_key <- desc$outputs[[1]] %||% "wavelet_coefficients"
  output_key <- desc$inputs[[1]] %||% "input_dense_mat"

  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(reco), output_key))
}

#' Default parameters for 'spat.haar_octwave'
#' @export
#' @keywords internal
lna_default.spat.haar_octwave <- function() {
  default_params("spat.haar_octwave")
}
