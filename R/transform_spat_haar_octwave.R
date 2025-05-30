#' Inverse step for the 'spat.haar_octwave' transform
#'
#' Placeholder implementation that loads coefficient datasets and performs a
#' very basic reconstruction using `perform_haar_lift_synthesis`. The real
#' inverse lifting logic is not yet implemented.
#' @keywords internal
invert_step.spat.haar_octwave <- function(type, desc, handle) {
  mask_arr <- handle$mask_info$mask
  if (is.null(mask_arr)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "invert_step.spat.haar_octwave:mask")
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

  levels <- desc$params$levels %||% length(detail_coeffs)
  reco <- perform_haar_lift_synthesis(list(root = root_coeffs,
                                          detail = detail_coeffs),
                                      mask_arr, levels)

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
