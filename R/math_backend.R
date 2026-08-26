#' Resolve Active Math Engine
#'
#' @return Character scalar, one of `"internal"` or `"fmrilatent"`.
#' @keywords internal
na_math_engine <- function() {
  eng <- lna_options("engine.math")[[1]] %||% "fmrilatent"
  eng <- as.character(eng[[1]])
  if (!(eng %in% c("internal", "fmrilatent"))) {
    abort_lna(
      sprintf("Unsupported engine.math value: '%s'", eng),
      .subclass = "lna_error_validation",
      location = "na_math_engine"
    )
  }
  eng
}

#' @keywords internal
na_use_fmrilatent_engine <- function() {
  identical(na_math_engine(), "fmrilatent")
}

#' @keywords internal
na_with_fmrilatent_compat <- function(expr) {
  withr::with_options(
    list(fmrilatent.compat = "neuroarchive_0.1.1"),
    expr
  )
}

#' @keywords internal
na_hrbf_centres_from_params <- function(params, mask_neurovol, h5_root = NULL) {
  if (!na_use_fmrilatent_engine()) {
    return(hrbf_centres_from_params_internal(
      params = params,
      mask_neurovol = mask_neurovol,
      h5_root = h5_root
    ))
  }

  if (!requireNamespace("fmrilatent", quietly = TRUE)) {
    abort_lna(
      "engine.math='fmrilatent' requires the fmrilatent package",
      .subclass = "lna_error_missing_package",
      location = "na_hrbf_centres_from_params"
    )
  }

  centres <- NULL
  sigmas <- NULL
  centres_stored <- isTRUE(params$centres_stored)
  centres_path <- params$centres_path
  sigma_vec_path <- params$sigma_vec_path
  if (centres_stored && !is.null(centres_path) && !is.null(sigma_vec_path)) {
    if (is.null(h5_root)) {
      abort_lna(
        "centres_stored requires h5_root to load centres/sigmas",
        .subclass = "lna_error_validation",
        location = "na_hrbf_centres_from_params"
      )
    }
    centres <- h5_read(h5_root, centres_path)
    sigmas <- as.numeric(h5_read(h5_root, sigma_vec_path))
  }

  fmri_ns <- asNamespace("fmrilatent")
  if (exists("lna_hrbf_centres_from_params", envir = fmri_ns, mode = "function")) {
    return(
      na_with_fmrilatent_compat(
        fmrilatent:::lna_hrbf_centres_from_params(
          params = params,
          mask = mask_neurovol,
          centres = centres,
          sigmas = sigmas,
          compat_profile = "neuroarchive_0.1.1"
        )
      )
    )
  }

  hrbf_centres_from_params_internal(
    params = params,
    mask_neurovol = mask_neurovol,
    h5_root = h5_root
  )
}

#' @keywords internal
na_hrbf_basis_from_params <- function(params, mask_neurovol, h5_root = NULL,
                                      mask_world_coords = NULL, mask_arr = NULL,
                                      mask_linear_indices = NULL) {
  if (!na_use_fmrilatent_engine()) {
    return(hrbf_basis_from_params_internal(
      params = params,
      mask_neurovol = mask_neurovol,
      h5_root = h5_root,
      mask_world_coords = mask_world_coords,
      mask_arr = mask_arr,
      mask_linear_indices = mask_linear_indices
    ))
  }

  if (!requireNamespace("fmrilatent", quietly = TRUE)) {
    abort_lna(
      "engine.math='fmrilatent' requires the fmrilatent package",
      .subclass = "lna_error_missing_package",
      location = "na_hrbf_basis_from_params"
    )
  }

  centres <- NULL
  sigmas <- NULL
  centres_stored <- isTRUE(params$centres_stored)
  centres_path <- params$centres_path
  sigma_vec_path <- params$sigma_vec_path
  if (centres_stored && !is.null(centres_path) && !is.null(sigma_vec_path)) {
    if (is.null(h5_root)) {
      abort_lna(
        "centres_stored requires h5_root to load centres/sigmas",
        .subclass = "lna_error_validation",
        location = "na_hrbf_basis_from_params"
      )
    }
    centres <- h5_read(h5_root, centres_path)
    sigmas <- as.numeric(h5_read(h5_root, sigma_vec_path))
  }

  na_with_fmrilatent_compat(
    fmrilatent::lna_hrbf_basis_from_params(
      params = params,
      mask = mask_neurovol,
      centres = centres,
      sigmas = sigmas,
      full_grid = TRUE,
      compat_profile = "neuroarchive_0.1.1",
      mask_world_coords = mask_world_coords,
      mask_arr = mask_arr,
      mask_linear_indices = mask_linear_indices
    )
  )
}

#' @keywords internal
na_haar_forward_matrix <- function(data_morton, mask_flat_morton, mask_dims,
                                   levels, scalings) {
  if (!na_use_fmrilatent_engine()) {
    return(NULL)
  }
  if (!requireNamespace("fmrilatent", quietly = TRUE)) {
    abort_lna(
      "engine.math='fmrilatent' requires the fmrilatent package",
      .subclass = "lna_error_missing_package",
      location = "na_haar_forward_matrix"
    )
  }
  na_with_fmrilatent_compat(
    fmrilatent::lna_forward_lift_matrix(
      data_morton = data_morton,
      mask_flat_morton = mask_flat_morton,
      mask_dims = mask_dims,
      levels = levels,
      scalings = scalings,
      compat_profile = "neuroarchive_0.1.1"
    )
  )
}

#' @keywords internal
na_haar_inverse_matrix <- function(root_coeff, detail_coeffs_by_level,
                                   mask_flat_morton, mask_dims, levels,
                                   scalings) {
  if (!na_use_fmrilatent_engine()) {
    return(NULL)
  }
  if (!requireNamespace("fmrilatent", quietly = TRUE)) {
    abort_lna(
      "engine.math='fmrilatent' requires the fmrilatent package",
      .subclass = "lna_error_missing_package",
      location = "na_haar_inverse_matrix"
    )
  }
  na_with_fmrilatent_compat(
    fmrilatent::lna_inverse_lift_matrix(
      root_coeff = root_coeff,
      detail_coeffs_by_level = detail_coeffs_by_level,
      mask_flat_morton = mask_flat_morton,
      mask_dims = mask_dims,
      levels = levels,
      scalings = scalings,
      compat_profile = "neuroarchive_0.1.1"
    )
  )
}
