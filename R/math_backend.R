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
na_get_fmrilatent_fn <- function(name, location, required = TRUE) {
  if (!requireNamespace("fmrilatent", quietly = TRUE)) {
    abort_lna(
      "HRBF transforms require the fmrilatent package",
      .subclass = "lna_error_missing_package",
      location = location
    )
  }

  fmri_ns <- asNamespace("fmrilatent")
  if (!exists(name, envir = fmri_ns, mode = "function")) {
    if (!required) {
      return(NULL)
    }
    abort_lna(
      sprintf("Installed fmrilatent does not provide `%s`; update fmrilatent", name),
      .subclass = "lna_error_missing_package",
      location = location
    )
  }

  get(name, envir = fmri_ns, mode = "function")
}

#' @keywords internal
na_hrbf_compat_params <- function(params) {
  out <- params
  # neuroarchive historically exposed the C6 polynomial under the
  # `wendland_c4` label.  Archives already contain that label, so translate it
  # at the fmrilatent boundary instead of silently changing their numerical
  # meaning now that fmrilatent distinguishes the two kernels.
  for (field in c("kernel_type", "kernel_type_fine_levels")) {
    if (identical(out[[field]], "wendland_c4")) {
      out[[field]] <- "wendland_c6"
    }
  }
  out
}

#' @keywords internal
na_hrbf_centres_from_params <- function(params, mask_neurovol, h5_root = NULL) {
  centres_fun <- na_get_fmrilatent_fn(
    "lna_hrbf_centres_from_params",
    "na_hrbf_centres_from_params",
    required = FALSE
  )

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

  if (is.null(centres_fun)) {
    return(list(
      centres = centres,
      sigmas = sigmas,
      level_vec = NULL,
      centres_stored = centres_stored,
      status = if (!is.null(centres) && !is.null(sigmas)) "present" else "unavailable_from_backend",
      source = if (!is.null(centres) && !is.null(sigmas)) "stored_inputs" else "backend_api_missing",
      authoritative = FALSE
    ))
  }

  out <- na_with_fmrilatent_compat(
    centres_fun(
      params = na_hrbf_compat_params(params),
      mask = mask_neurovol,
      centres = centres,
      sigmas = sigmas,
      compat_profile = "neuroarchive_0.1.1"
    )
  )
  out$status <- "present"
  out$source <- if (isTRUE(out$centres_stored)) "stored_inputs" else "fmrilatent"
  out$authoritative <- FALSE
  out
}

#' @keywords internal
na_hrbf_basis_from_params <- function(params, mask_neurovol, h5_root = NULL,
                                      mask_world_coords = NULL, mask_arr = NULL,
                                      mask_linear_indices = NULL) {
  basis_fun <- na_get_fmrilatent_fn(
    "lna_hrbf_basis_from_params",
    "na_hrbf_basis_from_params"
  )

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
    basis_fun(
      params = na_hrbf_compat_params(params),
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
  forward_fun <- na_get_fmrilatent_fn(
    "lna_forward_lift_matrix",
    "na_haar_forward_matrix"
  )
  na_with_fmrilatent_compat(
    forward_fun(
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
  inverse_fun <- na_get_fmrilatent_fn(
    "lna_inverse_lift_matrix",
    "na_haar_inverse_matrix"
  )
  na_with_fmrilatent_compat(
    inverse_fun(
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
