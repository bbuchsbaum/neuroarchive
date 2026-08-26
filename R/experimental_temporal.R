#' Experimental Temporal Basis Functions
#'
#' @description
#' User-facing functions for temporal basis decomposition and denoising.
#' These functions provide simplified interfaces for temporal compression,
#' trend removal, and spectral analysis of neuroimaging time series.
#'
#' @name experimental_temporal
NULL

#' Create temporal basis functions
#'
#' @description
#' Generate temporal basis functions for time series decomposition.
#' Supports multiple basis types commonly used in neuroimaging with
#' optional coefficient thresholding for denoising.
#'
#' @param n_timepoints Number of time points in the series
#' @param basis_type Type of basis: "dct", "polynomial", "bspline", "dpss", "wavelet"
#' @param n_components Number of basis components to generate
#' @param ... Additional parameters specific to each basis type
#'
#' @details
#' **Basis-specific parameters:**
#' \itemize{
#'   \item \strong{DCT}: \code{threshold_type} ("none", "energy", "hard"),
#'     \code{threshold_value}, \code{keep_energy} (for energy thresholding)
#'   \item \strong{Polynomial}: No additional parameters
#'   \item \strong{B-spline}: \code{order} (default: 3), \code{knot_spacing_method}
#'   \item \strong{DPSS}: \code{time_bandwidth_product} (default: 3), \code{n_tapers}
#'   \item \strong{Wavelet}: \code{wavelet} ("db4", "db8", "haar", etc.),
#'     \code{threshold_type} ("none", "soft", "hard", "adaptive"), \code{threshold_value}
#' }
#'
#' **For fMRI applications, consider using \code{\link{suggest_dpss_fmri}} for
#' DPSS parameter selection based on TR and study type.**
#'
#' @return Matrix with time points in rows, basis functions in columns
#'
#' @examples
#' \donttest{
#' # DCT basis for trend removal with energy thresholding
#' dct_basis <- basis_time(100, "dct", n_components = 10)
#'
#' # DPSS basis for spectral analysis with optimal fMRI parameters
#' dpss_params <- suggest_dpss_fmri(TR = 2.0, n_time = 300, study_type = "resting")
#' dpss_basis <- basis_time(300, "dpss",
#'                         n_components = dpss_params$n_basis,
#'                         time_bandwidth_product = dpss_params$time_bandwidth_product)
#'
#' # B-spline basis for flexible modeling
#' spline_basis <- basis_time(150, "bspline", n_components = 20, order = 3)
#'
#' # Wavelet basis with soft thresholding for denoising
#' wavelet_basis <- basis_time(256, "wavelet", n_components = 50,
#'                           wavelet = "db4",
#'                           threshold_type = "soft")
#' }
#'
#' @seealso \code{\link{suggest_dpss_fmri}}, \code{\link{temporal_project_denoise}}
#' @export
basis_time <- function(n_timepoints,
                          basis_type = c("dct", "polynomial", "bspline", "dpss", "wavelet", "modwt"),
                          n_components = NULL,
                          ...) {

  basis_type <- match.arg(basis_type)

  # Set sensible defaults for n_components
  if (is.null(n_components)) {
    n_components <- switch(basis_type,
      "dct" = min(20, n_timepoints %/% 5),
      "polynomial" = min(5, n_timepoints %/% 10),
      "bspline" = min(n_timepoints %/% 5, 30),
      "dpss" = 7,
      "wavelet" = min(n_timepoints %/% 4, 50),
      "modwt" = NULL  # Auto-computed based on levels
    )
  }

  # Call internal temporal basis generation
  temporal_basis(basis_type, n_timepoints, n_components, ...)
}

#' Project data onto temporal basis with optional denoising
#'
#' @description
#' High-level function that combines temporal basis generation and projection
#' with optional coefficient thresholding for denoising. This provides a
#' complete workflow for temporal dimensionality reduction.
#'
#' @param data Matrix with time in rows, variables (voxels) in columns
#' @param basis_type Type of temporal basis: "dct", "polynomial", "bspline", "dpss", "wavelet"
#' @param n_components Number of basis components (NULL for automatic selection)
#' @param threshold_type Type of coefficient thresholding: "none", "energy", "hard", "soft", "adaptive"
#' @param threshold_value Threshold value (automatic if NULL for adaptive methods)
#' @param keep_energy For energy thresholding: fraction of energy to preserve (0-1)
#' @param return_basis Logical; return the basis matrix for inspection
#' @param ... Additional parameters passed to basis generation
#'
#' @return List containing:
#' \itemize{
#'   \item \code{coefficients}: Projected coefficients (n_components x n_variables)
#'   \item \code{basis}: Temporal basis matrix (if \code{return_basis = TRUE})
#'   \item \code{compression_ratio}: Actual compression achieved
#'   \item \code{denoising_stats}: Statistics on thresholding effects (if applied)
#' }
#'
#' @examples
#' \donttest{
#' # Create example fMRI time series (200 TRs, 1000 voxels)
#' set.seed(123)
#' fmri_data <- matrix(rnorm(200 * 1000), nrow = 200, ncol = 1000)
#' # Add some temporal structure
#' trend <- outer(1:200, rep(1, 1000)) * 0.01
#' fmri_data <- fmri_data + trend
#'
#' # Basic DCT compression without thresholding
#' result_basic <- temporal_project_denoise(fmri_data, "dct", n_components = 50)
#'
#' # DCT with energy-based denoising (preserve 95% of energy)
#' result_denoised <- temporal_project_denoise(
#'   fmri_data, "dct",
#'   n_components = 50,
#'   threshold_type = "energy",
#'   keep_energy = 0.95
#' )
#'
#' # DPSS for resting-state fMRI with optimal parameters
#' dpss_params <- suggest_dpss_fmri(TR = 2.0, n_time = 200, study_type = "resting")
#' result_dpss <- temporal_project_denoise(
#'   fmri_data, "dpss",
#'   n_components = dpss_params$n_basis,
#'   time_bandwidth_product = dpss_params$time_bandwidth_product,
#'   return_basis = TRUE
#' )
#'
#' # Wavelet denoising with adaptive thresholding
#' # (requires power-of-2 length, so subsample)
#' fmri_subset <- fmri_data[1:128, ]
#' result_wavelet <- temporal_project_denoise(
#'   fmri_subset, "wavelet",
#'   wavelet = "db4",
#'   threshold_type = "adaptive"
#' )
#'
#' # Reconstruct denoised data
#' if (result_denoised$compression_ratio > 1) {
#'   reconstructed <- result_dpss$basis %*% result_dpss$coefficients
#'   cat("Reconstruction error:", mean((fmri_data - reconstructed)^2), "\n")
#' }
#' }
#'
#' @export
temporal_project_denoise <- function(data,
                                    basis_type = c("dct", "polynomial", "bspline", "dpss", "wavelet", "modwt"),
                                    n_components = NULL,
                                    threshold_type = c("none", "energy", "hard", "soft", "adaptive", "bayes_shrink", "sure", "fdr"),
                                    threshold_value = NULL,
                                    keep_energy = 0.99,
                                    return_basis = FALSE,
                                    ...) {

  basis_type <- match.arg(basis_type)
  threshold_type <- match.arg(threshold_type)

  stopifnot(is.matrix(data))
  n_time <- nrow(data)
  n_vars <- ncol(data)

  # Auto-select n_components if not specified
  if (is.null(n_components)) {
    n_components <- switch(basis_type,
      "dct" = min(n_time %/% 4, 50),
      "polynomial" = min(5, n_time %/% 10),
      "bspline" = min(n_time %/% 5, 30),
      "dpss" = min(7, n_time %/% 5),
      "wavelet" = min(n_time %/% 4, 50),
      "modwt" = NULL  # Auto-computed based on decomposition levels
    )
  }

  # Generate temporal basis
  basis_args <- list(
    n_timepoints = n_time,
    basis_type = basis_type,
    n_components = n_components,
    ...
  )

  # Add thresholding parameters for relevant basis types
  if (basis_type %in% c("dct", "wavelet", "modwt") && threshold_type != "none") {
    basis_args$threshold_type <- threshold_type
    if (!is.null(threshold_value)) {
      basis_args$threshold_value <- threshold_value
    }
    if (basis_type == "dct" && threshold_type == "energy") {
      basis_args$keep_energy <- keep_energy
    }
  }

  # Generate basis (except for MODWT which projects directly)
  if (basis_type != "modwt") {
    basis <- do.call(basis_time, basis_args)
  } else {
    basis <- NULL  # MODWT doesn't use pre-computed basis
  }

  # Project data onto basis with thresholding
  if (basis_type == "dct" && threshold_type != "none") {
    coeffs <- temporal_project(
      basis_type, basis, data,
      threshold_type = threshold_type,
      threshold_value = threshold_value,
      keep_energy = keep_energy
    )
  } else if (basis_type == "wavelet" && threshold_type != "none") {
    coeffs <- temporal_project(
      basis_type, basis, data,
      threshold_type = threshold_type,
      threshold_value = threshold_value
    )
  } else if (basis_type == "modwt") {
    # MODWT projection with all parameters from basis_args
    coeffs <- temporal_project(basis_type, NULL, data, ...)
  } else {
    coeffs <- temporal_project(basis_type, basis, data)
  }

  # Calculate compression ratio
  compression_ratio <- (n_time * n_vars) / (nrow(coeffs) * n_vars + nrow(basis) * ncol(basis))

  # Calculate denoising statistics if thresholding was applied
  denoising_stats <- NULL
  if (threshold_type != "none") {
    # Count zero coefficients
    n_total_coeffs <- length(coeffs)
    n_zero_coeffs <- sum(coeffs == 0)
    sparsity <- n_zero_coeffs / n_total_coeffs

    denoising_stats <- list(
      threshold_type = threshold_type,
      sparsity = sparsity,
      n_zero_coeffs = n_zero_coeffs,
      n_total_coeffs = n_total_coeffs
    )

    if (!is.null(threshold_value)) {
      denoising_stats$threshold_value <- threshold_value
    }
    if (basis_type == "dct" && threshold_type == "energy") {
      denoising_stats$keep_energy <- keep_energy
    }
  }

  # Prepare result
  result <- list(
    coefficients = coeffs,
    compression_ratio = compression_ratio,
    basis_type = basis_type,
    n_components = ncol(basis),
    denoising_stats = denoising_stats
  )

  if (return_basis) {
    result$basis <- basis
  }

  result
}

#' Apply delta (temporal difference) transform
#'
#' @description
#' Compute temporal differences to remove slow trends and emphasize
#' rapid changes. This is commonly used for preprocessing fMRI data.
#'
#' @param data Matrix with time in rows, variables in columns
#' @param order Difference order: 1 (first differences), 2 (second differences), etc.
#' @param method Differencing method: "forward", "backward", "central"
#' @param remove_mean Remove mean before differencing (default: TRUE)
#'
#' @return Matrix of temporal differences (n_timepoints - order rows)
#'
#' @examples
#' \donttest{
#' # Create example time series data (50 timepoints, 100 variables)
#' set.seed(123)
#' my_timeseries <- matrix(rnorm(50 * 100), nrow = 50, ncol = 100)
#' # Add some trend
#' my_timeseries <- my_timeseries + outer(1:50, rep(1, 100))
#'
#' # First differences for trend removal
#' diff_data <- delta_transform(my_timeseries, order = 1)
#'
#' # Second differences for acceleration analysis
#' accel_data <- delta_transform(my_timeseries, order = 2, method = "central")
#' }
#'
#' @export
delta_transform <- function(data,
                           order = 1L,
                           method = c("forward", "backward", "central"),
                           remove_mean = TRUE) {

  method <- match.arg(method)

  if (remove_mean) {
    data <- scale(data, scale = FALSE)
  }

  # Apply differencing
  result <- data
  for (i in seq_len(order)) {
    result <- switch(method,
      "forward" = diff(result, lag = 1),
      "backward" = diff(result, lag = 1, differences = 1),
      "central" = {
        n <- nrow(result)
        if (n < 3) stop("Central differences require at least 3 time points")
        (result[3:n, , drop = FALSE] - result[1:(n-2), , drop = FALSE]) / 2
      }
    )
  }

  result
}
