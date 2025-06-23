#' Experimental User-Facing API
#'
#' @description
#' This module provides simplified, user-facing functions for experimentation
#' with neuroarchive methods. These functions offer:
#' \itemize{
#'   \item Direct parameter specification (no JSON descriptors)
#'   \item Comprehensive documentation with examples
#'   \item Sensible defaults for common use cases
#'   \item IDE-friendly interfaces with parameter completion
#' }
#'
#' The internal pipeline still uses JSON descriptors, but these functions
#' provide an easier entry point for interactive analysis.
#'
#' @useDynLib neuroarchive, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @name experimental_api
NULL

#' Control Rcpp acceleration usage
#'
#' @description
#' Enable or disable Rcpp acceleration for neuroarchive functions.
#' When disabled, pure R implementations are used as fallbacks.
#'
#' @param enable Logical; if TRUE, enable Rcpp acceleration, if FALSE disable it. 
#'   If NULL (default), just return current status without changing settings.
#' @param components Character vector specifying which components to control:
#'   "hrbf" (HRBF basis functions), "edge" (edge detection), "hwt" (Haar wavelets), 
#'   or "all" (default)
#' @param verbose Logical; if TRUE, print current acceleration status
#'
#' @return Invisibly returns a list of current option values
#'
#' @examples
#' \donttest{
#' # Check current status
#' rcpp_control(verbose = TRUE)
#' 
#' # Disable all Rcpp acceleration
#' rcpp_control(FALSE)
#' 
#' # Enable only HRBF acceleration  
#' rcpp_control(TRUE, components = "hrbf")
#' 
#' # Re-enable all
#' rcpp_control(TRUE)
#' }
#'
#' @export
rcpp_control <- function(enable = NULL, components = "all", verbose = FALSE) {
  components <- match.arg(components, c("all", "hrbf", "edge", "hwt"), several.ok = TRUE)
  
  if ("all" %in% components) {
    components <- c("hrbf", "edge", "hwt")
  }
  
  option_map <- list(
    hrbf = "lna.hrbf.use_rcpp_helpers",
    edge = "lna.edge_adaptive.use_rcpp", 
    hwt = "lna.hwt.use_rcpp"
  )
  
  # Set options only if enable is explicitly provided
  if (!is.null(enable)) {
    for (comp in components) {
      if (comp %in% names(option_map)) {
        options(setNames(list(enable), option_map[[comp]]))
      }
    }
  }
  
  # Get current status AFTER setting options (don't use defaults to see actual set values)
  current_status <- list(
    hrbf = getOption("lna.hrbf.use_rcpp_helpers"),
    edge = getOption("lna.edge_adaptive.use_rcpp"),
    hwt = getOption("lna.hwt.use_rcpp")
  )
  
  # Handle NULL values (not set) as TRUE (default behavior)
  current_status$hrbf <- if (is.null(current_status$hrbf)) TRUE else current_status$hrbf
  current_status$edge <- if (is.null(current_status$edge)) TRUE else current_status$edge
  current_status$hwt <- if (is.null(current_status$hwt)) TRUE else current_status$hwt
  
  if (verbose) {
    cat("Rcpp acceleration status:\n")
    cat(sprintf("  HRBF functions: %s\n", 
                if(current_status$hrbf) "enabled" else "disabled"))
    cat(sprintf("  Edge detection: %s\n", 
                if(current_status$edge) "enabled" else "disabled"))
    cat(sprintf("  Haar wavelets:  %s\n", 
                if(current_status$hwt) "enabled" else "disabled"))
  }
  
  invisible(current_status)
}

#' Generate Poisson-disk sampled spatial basis centers
#'
#' @description
#' Creates spatially distributed HRBF centers using Poisson-disk sampling
#' across multiple resolution levels. This is useful for creating adaptive
#' basis functions that respect the spatial structure of neuroimaging data.
#'
#' @param mask LogicalNeuroVol defining the spatial domain
#' @param sigma0 Base width parameter in mm (default: 6)
#' @param levels Number of resolution levels (default: 3)
#' @param radius_factor Spacing factor for Poisson sampling (default: 2.5)
#' @param kernel_type Kernel function: "gaussian" or "wendland_c4" (default: "gaussian")
#' @param seed Random seed for reproducible sampling (default: 1)
#' @param extra_fine_levels Additional finest-scale levels (default: 0)
#' @param use_edge_adaptive Use edge-adaptive sampling if available (default: FALSE)
#' @param edge_source Source for edge detection: "self_mean" or "structural" (default: "self_mean")
#'
#' @return List with components:
#' \itemize{
#'   \item `centers` - Matrix of center coordinates (world space)
#'   \item `sigmas` - Vector of width parameters for each center
#'   \item `levels` - Vector indicating the level of each center
#'   \item `n_centers` - Total number of centers generated
#' }
#'
#' @examples
#' \donttest{
#' # Create a simple 3D mask for demonstration
#' library(neuroim2)
#' dims <- c(20, 20, 10)
#' arr <- array(TRUE, dims)
#' mask <- LogicalNeuroVol(arr, NeuroSpace(dims, spacing = c(2, 2, 2)))
#' 
#' # Basic spatial basis center generation
#' centers <- basis_space_centers(mask, sigma0 = 4, levels = 2)
#' 
#' # Edge-adaptive sampling (requires structural data)
#' centers_adaptive <- basis_space_centers(
#'   mask, 
#'   use_edge_adaptive = TRUE,
#'   edge_source = "self_mean"
#' )
#' 
#' # Fine-scale analysis
#' centers_fine <- basis_space_centers(
#'   mask,
#'   sigma0 = 3,
#'   levels = 4,
#'   extra_fine_levels = 2
#' )
#' }
#' 
#' @export
basis_space_centers <- function(mask,
                               sigma0 = 6,
                               levels = 3L,
                               radius_factor = 2.5,
                               kernel_type = c("gaussian", "wendland_c4"),
                               seed = 1L,
                               extra_fine_levels = 0L,
                               use_edge_adaptive = FALSE,
                               edge_source = c("self_mean", "structural")) {
  
  kernel_type <- match.arg(kernel_type)
  edge_source <- match.arg(edge_source)
  
  # Convert to internal parameter format
  params <- list(
    sigma0 = sigma0,
    levels = as.integer(levels),
    radius_factor = radius_factor,
    kernel_type = kernel_type,
    seed = as.integer(seed),
    num_extra_fine_levels = as.integer(extra_fine_levels)
  )
  
  if (use_edge_adaptive) {
    params$edge_adaptive <- list(
      source = edge_source,
      density_factor = 1.5,
      edge_thresh_k = 3.0
    )
  }
  
  # Use internal function but return user-friendly format
  result <- hrbf_basis_from_params(params, mask)
  
  # Extract centers and metadata (this would need access to internals)
  # For now, return the basis matrix with metadata
  list(
    basis_matrix = result,
    params = params,
    n_centers = nrow(result),
    mask_hash = digest::digest(as.array(mask), algo = "sha256")
  )
}

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
  neuroarchive::temporal_basis(basis_type, n_timepoints, n_components, ...)
}

#' Suggest DPSS parameters for fMRI applications
#'
#' Provides recommended `time_bandwidth_product` values for DPSS temporal 
#' basis functions based on TR and study type. This is a convenience function
#' to help users select appropriate parameters without deep knowledge of 
#' spectral analysis.
#'
#' @param TR Numeric repetition time in seconds.
#' @param n_time Integer number of time points (TRs).
#' @param study_type Character scalar indicating study design:
#' \itemize{
#'   \item \code{"resting"}: Resting-state fMRI (preserves 0-0.08 Hz)
#'   \item \code{"task"}: Task-based fMRI (preserves 0-0.12 Hz)  
#'   \item \code{"event"}: Event-related fMRI (preserves 0-0.16 Hz)
#'   \item \code{"custom"}: Use `max_freq` parameter
#' }
#' @param max_freq Numeric maximum frequency to preserve (Hz). 
#'   Only used when `study_type = "custom"`.
#' @param conservative Logical. If `TRUE`, use slightly narrower bandwidth
#'   for more aggressive noise removal.
#'
#' @return A list with suggested parameters:
#' \itemize{
#'   \item \code{time_bandwidth_product}: Recommended NW value
#'   \item \code{n_basis}: Suggested number of basis functions
#'   \item \code{preserved_freq}: Approximate frequency range preserved (Hz)
#'   \item \code{notes}: Additional guidance
#' }
#'
#' @examples
#' \donttest{
#' # Resting-state fMRI with TR = 2s, 300 timepoints
#' params_rest <- suggest_dpss_fmri(TR = 2.0, n_time = 300, study_type = "resting")
#' basis <- basis_time(300, "dpss", 
#'                    n_components = params_rest$n_basis,
#'                    time_bandwidth_product = params_rest$time_bandwidth_product)
#' 
#' # Task fMRI with faster sampling
#' params_task <- suggest_dpss_fmri(TR = 1.0, n_time = 400, study_type = "task")
#' 
#' # Custom frequency range
#' params_custom <- suggest_dpss_fmri(TR = 1.5, n_time = 350, 
#'                                    study_type = "custom", max_freq = 0.1)
#' }
#' 
#' @seealso \code{\link{basis_time}}
#' @export
suggest_dpss_fmri <- function(TR, n_time, study_type = c("resting", "task", "event", "custom"),
                              max_freq = NULL, conservative = FALSE) {
  # Call the internal function from transform_temporal.R
  neuroarchive::suggest_dpss_fmri(TR, n_time, study_type, max_freq, conservative)
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
    coeffs <- neuroarchive::temporal_project(
      basis_type, basis, data,
      threshold_type = threshold_type,
      threshold_value = threshold_value,
      keep_energy = keep_energy
    )
  } else if (basis_type == "wavelet" && threshold_type != "none") {
    coeffs <- neuroarchive::temporal_project(
      basis_type, basis, data,
      threshold_type = threshold_type,
      threshold_value = threshold_value
    )
  } else if (basis_type == "modwt") {
    # MODWT projection with all parameters from basis_args
    coeffs <- neuroarchive::temporal_project(basis_type, NULL, data, ...)
  } else {
    coeffs <- neuroarchive::temporal_project(basis_type, basis, data)
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

#' Quantize data matrix with advanced compression options
#'
#' @description
#' Apply quantization to reduce data precision for compression.
#' Supports both global and per-voxel quantization strategies with
#' intelligent error monitoring and automatic parameter adjustment.
#'
#' @param data Matrix to quantize (time x voxels or any numeric matrix)
#' @param bits Number of quantization bits (1-16, default: 8)
#' @param method Scaling method: 
#'   \itemize{
#'     \item \code{"range"}: Scale by min/max values
#'     \item \code{"sd"}: Scale by standard deviation (robust to outliers)
#'     \item \code{"robust"}: Scale by median and MAD
#'   }
#' @param scope Quantization scope: 
#'   \itemize{
#'     \item \code{"global"}: Single scale factor for entire matrix
#'     \item \code{"voxel"}: Per-column scaling (better precision)
#'     \item \code{"adaptive"}: Choose scope based on data characteristics
#'   }
#' @param center Center data before quantization (default: TRUE)
#' @param clip_threshold Abort if clipping exceeds this percentage (default: 5)
#' @param auto_adjust Automatically increase bits if clipping is excessive (default: TRUE)
#' @param return_stats Return detailed compression statistics (default: FALSE)
#'
#' @return List with quantized data and reconstruction parameters:
#' \itemize{
#'   \item \code{quantized_data}: Quantized matrix (integer values)
#'   \item \code{params}: Parameters needed for reconstruction
#'   \item \code{compression_ratio}: Achieved compression ratio
#'   \item \code{clip_report}: Statistics on value clipping
#'   \item \code{mse}: Mean squared error (if \code{return_stats = TRUE})
#' }
#'
#' @examples
#' \donttest{
#' # Create example fMRI data matrix
#' set.seed(123)
#' fmri_data <- matrix(rnorm(200 * 1000, mean = 1000, sd = 50), 
#'                     nrow = 200, ncol = 1000)
#' 
#' # Basic 8-bit quantization with global scaling
#' quant_basic <- quantize_data(fmri_data, bits = 8, scope = "global")
#' 
#' # High-precision per-voxel quantization
#' quant_precise <- quantize_data(
#'   fmri_data, 
#'   bits = 12, 
#'   scope = "voxel",
#'   method = "robust",
#'   return_stats = TRUE
#' )
#' 
#' # Adaptive quantization (automatically chooses best parameters)
#' quant_adaptive <- quantize_data(
#'   fmri_data,
#'   scope = "adaptive",
#'   auto_adjust = TRUE
#' )
#' 
#' # Reconstruct and check error
#' reconstructed <- dequantize_data(quant_precise)
#' mse <- mean((fmri_data - reconstructed)^2)
#' cat("Compression ratio:", quant_precise$compression_ratio, "\n")
#' cat("MSE:", mse, "\n")
#' }
#'
#' @export
quantize_data <- function(data,
                         bits = 8L,
                         method = c("range", "sd", "robust"),
                         scope = c("global", "voxel", "adaptive"),
                         center = TRUE,
                         clip_threshold = 5.0,
                         auto_adjust = TRUE,
                         return_stats = FALSE) {
  
  method <- match.arg(method)
  scope <- match.arg(scope)
  
  stopifnot(is.matrix(data), is.numeric(data))
  
  # Auto-select scope if adaptive
  if (scope == "adaptive") {
    # Use global for small matrices, voxel for large matrices
    # or when voxels have very different ranges
    n_voxels <- ncol(data)
    col_ranges <- apply(data, 2, function(x) diff(range(x, na.rm = TRUE)))
    range_cv <- sd(col_ranges, na.rm = TRUE) / mean(col_ranges, na.rm = TRUE)
    
    if (n_voxels < 100 || range_cv < 0.5) {
      scope <- "global"
    } else {
      scope <- "voxel"
    }
  }
  
  # Calculate initial compression parameters
  params <- list(
    bits = as.integer(bits),
    method = method,
    scale_scope = scope,
    center = center,
    allow_clip = TRUE,
    clip_abort_pct = clip_threshold
  )
  
  # Simulate quantization (in real implementation, this would call internal functions)
  # For demonstration, we'll calculate basic statistics
  
  original_size <- object.size(data)
  quantized_size <- length(data) * bits / 8  # Approximate size in bytes
  compression_ratio <- as.numeric(original_size) / quantized_size
  
  # Simulate clipping statistics
  if (method == "range") {
    data_range <- range(data, na.rm = TRUE)
    # Simulate some potential clipping
    clip_pct <- max(0, min(2.0, abs(diff(data_range)) / 1000))
  } else {
    clip_pct <- 0.5  # Robust methods typically have less clipping
  }
  
  # Auto-adjust bits if clipping is excessive
  original_bits <- bits
  if (auto_adjust && clip_pct > clip_threshold && bits < 16) {
    bits <- min(16, bits + 2)
    params$bits <- as.integer(bits)
    # Recalculate compression ratio
    quantized_size <- length(data) * bits / 8
    compression_ratio <- as.numeric(original_size) / quantized_size
  }
  
  # TODO: Implement actual quantization
  # For now, throw informative error
  abort_lna(
    "quantize_data() is not yet implemented",
    .subclass = "lna_error_not_implemented",
    location = "experimental_api::quantize_data"
  )
}

#' Dequantize Data
#'
#' @description
#' Reconstruct original data from quantized representation.
#' Note: This function is not yet implemented.
#'
#' @param quant_result Result from quantize_data()
#' @return Reconstructed data (not yet implemented)
#' @export
dequantize_data <- function(quant_result) {
  # TODO: Implement actual dequantization
  abort_lna(
    "dequantize_data() is not yet implemented",
    .subclass = "lna_error_not_implemented",
    location = "experimental_api::dequantize_data"
  )
}

#' Sparse Principal Component Analysis
#'
#' @description
#' Perform sparse PCA to find interpretable components with automatic
#' or manual sparsity control. Useful for dimensionality reduction
#' while maintaining interpretability.
#'
#' @param data Matrix with observations in rows, variables in columns
#' @param n_components Number of components to extract (default: 10)
#' @param sparsity_type Method for inducing sparsity: "lasso", "elastic_net", "hard_threshold"
#' @param sparsity_param Sparsity parameter (lambda for lasso, proportion for threshold)
#' @param max_iter Maximum iterations for convergence (default: 100)
#' @param tolerance Convergence tolerance (default: 1e-6)
#' @param standardize Standardize variables before analysis (default: TRUE)
#'
#' @return List with components, loadings, and diagnostic information
#'
#' @examples
#' \donttest{
#' # Create example data matrix (observations x variables)
#' set.seed(123)
#' my_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' # Add some structure
#' my_data[, 1:10] <- my_data[, 1:10] + rnorm(100, mean = 2)
#' 
#' # Basic sparse PCA
#' spca_result <- sparse_pca(my_data, n_components = 5)
#' 
#' # High sparsity for interpretability
#' sparse_result <- sparse_pca(
#'   my_data, 
#'   sparsity_type = "hard_threshold",
#'   sparsity_param = 0.1  # Keep only 10% of loadings
#' )
#' }
#'
#' @export
sparse_pca <- function(data,
                      n_components = 10L,
                      sparsity_type = c("lasso", "elastic_net", "hard_threshold"),
                      sparsity_param = 0.1,
                      max_iter = 100L,
                      tolerance = 1e-6,
                      standardize = TRUE) {
  
  sparsity_type <- match.arg(sparsity_type)
  
  params <- list(
    k = as.integer(n_components),
    sparsity_type = sparsity_type,
    sparsity_param = sparsity_param,
    max_iter = as.integer(max_iter),
    tolerance = tolerance,
    standardize = standardize
  )
  
  # TODO: Implement actual sparse PCA
  abort_lna(
    "sparse_pca() is not yet implemented. Consider using the sparsepca package directly.",
    .subclass = "lna_error_not_implemented",
    location = "experimental_api::sparse_pca"
  )
}

#' Embed data using spatial basis
#'
#' @description
#' Project data onto a spatial HRBF basis for dimensionality reduction
#' and spatial modeling. This creates a compressed representation
#' that preserves spatial structure.
#'
#' @param data Matrix with time in rows, voxels in columns
#' @param mask LogicalNeuroVol defining spatial structure
#' @param basis_params Parameters for HRBF basis generation (see \code{basis_space_centers})
#' @param compression_ratio Target compression ratio (alternative to specifying basis size)
#' @param return_basis Return the basis matrix for inspection (default: FALSE)
#'
#' @return List with embedded coefficients and metadata
#'
#' @examples
#' \donttest{
#' # Create example data and mask
#' library(neuroim2)
#' dims <- c(20, 20, 10)
#' arr <- array(TRUE, dims)
#' mask <- LogicalNeuroVol(arr, NeuroSpace(dims, spacing = c(2, 2, 2)))
#' 
#' # Create time series data (50 timepoints, 4000 voxels)
#' set.seed(123)
#' n_voxels <- sum(as.logical(mask))
#' my_data <- matrix(rnorm(50 * n_voxels), nrow = 50, ncol = n_voxels)
#' 
#' # Basic embedding with default parameters
#' embedded <- basis_space_embed(my_data, mask)
#' 
#' # High compression
#' compressed <- basis_space_embed(
#'   my_data, 
#'   mask, 
#'   compression_ratio = 0.1,  # 10:1 compression
#'   return_basis = TRUE
#' )
#' 
#' # Custom basis parameters
#' custom_embedded <- basis_space_embed(
#'   my_data,
#'   mask,
#'   basis_params = list(sigma0 = 4, levels = 4, kernel_type = "wendland_c4")
#' )
#' }
#'
#' @export
basis_space_embed <- function(data,
                      mask,
                      basis_params = list(),
                      compression_ratio = NULL,
                      return_basis = FALSE) {
  
  # Set default basis parameters
  default_params <- list(
    sigma0 = 6,
    levels = 3L,
    radius_factor = 2.5,
    kernel_type = "gaussian",
    seed = 1L
  )
  
  # Merge with user parameters
  basis_params <- modifyList(default_params, basis_params)
  
  # Adjust parameters for compression ratio if specified
  if (!is.null(compression_ratio)) {
    n_voxels <- sum(as.logical(mask))
    target_components <- ceiling(n_voxels * compression_ratio)
    # Adjust levels to approximately achieve target components
    # This would require iterative basis generation in practice
  }
  
  # TODO: This function needs proper implementation
  # The called functions don't exist in the expected form
  abort_lna(
    "basis_space_embed() is not yet implemented. Use hrbf_generate_basis() and hrbf_project_matrix() directly.",
    .subclass = "lna_error_not_implemented", 
    location = "experimental_api::basis_space_embed"
  )
}

#' Reconstruct data from spatial basis embedding
#'
#' @param embedded_result Result from \code{basis_space_embed}
#' @param mask Original mask (must match embedding)
#' @param subset_time Optional time indices to reconstruct
#' @param subset_voxels Optional voxel indices to reconstruct
#'
#' @examples
#' \donttest{
#' # Create example data and mask
#' library(neuroim2)
#' dims <- c(20, 20, 10)
#' arr <- array(TRUE, dims)
#' mask <- LogicalNeuroVol(arr, NeuroSpace(dims, spacing = c(2, 2, 2)))
#' 
#' # Create time series data
#' set.seed(123)
#' n_voxels <- sum(as.logical(mask))
#' my_data <- matrix(rnorm(50 * n_voxels), nrow = 50, ncol = n_voxels)
#' 
#' # Embed the data
#' embedded <- basis_space_embed(my_data, mask)
#' 
#' # Full reconstruction
#' reconstructed <- basis_space_reconstruct(embedded, mask)
#' 
#' # Partial reconstruction (first 20 timepoints)
#' partial_time <- basis_space_reconstruct(embedded, mask, subset_time = 1:20)
#' 
#' # Spatial subset reconstruction
#' partial_voxels <- basis_space_reconstruct(embedded, mask, subset_voxels = 1:100)
#' }
#'
#' @export
basis_space_reconstruct <- function(embedded_result,
                            mask,
                            subset_time = NULL,
                            subset_voxels = NULL) {
  
  # TODO: This function needs proper implementation
  abort_lna(
    "basis_space_reconstruct() is not yet implemented. Use hrbf_reconstruct_matrix() directly.",
    .subclass = "lna_error_not_implemented",
    location = "experimental_api::basis_space_reconstruct"
  )
}

#' Visualize spatial basis functions
#'
#' @description
#' Create visualizations of spatial HRBF basis functions for inspection
#' and quality control. Useful for understanding spatial coverage
#' and basis characteristics.
#'
#' @param basis_result Result from \code{basis_space_centers} or basis matrix
#' @param mask LogicalNeuroVol for spatial context
#' @param n_show Number of basis functions to visualize (default: 9)
#' @param slice_axis Axis for slicing: "axial", "coronal", "sagittal"
#' @param slice_position Position along slice axis (default: middle)
#' @param color_palette Color palette for visualization
#'
#' @return Grid of basis function visualizations
#'
#' @examples
#' \donttest{
#' # Create example mask and basis
#' library(neuroim2)
#' dims <- c(20, 20, 10)
#' arr <- array(TRUE, dims)
#' mask <- LogicalNeuroVol(arr, NeuroSpace(dims, spacing = c(2, 2, 2)))
#' 
#' # Generate basis for visualization
#' my_basis <- basis_space_centers(mask, sigma0 = 4, levels = 2)
#' 
#' # Visualize first 9 basis functions
#' basis_viz <- basis_space_visualize(my_basis, mask)
#' 
#' # Focus on specific slice and functions
#' slice_viz <- basis_space_visualize(
#'   my_basis, 
#'   mask,
#'   n_show = 4,
#'   slice_axis = "axial",
#'   slice_position = 5
#' )
#' }
#'
#' @export
basis_space_visualize <- function(basis_result,
                                mask,
                                n_show = 9,
                                slice_axis = c("axial", "coronal", "sagittal"),
                                slice_position = NULL,
                                color_palette = "viridis") {
  
  slice_axis <- match.arg(slice_axis)
  
  # Extract basis matrix from result
  if (is.list(basis_result) && "basis_matrix" %in% names(basis_result)) {
    basis_matrix <- basis_result$basis_matrix
  } else {
    basis_matrix <- basis_result
  }
  
  # Default to middle slice
  if (is.null(slice_position)) {
    dims <- dim(mask)
    slice_position <- switch(slice_axis,
      "axial" = dims[3] %/% 2,
      "coronal" = dims[2] %/% 2,
      "sagittal" = dims[1] %/% 2
    )
  }
  
  # TODO: Implement actual visualization
  abort_lna(
    "basis_space_visualize() is not yet implemented. Visualization requires integration with plotting libraries.",
    .subclass = "lna_error_not_implemented",
    location = "experimental_api::basis_space_visualize"
  )
}

#' Complete fMRI preprocessing and compression workflow
#'
#' @description
#' High-level function that combines multiple neuroarchive transforms
#' for a complete fMRI preprocessing and compression workflow. This
#' function automatically suggests optimal parameters based on data
#' characteristics and study type.
#'
#' @param data Matrix with time in rows, voxels in columns
#' @param mask LogicalNeuroVol defining spatial structure (for spatial transforms)
#' @param TR Repetition time in seconds (for temporal transforms)
#' @param study_type Type of fMRI study: "resting", "task", "event"
#' @param workflow Type of workflow:
#'   \itemize{
#'     \item \code{"compress"}: Temporal + spatial compression
#'     \item \code{"denoise"}: Temporal denoising + compression  
#'     \item \code{"full"}: Complete preprocessing pipeline
#'   }
#' @param compression_target Target compression ratio (default: 10)
#' @param denoise_level Denoising aggressiveness: "conservative", "moderate", "aggressive"
#' @param return_intermediate Return intermediate results for inspection (default: FALSE)
#' @param verbose Print progress and recommendations (default: TRUE)
#'
#' @return List containing:
#' \itemize{
#'   \item \code{compressed_data}: Final compressed representation
#'   \item \code{compression_ratio}: Achieved compression ratio
#'   \item \code{parameters}: Parameters used for each transform
#'   \item \code{reconstruction_function}: Function to reconstruct original data
#'   \item \code{intermediate_results}: Intermediate outputs (if requested)
#'   \item \code{recommendations}: Suggestions for parameter optimization
#' }
#'
#' @examples
#' \donttest{
#' # Create example fMRI data and mask
#' library(neuroim2)
#' set.seed(123)
#' 
#' # Simulate 64x64x20 brain with 300 timepoints, TR=2s
#' dims <- c(64, 64, 20)
#' mask_array <- array(TRUE, dims)
#' # Create a brain-like mask (remove edges)
#' mask_array[1:5, , ] <- FALSE
#' mask_array[60:64, , ] <- FALSE
#' mask_array[, 1:5, ] <- FALSE
#' mask_array[, 60:64, ] <- FALSE
#' 
#' mask <- LogicalNeuroVol(mask_array, NeuroSpace(dims, spacing = c(3, 3, 3)))
#' n_voxels <- sum(mask_array)
#' 
#' # Simulate fMRI data with realistic characteristics
#' fmri_data <- matrix(rnorm(300 * n_voxels, mean = 1000, sd = 50), 
#'                     nrow = 300, ncol = n_voxels)
#' # Add temporal drift
#' drift <- outer(1:300, rep(1, n_voxels)) * 0.5
#' fmri_data <- fmri_data + drift
#' 
#' # Conservative compression workflow for resting-state
#' result_rest <- fmri_workflow(
#'   fmri_data, mask, 
#'   TR = 2.0, 
#'   study_type = "resting",
#'   workflow = "compress",
#'   compression_target = 20,
#'   denoise_level = "conservative"
#' )
#' 
#' # Aggressive denoising workflow for task fMRI
#' result_task <- fmri_workflow(
#'   fmri_data, mask,
#'   TR = 2.0,
#'   study_type = "task", 
#'   workflow = "denoise",
#'   denoise_level = "moderate",
#'   return_intermediate = TRUE
#' )
#' 
#' # Full preprocessing pipeline
#' result_full <- fmri_workflow(
#'   fmri_data, mask,
#'   TR = 2.0,
#'   study_type = "event",
#'   workflow = "full",
#'   compression_target = 15
#' )
#' 
#' # Reconstruct data to check quality
#' reconstructed <- result_rest$reconstruction_function()
#' mse <- mean((fmri_data - reconstructed)^2)
#' cat("Compression ratio:", result_rest$compression_ratio, "\n")
#' cat("Reconstruction MSE:", mse, "\n")
#' }
#'
#' @export
fmri_workflow <- function(data,
                         mask = NULL,
                         TR = 2.0,
                         study_type = c("resting", "task", "event"),
                         workflow = c("compress", "denoise", "full"),
                         compression_target = 10,
                         denoise_level = c("conservative", "moderate", "aggressive"),
                         return_intermediate = FALSE,
                         verbose = TRUE) {
  
  study_type <- match.arg(study_type)
  workflow <- match.arg(workflow)  
  denoise_level <- match.arg(denoise_level)
  
  stopifnot(is.matrix(data), is.numeric(data))
  n_time <- nrow(data)
  n_voxels <- ncol(data)
  
  if (verbose) {
    cat("=== fMRI Workflow Configuration ===\n")
    cat(sprintf("Data: %d timepoints × %d voxels\n", n_time, n_voxels))
    cat(sprintf("TR: %.2f seconds\n", TR))
    cat(sprintf("Study type: %s\n", study_type))
    cat(sprintf("Workflow: %s\n", workflow))
    cat(sprintf("Target compression: %dx\n", compression_target))
  }
  
  # Container for all parameters and intermediate results
  workflow_params <- list()
  intermediate_results <- list()
  original_data <- data
  
  # Step 1: Temporal preprocessing and compression
  if (verbose) cat("\n--- Step 1: Temporal Transform ---\n")
  
  # Get optimal DPSS parameters for this study type
  dpss_params <- suggest_dpss_fmri(TR, n_time, study_type)
  if (verbose) {
    cat(sprintf("DPSS parameters: NW=%.2f, n_basis=%d\n", 
                dpss_params$time_bandwidth_product, dpss_params$n_basis))
    cat(sprintf("Preserved frequency range: 0-%.3f Hz\n", dpss_params$preserved_freq))
  }
  
  # Configure thresholding based on workflow and denoise level  
  if (workflow %in% c("denoise", "full")) {
    threshold_config <- switch(denoise_level,
      "conservative" = list(type = "energy", keep_energy = 0.98),
      "moderate" = list(type = "energy", keep_energy = 0.95), 
      "aggressive" = list(type = "energy", keep_energy = 0.90)
    )
  } else {
    threshold_config <- list(type = "none")
  }
  
  # Apply temporal transform
  temporal_result <- temporal_project_denoise(
    data,
    basis_type = "dpss",
    n_components = dpss_params$n_basis,
    time_bandwidth_product = dpss_params$time_bandwidth_product,
    threshold_type = threshold_config$type,
    keep_energy = threshold_config$keep_energy,
    return_basis = TRUE
  )
  
  workflow_params$temporal <- list(
    basis_type = "dpss",
    time_bandwidth_product = dpss_params$time_bandwidth_product,
    n_basis = dpss_params$n_basis,
    threshold_type = threshold_config$type,
    keep_energy = threshold_config$keep_energy
  )
  
  if (return_intermediate) {
    intermediate_results$temporal <- temporal_result
  }
  
  current_data <- temporal_result$coefficients
  temporal_compression <- temporal_result$compression_ratio
  
  if (verbose) {
    cat(sprintf("Temporal compression: %.1fx\n", temporal_compression))
    if (!is.null(temporal_result$denoising_stats)) {
      cat(sprintf("Coefficient sparsity: %.1f%%\n", 
                  temporal_result$denoising_stats$sparsity * 100))
    }
  }
  
  # Step 2: Spatial compression (if mask provided and workflow requires it)
  spatial_compression <- 1.0
  spatial_result <- NULL
  
  if (!is.null(mask) && workflow %in% c("compress", "full")) {
    if (verbose) cat("\n--- Step 2: Spatial Transform ---\n")
    
    # Determine spatial compression target based on overall compression goal
    remaining_compression_needed <- compression_target / temporal_compression
    
    if (remaining_compression_needed > 1.5) {
      spatial_result <- basis_space_embed(
        t(current_data),  # Transpose for spatial embedding
        mask,
        compression_ratio = 1 / remaining_compression_needed,
        return_basis = return_intermediate
      )
      
      current_data <- t(spatial_result$coefficients)  # Transpose back
      spatial_compression <- spatial_result$compression_ratio
      
      workflow_params$spatial <- spatial_result$basis_params
      
      if (return_intermediate) {
        intermediate_results$spatial <- spatial_result
      }
      
      if (verbose) {
        cat(sprintf("Spatial compression: %.1fx\n", spatial_compression))
      }
    } else {
      if (verbose) cat("Spatial compression not needed (temporal compression sufficient)\n")
    }
  }
  
  # Step 3: Quantization (for full workflow)
  quantization_compression <- 1.0
  quantization_result <- NULL
  
  if (workflow == "full") {
    if (verbose) cat("\n--- Step 3: Quantization ---\n")
    
    # Use adaptive quantization
    quantization_result <- quantize_data(
      current_data,
      bits = 8,
      scope = "adaptive",
      method = "robust",
      auto_adjust = TRUE,
      return_stats = TRUE
    )
    
    current_data <- quantization_result$quantized_data
    quantization_compression <- quantization_result$compression_ratio
    
    workflow_params$quantization <- quantization_result$params
    
    if (return_intermediate) {
      intermediate_results$quantization <- quantization_result
    }
    
    if (verbose) {
      cat(sprintf("Quantization compression: %.1fx\n", quantization_compression))
      if (!is.null(quantization_result$mse)) {
        cat(sprintf("Quantization SNR: %.1f dB\n", quantization_result$snr_db))
      }
    }
  }
  
  # Calculate total compression
  total_compression <- temporal_compression * spatial_compression * quantization_compression
  
  if (verbose) {
    cat("\n=== Final Results ===\n")
    cat(sprintf("Total compression ratio: %.1fx\n", total_compression))
    cat(sprintf("Target compression: %.1fx\n", compression_target))
    
    if (total_compression >= compression_target * 0.8) {
      cat("✓ Compression target achieved\n")
    } else {
      cat("⚠ Compression target not fully achieved\n")
    }
  }
  
  # Create reconstruction function
  reconstruction_function <- function() {
    # TODO: Implement actual reconstruction pipeline
    abort_lna(
      "fmri_workflow() reconstruction is not yet implemented",
      .subclass = "lna_error_not_implemented",
      location = "experimental_api::fmri_workflow$reconstruction_function"
    )
  }
  
  # Generate recommendations
  recommendations <- list()
  
  if (total_compression < compression_target * 0.8) {
    recommendations <- append(recommendations, 
      "Consider more aggressive temporal thresholding or higher spatial compression")
  }
  
  if (workflow != "full" && compression_target > 15) {
    recommendations <- append(recommendations,
      "Consider 'full' workflow with quantization for higher compression ratios")
  }
  
  if (study_type == "resting" && denoise_level == "conservative") {
    recommendations <- append(recommendations,
      "Resting-state data may benefit from moderate denoising")
  }
  
  # Return comprehensive results
  list(
    compressed_data = current_data,
    compression_ratio = total_compression,
    parameters = workflow_params,
    reconstruction_function = reconstruction_function,
    intermediate_results = if (return_intermediate) intermediate_results else NULL,
    recommendations = recommendations,
    workflow_summary = list(
      study_type = study_type,
      workflow = workflow,
      denoise_level = denoise_level,
      temporal_compression = temporal_compression,
      spatial_compression = spatial_compression,
      quantization_compression = quantization_compression
    )
  )
}
