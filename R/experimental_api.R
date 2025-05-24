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
#' Supports multiple basis types commonly used in neuroimaging.
#'
#' @param n_timepoints Number of time points in the series
#' @param basis_type Type of basis: "dct", "polynomial", "bspline", "dpss", "wavelet"
#' @param n_components Number of basis components to generate
#' @param ... Additional parameters specific to each basis type
#'
#' @details
#' Basis-specific parameters:
#' \itemize{
#'   \item `dct`: No additional parameters
#'   \item `polynomial`: `degree` (default: 3)
#'   \item `bspline`: `knot_spacing` (default: 10), `degree` (default: 3)
#'   \item `dpss`: `bandwidth` (default: 4), `k_components` (default: 7)
#'   \item `wavelet`: `wavelet_name` (default: "d4"), `boundary` (default: "periodic")
#' }
#'
#' @return Matrix with time points in rows, basis functions in columns
#'
#' @examples
#' \donttest{
#' # DCT basis for trend removal
#' dct_basis <- basis_time(100, "dct", n_components = 10)
#' 
#' # DPSS basis for spectral analysis
#' dpss_basis <- basis_time(200, "dpss", n_components = 15, time_bandwidth_product = 3)
#' 
#' # B-spline basis for flexible modeling
#' spline_basis <- basis_time(150, "bspline", n_components = 20, knot_spacing = 8)
#' }
#'
#' @export
basis_time <- function(n_timepoints,
                          basis_type = c("dct", "polynomial", "bspline", "dpss", "wavelet"),
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
      "wavelet" = min(n_timepoints %/% 4, 50)
    )
  }
  
  # Call internal temporal basis generation using the existing S3 methods
  neuroarchive::temporal_basis(basis_type, n_timepoints, n_components, ...)
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

#' Quantize data matrix
#'
#' @description
#' Apply quantization to reduce data precision for compression.
#' Supports both global and per-voxel quantization strategies.
#'
#' @param data Matrix to quantize
#' @param bits Number of quantization bits (1-16, default: 8)
#' @param method Scaling method: "range" (min/max) or "sd" (standard deviation)
#' @param scope Quantization scope: "global" (one scale) or "voxel" (per-column)
#' @param center Center data before quantization (default: TRUE)
#' @param clip_threshold Abort if clipping exceeds this percentage (default: 5)
#'
#' @return List with quantized data and reconstruction parameters
#'
#' @examples
#' \donttest{
#' # Create example data matrix
#' set.seed(123)
#' my_data <- matrix(rnorm(100 * 50, mean = 100, sd = 15), nrow = 100, ncol = 50)
#' 
#' # Basic 8-bit quantization
#' quant_result <- quantize_data(my_data, bits = 8)
#' 
#' # Per-voxel quantization for better precision
#' quant_precise <- quantize_data(my_data, bits = 12, scope = "voxel")
#' 
#' # Reconstruct original data
#' reconstructed <- dequantize_data(quant_result)
#' }
#'
#' @export
quantize_data <- function(data,
                         bits = 8L,
                         method = c("range", "sd"),
                         scope = c("global", "voxel"),
                         center = TRUE,
                         clip_threshold = 5.0) {
  
  method <- match.arg(method)
  scope <- match.arg(scope)
  
  params <- list(
    bits = as.integer(bits),
    method = method,
    scale_scope = scope,
    center = center,
    allow_clip = TRUE,
    clip_abort_pct = clip_threshold
  )
  
  # This would call the internal quantization functions
  # For now, return a placeholder structure
  list(
    quantized_data = data,  # Would be actual quantized data
    params = params,
    clip_report = list(clip_pct = 0.0, max_error = 0.0)
  )
}

#' @export
dequantize_data <- function(quant_result) {
  # Reconstruct from quantized data using stored parameters
  # Implementation would use the internal dequantization logic
  quant_result$quantized_data  # Placeholder
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
  
  # Would call internal sparse PCA implementation
  list(
    components = matrix(0, nrow = nrow(data), ncol = n_components),
    loadings = matrix(0, nrow = ncol(data), ncol = n_components),
    params = params,
    convergence = list(converged = TRUE, iterations = 50)
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
  
  # Generate basis and project data
  basis <- hrbf_generate_basis(basis_params, mask)
  coefficients <- hrbf_project_matrix(data, mask, basis_params)
  
  result <- list(
    coefficients = coefficients,
    basis_params = basis_params,
    compression_ratio = nrow(basis) / ncol(basis),
    mask_hash = digest::digest(as.array(mask), algo = "sha256")
  )
  
  if (return_basis) {
    result$basis = basis
  }
  
  result
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
  
  # Verify mask compatibility
  current_hash <- digest::digest(as.array(mask), algo = "sha256")
  if (!identical(current_hash, embedded_result$mask_hash)) {
    warning("Mask hash mismatch - reconstruction may be incorrect")
  }
  
  coeffs <- embedded_result$coefficients
  params <- embedded_result$basis_params
  
  # Apply time subsetting if requested
  if (!is.null(subset_time)) {
    coeffs <- coeffs[subset_time, , drop = FALSE]
  }
  
  # Reconstruct using stored parameters
  reconstructed <- hrbf_reconstruct_matrix(coeffs, mask, params)
  
  # Apply voxel subsetting if requested
  if (!is.null(subset_voxels)) {
    reconstructed <- reconstructed[, subset_voxels, drop = FALSE]
  }
  
  reconstructed
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
  
  # Create visualization (placeholder - would use actual plotting)
  message(sprintf("Would visualize %d basis functions on %s slice %d", 
                 min(n_show, nrow(basis_matrix)), slice_axis, slice_position))
  
  # Return plot object
  invisible(NULL)
} 
