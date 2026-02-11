#' Experimental Basis Space Functions
#'
#' @description
#' User-facing functions for spatial basis decomposition using HRBF
#' (Hemodynamic Response Basis Functions). These functions provide
#' simplified interfaces for spatial compression and analysis.
#'
#' @name experimental_basis
NULL

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

  stopifnot(is.matrix(data), is.numeric(data))

  # Set default basis parameters
  default_params <- list(
    sigma0 = 6,
    levels = 3L,
    radius_factor = 2.5,
    kernel_type = "gaussian",
    seed = 1L
  )

  # Merge with user parameters
  params <- modifyList(default_params, basis_params)

  # Generate the HRBF basis (atoms x voxels sparse matrix)
  B <- hrbf_generate_basis(params, mask)

  # Project: coefficients = data %*% t(B)  (time x atoms)
  coefficients <- as.matrix(Matrix::tcrossprod(as_dense_mat(data), B))

  result <- list(
    coefficients = coefficients,
    params = params,
    n_atoms = nrow(B),
    n_voxels = ncol(data),
    n_time = nrow(data)
  )
  if (return_basis) {
    result$basis <- B
  }
  class(result) <- "lna_basis_embed_result"
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

  stopifnot(inherits(embedded_result, "lna_basis_embed_result"))

  coeff <- embedded_result$coefficients
  params <- embedded_result$params

  if (!is.null(subset_time)) {
    coeff <- coeff[subset_time, , drop = FALSE]
  }

  # Regenerate basis and reconstruct: data = coeff %*% B
  B <- if (!is.null(embedded_result$basis)) {
    embedded_result$basis
  } else {
    hrbf_generate_basis(params, mask)
  }

  Xhat <- as.matrix(coeff %*% B)

  if (!is.null(subset_voxels)) {
    Xhat <- Xhat[, subset_voxels, drop = FALSE]
  }

  Xhat
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
#' @keywords internal
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
    location = "experimental_basis::basis_space_visualize"
  )
}
