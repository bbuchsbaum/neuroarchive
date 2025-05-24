#' HRBF Core Utilities
#'
#' These helpers expose analytic HRBF basis generation and simple
#' projection/reconstruction outside of a full LNA pipeline.
#'
#' @param params List of HRBF parameters. Same fields as the
#'   `spat.hrbf` transform descriptor.
#' @param mask A `LogicalNeuroVol` mask defining the voxel grid.
#' @return For `hrbf_generate_basis`, a sparse matrix with one row per
#'   HRBF atom and columns matching mask voxels.  `hrbf_project_matrix`
#'   returns the coefficient matrix and `hrbf_reconstruct_matrix` returns
#'   a dense numeric matrix.
#' @export
hrbf_generate_basis <- function(params, mask) {
  # Clean interface for analytic HRBF basis generation
  # Force analytic generation by ensuring centres_stored = FALSE
  params_clean <- params
  params_clean$centres_stored <- FALSE
  
  # Call internal function without HDF5 parameters
  hrbf_basis_from_params(params_clean, mask)
}

#' @rdname hrbf_generate_basis
#' @param X Numeric matrix with time points in rows and voxels in columns.
#' @export
hrbf_project_matrix <- function(X, mask, params) {
  B <- hrbf_generate_basis(params, mask)
  X_dense <- as_dense_mat(X)
  Matrix::tcrossprod(X_dense, B)  # Use Matrix operations to match transform behavior
}

#' @rdname hrbf_generate_basis
#' @param coeff Coefficient matrix with rows corresponding to time points.
#' @export
hrbf_reconstruct_matrix <- function(coeff, mask, params) {
  B <- hrbf_generate_basis(params, mask)
  coeff_dense <- as_dense_mat(coeff)
  coeff_dense %*% B  # Matrix multiplication will preserve Matrix objects
}
