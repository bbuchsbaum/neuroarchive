#' HRBF Core Utilities
#'
#' These helpers expose analytic HRBF basis generation and simple
#' projection/reconstruction outside of a full LNA pipeline.
#'
#' @param params List of HRBF parameters. Same fields as the
#'   `spat.hrbf` transform descriptor.
#' @param mask A `LogicalNeuroVol` mask defining the voxel grid.
#' @param h5_root Optional HDF5 group if centres are stored in a file.
#' @return For `hrbf_generate_basis`, a sparse matrix with one row per
#'   HRBF atom and columns matching mask voxels.  `hrbf_project_matrix`
#'   returns the coefficient matrix and `hrbf_reconstruct_matrix` returns
#'   a dense numeric matrix.
#' @export
hrbf_generate_basis <- function(params, mask, h5_root = NULL) {
  hrbf_basis_from_params(params, mask, h5_root)
}

#' @rdname hrbf_generate_basis
#' @param X Numeric matrix with time points in rows and voxels in columns.
#' @export
hrbf_project_matrix <- function(X, mask, params, h5_root = NULL) {
  B <- hrbf_generate_basis(params, mask, h5_root)
  X_dense <- as_dense_mat(X)
  tcrossprod(X_dense, B)
}

#' @rdname hrbf_generate_basis
#' @param coeff Coefficient matrix with rows corresponding to time points.
#' @export
hrbf_reconstruct_matrix <- function(coeff, mask, params, h5_root = NULL) {
  B <- hrbf_generate_basis(params, mask, h5_root)
  coeff_dense <- as_dense_mat(coeff)
  coeff_dense %*% B
}
