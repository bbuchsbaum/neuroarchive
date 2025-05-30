#' Placeholder synthesis for Haar octwave coefficients
#'
#' This helper is a stub and simply returns a zero matrix of the
#' appropriate dimensions. It will be replaced by a real implementation
#' of the inverse lifting scheme.
#' @keywords internal
perform_haar_lift_synthesis <- function(coeff_list, mask_3d_array, levels) {
  n_vox <- sum(mask_3d_array)
  if (is.null(coeff_list$root)) return(matrix(0, 0, n_vox))
  n_time <- nrow(coeff_list$root)
  matrix(0, nrow = n_time, ncol = n_vox)
}
