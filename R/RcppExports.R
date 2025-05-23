# This file was generated manually because Rcpp::compileAttributes is unavailable
#' @useDynLib neuroarchive, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @export
label_components_6N_rcpp <- function(flat_mask, dims) {
  .Call(`_neuroarchive_label_components_6N_rcpp`, flat_mask, dims)
}

#' @export
poisson_disk_sample_component_rcpp <- function(component_vox_coords_0based,
                                              radius_vox_sq, component_seed) {
  .Call(`_neuroarchive_poisson_disk_sample_component_rcpp`,
        component_vox_coords_0based, radius_vox_sq, component_seed)
}
