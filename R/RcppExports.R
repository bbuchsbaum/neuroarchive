# This file was generated manually because Rcpp::compileAttributes is unavailable
#' @useDynLib neuroarchive, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @export
label_components_6N_rcpp <- function(flat_mask, dims) {
  .Call(`_neuroarchive_label_components_6N_rcpp`, flat_mask, dims)
}
