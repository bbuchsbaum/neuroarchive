#include <Rcpp.h>

using namespace Rcpp;

// Declarations from lna_hrbf_helpers.cpp
IntegerVector label_components_6N_rcpp(LogicalVector flat_mask, IntegerVector dims);
IntegerMatrix poisson_disk_sample_component_rcpp(IntegerMatrix component_vox_coords_0based, double radius_vox_sq, int component_seed);

// Register routines
// [[Rcpp::export]]
RcppExport SEXP _neuroarchive_label_components_6N_rcpp(SEXP flat_maskSEXP, SEXP dimsSEXP) {
  BEGIN_RCPP;
  Rcpp::RObject rcpp_result_gen;
  try {
    Rcpp::traits::input_parameter< LogicalVector >::type flat_mask(flat_maskSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(label_components_6N_rcpp(flat_mask, dims));
  } catch(Rcpp::exception &e) {
    Rcpp::forward_exception_to_r(e);
  }
  return rcpp_result_gen;
  END_RCPP;
}

// [[Rcpp::export]]
RcppExport SEXP _neuroarchive_poisson_disk_sample_component_rcpp(SEXP component_vox_coords_0basedSEXP, SEXP radius_vox_sqSEXP, SEXP component_seedSEXP) {
  BEGIN_RCPP;
  Rcpp::RObject rcpp_result_gen;
  try {
    Rcpp::traits::input_parameter< IntegerMatrix >::type component_vox_coords_0based(component_vox_coords_0basedSEXP);
    Rcpp::traits::input_parameter< double >::type radius_vox_sq(radius_vox_sqSEXP);
    Rcpp::traits::input_parameter< int >::type component_seed(component_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(poisson_disk_sample_component_rcpp(component_vox_coords_0based, radius_vox_sq, component_seed));
  } catch(Rcpp::exception &e) {
    Rcpp::forward_exception_to_r(e);
  }
  return rcpp_result_gen;
  END_RCPP;
}

static const R_CallMethodDef CallEntries[] = {
  {"_neuroarchive_label_components_6N_rcpp", (DL_FUNC) &_neuroarchive_label_components_6N_rcpp, 2},
  {"_neuroarchive_poisson_disk_sample_component_rcpp", (DL_FUNC) &_neuroarchive_poisson_disk_sample_component_rcpp, 3},
  {NULL, NULL, 0}
};

RcppExport void R_init_neuroarchive(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
