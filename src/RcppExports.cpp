#include <Rcpp.h>

using namespace Rcpp;

// Declarations from lna_hrbf_helpers.cpp
IntegerVector label_components_6N_rcpp(LogicalVector flat_mask, IntegerVector dims);

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

static const R_CallMethodDef CallEntries[] = {
  {"_neuroarchive_label_components_6N_rcpp", (DL_FUNC) &_neuroarchive_label_components_6N_rcpp, 2},
  {NULL, NULL, 0}
};

RcppExport void R_init_neuroarchive(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
