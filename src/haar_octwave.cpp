#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List forward_lift_rcpp(NumericVector data_masked_morton_ordered,
                       LogicalVector mask_flat_morton_ordered,
                       IntegerVector mask_dims,
                       int levels,
                       List scaling_factors_per_level) {
  int nvox = data_masked_morton_ordered.size();
  NumericVector root_coeff(1); // placeholder root coefficient
  root_coeff[0] = 0.0;
  List detail(levels);
  for (int l = 0; l < levels; ++l) {
    detail[l] = NumericVector(nvox, 0.0);
  }
  return List::create(_["root_coeff"] = root_coeff,
                      _["detail_coeffs_by_level"] = detail);
}

// [[Rcpp::export]]
NumericVector inverse_lift_rcpp(double root_coeff,
                               List detail_coeffs_by_level,
                               LogicalVector mask_flat_morton_ordered,
                               IntegerVector mask_dims,
                               int levels,
                               List scaling_factors_per_level) {
  int nvox = 0;
  for (int i = 0; i < mask_flat_morton_ordered.size(); ++i) {
    if (mask_flat_morton_ordered[i]) ++nvox;
  }
  NumericVector result(nvox, 0.0);
  return result;
}
