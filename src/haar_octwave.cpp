#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List forward_lift_rcpp(NumericVector data_masked_morton_ordered,
                             LogicalVector mask_flat_morton_ordered,
                             IntegerVector mask_dims,
                             int levels,
                             List scaling_factors_per_level) {
  int nvox = data_masked_morton_ordered.size();
  double sum = 0.0;
  for (int i = 0; i < nvox; ++i) sum += data_masked_morton_ordered[i];
  double root_coeff = nvox > 0 ? (sum / nvox) * std::sqrt((double)nvox) : NA_REAL;

  List detail_list(levels);
  for (int lvl = 0; lvl < levels; ++lvl) {
    detail_list[lvl] = NumericVector(nvox, 0.0);
  }

  return Rcpp::List::create(Rcpp::Named("root_coeff") = root_coeff,
                            Rcpp::Named("detail_coeffs_by_level") = detail_list);
}

