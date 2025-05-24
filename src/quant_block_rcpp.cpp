#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List quantize_voxel_block_rcpp(NumericVector block,
                                     int bits,
                                     std::string method,
                                     bool center) {
  IntegerVector dims = block.attr("dim");
  if (dims.size() != 4) stop("Input block must be a 4D array");
  const int nx = dims[0];
  const int ny = dims[1];
  const int nz = dims[2];
  const int nt = dims[3];
  const int vox = nx * ny * nz;
  const int qmax = (1 << bits) - 1;

  NumericVector scale(vox);
  NumericVector offset(vox);
  IntegerVector q(vox * nt);
  q.attr("dim") = dims;

  const double *xptr = REAL(block);
  double *scptr = REAL(scale);
  double *ofptr = REAL(offset);
  int *qptr = INTEGER(q);

  for (int v = 0; v < vox; ++v) {
    double sum = 0.0, sumsq = 0.0;
    double minv = R_PosInf, maxv = R_NegInf;
    for (int t = 0; t < nt; ++t) {
      double val = xptr[v + vox * t];
      sum += val;
      sumsq += val * val;
      if (val < minv) minv = val;
      if (val > maxv) maxv = val;
    }
    double mean = sum / nt;
    double sd = std::sqrt(sumsq / nt - mean * mean);

    double sc = 0.0, off = 0.0;
    if (center) {
      double max_abs = (method == "sd") ? 3.0 * sd : std::max(std::abs(maxv - mean), std::abs(minv - mean));
      sc = (2.0 * max_abs) / qmax;
      off = mean - max_abs;
    } else {
      double lo, hi;
      if (method == "sd") {
        lo = mean - 3.0 * sd;
        hi = mean + 3.0 * sd;
      } else {
        lo = minv;
        hi = maxv;
      }
      sc = (hi - lo) / qmax;
      off = lo;
    }

    if (sc == 0.0) {
      sc = 1.0;
      for (int t = 0; t < nt; ++t) qptr[v + vox * t] = 0;
    } else {
      for (int t = 0; t < nt; ++t) {
        double val = xptr[v + vox * t];
        int qi = static_cast<int>(std::llrint((val - off) / sc));
        qptr[v + vox * t] = qi;
      }
    }
    scptr[v] = sc;
    ofptr[v] = off;
  }

  int nclip = 0;
  for (int i = 0; i < vox * nt; ++i) {
    if (qptr[i] < 0) { qptr[i] = 0; ++nclip; }
    else if (qptr[i] > qmax) { qptr[i] = qmax; ++nclip; }
  }

  scale.attr("dim") = IntegerVector::create(nx, ny, nz);
  offset.attr("dim") = IntegerVector::create(nx, ny, nz);

  return Rcpp::List::create(
    Rcpp::Named("q") = q,
    Rcpp::Named("scale") = scale,
    Rcpp::Named("offset") = offset,
    Rcpp::Named("n_clipped_total") = nclip
  );
}

