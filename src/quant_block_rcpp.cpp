#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {
inline double safe_sd(double sumsq_over_n, double mean) {
  double var = sumsq_over_n - mean * mean;
  if (var < 0.0) var = 0.0;
  return std::sqrt(var);
}

struct SlabRegion {
  int nx;
  int ny;
  int nz;
  int nt;
  int nxyz;
  int sx0;
  int sy0;
  int sz0;
  int snx;
  int sny;
  int snz;
  int vox;
};

SlabRegion parse_slab_region(const NumericVector& x,
                             const IntegerVector& slab_start_xyz_1based,
                             const IntegerVector& slab_dims_xyz) {
  IntegerVector dims = x.attr("dim");
  if (dims.size() != 4) stop("Input x must be a 4D array");
  if (slab_start_xyz_1based.size() != 3 || slab_dims_xyz.size() != 3) {
    stop("slab_start_xyz_1based and slab_dims_xyz must be integer vectors of length 3");
  }

  SlabRegion r;
  r.nx = dims[0];
  r.ny = dims[1];
  r.nz = dims[2];
  r.nt = dims[3];
  r.nxyz = r.nx * r.ny * r.nz;

  r.sx0 = slab_start_xyz_1based[0] - 1;
  r.sy0 = slab_start_xyz_1based[1] - 1;
  r.sz0 = slab_start_xyz_1based[2] - 1;
  r.snx = slab_dims_xyz[0];
  r.sny = slab_dims_xyz[1];
  r.snz = slab_dims_xyz[2];

  if (r.snx < 1 || r.sny < 1 || r.snz < 1) stop("slab dims must be >= 1");
  if (r.sx0 < 0 || r.sy0 < 0 || r.sz0 < 0) stop("slab start indices must be >= 1");
  if (r.sx0 + r.snx > r.nx || r.sy0 + r.sny > r.ny || r.sz0 + r.snz > r.nz) {
    stop("slab region extends beyond x spatial dimensions");
  }

  r.vox = r.snx * r.sny * r.snz;
  return r;
}

void validate_output_buffers(const SlabRegion& r,
                             const IntegerVector& q_out,
                             const NumericVector& scale_out,
                             const NumericVector& offset_out) {
  IntegerVector q_dims = q_out.attr("dim");
  if (q_dims.size() != 4 || q_dims[0] != r.snx || q_dims[1] != r.sny ||
      q_dims[2] != r.snz || q_dims[3] != r.nt) {
    stop("q_out must have dim = c(slab_x, slab_y, slab_z, time)");
  }
  if (q_out.size() != r.vox * r.nt) {
    stop("q_out length does not match expected slab volume * time");
  }

  IntegerVector scale_dims = scale_out.attr("dim");
  if (scale_dims.size() != 3 || scale_dims[0] != r.snx || scale_dims[1] != r.sny ||
      scale_dims[2] != r.snz) {
    stop("scale_out must have dim = c(slab_x, slab_y, slab_z)");
  }
  if (scale_out.size() != r.vox) {
    stop("scale_out length does not match expected slab volume");
  }

  IntegerVector offset_dims = offset_out.attr("dim");
  if (offset_dims.size() != 3 || offset_dims[0] != r.snx || offset_dims[1] != r.sny ||
      offset_dims[2] != r.snz) {
    stop("offset_out must have dim = c(slab_x, slab_y, slab_z)");
  }
  if (offset_out.size() != r.vox) {
    stop("offset_out length does not match expected slab volume");
  }
}

int quantize_slab_into(const NumericVector& x,
                       const SlabRegion& r,
                       int bits,
                       const std::string& method,
                       bool center,
                       int* qptr,
                       double* scptr,
                       double* ofptr) {
  if (bits < 1 || bits > 30) stop("bits must be between 1 and 30");
  const int qmax = (1 << bits) - 1;
  const double* xptr = REAL(x);
  int nclip = 0;

  for (int zz = 0; zz < r.snz; ++zz) {
    for (int yy = 0; yy < r.sny; ++yy) {
      for (int xx = 0; xx < r.snx; ++xx) {
        const int v_slab = xx + r.snx * (yy + r.sny * zz);
        const int x_full = r.sx0 + xx;
        const int y_full = r.sy0 + yy;
        const int z_full = r.sz0 + zz;
        const int v_full = x_full + r.nx * (y_full + r.ny * z_full);

        double sum = 0.0, sumsq = 0.0;
        double minv = R_PosInf, maxv = R_NegInf;
        for (int t = 0; t < r.nt; ++t) {
          double val = xptr[v_full + r.nxyz * t];
          sum += val;
          sumsq += val * val;
          if (val < minv) minv = val;
          if (val > maxv) maxv = val;
        }
        double mean = sum / r.nt;
        double sd = safe_sd(sumsq / r.nt, mean);

        double sc = 0.0, off = 0.0;
        if (center) {
          double max_abs = (method == "sd")
                             ? 3.0 * sd
                             : std::max(std::abs(maxv - mean), std::abs(minv - mean));
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
          for (int t = 0; t < r.nt; ++t) qptr[v_slab + r.vox * t] = 0;
        } else {
          for (int t = 0; t < r.nt; ++t) {
            double val = xptr[v_full + r.nxyz * t];
            int qi = static_cast<int>(std::llrint((val - off) / sc));
            if (qi < 0) {
              qi = 0;
              ++nclip;
            } else if (qi > qmax) {
              qi = qmax;
              ++nclip;
            }
            qptr[v_slab + r.vox * t] = qi;
          }
        }
        scptr[v_slab] = sc;
        ofptr[v_slab] = off;
      }
    }
  }

  return nclip;
}
}  // namespace

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

// [[Rcpp::export]]
Rcpp::List quantize_voxel_slab_rcpp(NumericVector x,
                                    IntegerVector slab_start_xyz_1based,
                                    IntegerVector slab_dims_xyz,
                                    int bits,
                                    std::string method,
                                    bool center) {
  SlabRegion region = parse_slab_region(x, slab_start_xyz_1based, slab_dims_xyz);
  NumericVector scale(region.vox);
  NumericVector offset(region.vox);
  IntegerVector q(region.vox * region.nt);
  q.attr("dim") = IntegerVector::create(region.snx, region.sny, region.snz, region.nt);
  scale.attr("dim") = IntegerVector::create(region.snx, region.sny, region.snz);
  offset.attr("dim") = IntegerVector::create(region.snx, region.sny, region.snz);

  int nclip = quantize_slab_into(
    x,
    region,
    bits,
    method,
    center,
    INTEGER(q),
    REAL(scale),
    REAL(offset)
  );

  return Rcpp::List::create(
    Rcpp::Named("q") = q,
    Rcpp::Named("scale") = scale,
    Rcpp::Named("offset") = offset,
    Rcpp::Named("n_clipped_total") = nclip
  );
}

// [[Rcpp::export]]
int quantize_voxel_slab_fill_rcpp(NumericVector x,
                                  IntegerVector slab_start_xyz_1based,
                                  IntegerVector slab_dims_xyz,
                                  int bits,
                                  std::string method,
                                  bool center,
                                  IntegerVector q_out,
                                  NumericVector scale_out,
                                  NumericVector offset_out) {
  SlabRegion region = parse_slab_region(x, slab_start_xyz_1based, slab_dims_xyz);
  validate_output_buffers(region, q_out, scale_out, offset_out);

  return quantize_slab_into(
    x,
    region,
    bits,
    method,
    center,
    INTEGER(q_out),
    REAL(scale_out),
    REAL(offset_out)
  );
}
