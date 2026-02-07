#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {
struct LevelInfo {
  std::vector<double> sqrt_nvalid;
  std::vector<double> sqrt_nvalid_div_8;
  std::vector<int> counts;
  int n_blocks;
  int total;
};

std::vector<LevelInfo> build_level_info(int levels, List scaling_factors_per_level) {
  if (levels < 1) {
    stop("levels must be >= 1");
  }
  if (scaling_factors_per_level.size() < levels) {
    stop("scaling_factors_per_level must have length >= levels");
  }

  std::vector<LevelInfo> infos;
  infos.reserve(levels);
  for (int lvl = 0; lvl < levels; ++lvl) {
    List sc = scaling_factors_per_level[lvl];
    NumericVector sqrt_nvalid_r = sc["sqrt_nvalid"];
    NumericVector sqrt_nvalid_div_8_r = sc["sqrt_nvalid_div_8"];
    int n_blocks = sqrt_nvalid_r.size();
    if (sqrt_nvalid_div_8_r.size() != n_blocks) {
      stop("sqrt_nvalid_div_8 must match sqrt_nvalid length at level %d", lvl + 1);
    }

    std::vector<double> sqrt_nvalid;
    std::vector<double> sqrt_nvalid_div_8;
    sqrt_nvalid.resize(n_blocks);
    sqrt_nvalid_div_8.resize(n_blocks);
    for (int b = 0; b < n_blocks; ++b) {
      double snv = sqrt_nvalid_r[b];
      double snv8 = sqrt_nvalid_div_8_r[b];
      if (!std::isfinite(snv) || !std::isfinite(snv8)) {
        stop("Non-finite scaling factor encountered at level %d", lvl + 1);
      }
      sqrt_nvalid[b] = snv;
      sqrt_nvalid_div_8[b] = snv8;
    }

    std::vector<int> counts;
    counts.resize(n_blocks);
    int total = 0;
    for (int b = 0; b < n_blocks; ++b) {
      double snv = sqrt_nvalid[b];
      int nv = static_cast<int>(std::llround(snv * snv));
      if (nv <= 0) {
        stop("Encountered non-positive block count in scaling factors.");
      }
      counts[b] = nv;
      total += nv;
    }
    // Internal consistency check
    int total_check = 0;
    for (int b = 0; b < n_blocks; ++b) total_check += counts[b];
    if (total_check != total) {
      stop("Internal error computing block counts at level %d", lvl + 1);
    }

    LevelInfo info;
    info.sqrt_nvalid = std::move(sqrt_nvalid);
    info.sqrt_nvalid_div_8 = std::move(sqrt_nvalid_div_8);
    info.counts = std::move(counts);
    info.n_blocks = n_blocks;
    info.total = total;
    infos.push_back(info);
  }

  return infos;
}
}  // namespace

// [[Rcpp::export]]
List forward_lift_rcpp(NumericVector data_masked_morton_ordered,
                       LogicalVector mask_flat_morton_ordered,
                       IntegerVector mask_dims,
                       int levels,
                       List scaling_factors_per_level) {
  if (levels < 1) {
    stop("levels must be >= 1");
  }

  NumericVector current = clone(data_masked_morton_ordered);
  List detail(levels);

  for (int lvl = 0; lvl < levels; ++lvl) {
    List sc = scaling_factors_per_level[lvl];
    NumericVector sqrt_nvalid = sc["sqrt_nvalid"];
    NumericVector sqrt_nvalid_div_8 = sc["sqrt_nvalid_div_8"];

    int n_blocks = sqrt_nvalid.size();
    NumericVector next_data(n_blocks);
    NumericVector dvec(current.size());

    int idx = 0;
    for (int b = 0; b < n_blocks; ++b) {
      double snv = sqrt_nvalid[b];
      int nv = static_cast<int>(std::llround(snv * snv));
      if (nv <= 0) {
        stop("Encountered non-positive block count in scaling factors (packed scalings should exclude empty blocks).");
      }
      if (idx + nv > current.size()) {
        stop("Scaling factors inconsistent with current vector length (block overrun).");
      }

      double sum = 0.0;
      for (int j = 0; j < nv; ++j) {
        sum += current[idx + j];
      }
      double avg = sum / static_cast<double>(nv);

      next_data[b] = avg * snv;

      double sc_detail = sqrt_nvalid_div_8[b];
      for (int j = 0; j < nv; ++j) {
        dvec[idx + j] = (current[idx + j] - avg) * sc_detail;
      }
      idx += nv;
    }

    if (idx != current.size()) {
      stop("Scaling factors inconsistent with current vector length (block underrun).");
    }

    detail[lvl] = dvec;
    current = next_data;
  }

  return List::create(_["root_coeff"] = current,
                      _["detail_coeffs_by_level"] = detail);
}

// [[Rcpp::export]]
List forward_lift_matrix_rcpp(NumericMatrix data_T_x_Nmask_morton_ordered,
                              LogicalVector mask_flat_morton_ordered,
                              IntegerVector mask_dims,
                              int levels,
                              List scaling_factors_per_level) {
  std::vector<LevelInfo> infos = build_level_info(levels, scaling_factors_per_level);

  int n_time = data_T_x_Nmask_morton_ordered.nrow();
  if (n_time < 1) {
    stop("data_T_x_Nmask_morton_ordered must have at least one row");
  }

  int n_mask = infos[0].total;
  if (data_T_x_Nmask_morton_ordered.ncol() != n_mask) {
    stop("data_T_x_Nmask_morton_ordered must have %d columns (masked voxels in Morton order)", n_mask);
  }

  NumericMatrix root_coeff(n_time, infos[levels - 1].n_blocks);
  std::vector<NumericMatrix> detail_mats;
  detail_mats.reserve(levels);
  List detail(levels);
  for (int lvl = 0; lvl < levels; ++lvl) {
    NumericMatrix mat(n_time, infos[lvl].total);
    detail_mats.push_back(mat);
    detail[lvl] = mat;
  }

  const double* data_ptr = data_T_x_Nmask_morton_ordered.begin();
  double* root_ptr = root_coeff.begin();

  for (int tt = 0; tt < n_time; ++tt) {
    std::vector<double> current;
    current.resize(n_mask);
    for (int k = 0; k < n_mask; ++k) {
      current[k] = data_ptr[tt + n_time * k];
    }

    for (int lvl = 0; lvl < levels; ++lvl) {
      const LevelInfo& info = infos[lvl];
      std::vector<double> next;
      next.resize(info.n_blocks);

      double* det_ptr = detail_mats[lvl].begin();

      int idx = 0;
      for (int b = 0; b < info.n_blocks; ++b) {
        int nv = info.counts[b];
        if (idx + nv > static_cast<int>(current.size())) {
          stop("Scaling factors inconsistent with current vector length (block overrun) at level %d.", lvl + 1);
        }

        double sum = 0.0;
        for (int j = 0; j < nv; ++j) {
          sum += current[idx + j];
        }
        double avg = sum / static_cast<double>(nv);

        double snv = info.sqrt_nvalid[b];
        next[b] = avg * snv;

        double sc_detail = info.sqrt_nvalid_div_8[b];
        for (int j = 0; j < nv; ++j) {
          det_ptr[tt + n_time * (idx + j)] = (current[idx + j] - avg) * sc_detail;
        }

        idx += nv;
      }

      if (idx != static_cast<int>(current.size())) {
        stop("Scaling factors inconsistent with current vector length (underrun) at level %d.", lvl + 1);
      }

      current.swap(next);
    }

    if (static_cast<int>(current.size()) != infos[levels - 1].n_blocks) {
      stop("Internal error: root coefficient width mismatch");
    }
    for (int k = 0; k < infos[levels - 1].n_blocks; ++k) {
      root_ptr[tt + n_time * k] = current[k];
    }
  }

  return List::create(_["root_coeff"] = root_coeff,
                      _["detail_coeffs_by_level"] = detail);
}

// [[Rcpp::export]]
NumericVector inverse_lift_rcpp(NumericVector root_coeff,
                               List detail_coeffs_by_level,
                               LogicalVector mask_flat_morton_ordered,
                               IntegerVector mask_dims,
                               int levels,
                               List scaling_factors_per_level) {
  if (levels < 1) {
    stop("levels must be >= 1");
  }
  if (detail_coeffs_by_level.size() < levels) {
    stop("detail_coeffs_by_level must have length >= levels");
  }

  NumericVector current = clone(root_coeff);

  for (int lvl = levels - 1; lvl >= 0; --lvl) {
    List sc = scaling_factors_per_level[lvl];
    NumericVector sqrt_nvalid = sc["sqrt_nvalid"];
    NumericVector sqrt_nvalid_div_8 = sc["sqrt_nvalid_div_8"];

    int n_blocks = sqrt_nvalid.size();
    if (current.size() != n_blocks) {
      stop("Current coefficient vector length does not match number of blocks at level %d.", lvl + 1);
    }

    int total = 0;
    for (int b = 0; b < n_blocks; ++b) {
      double snv = sqrt_nvalid[b];
      int nv = static_cast<int>(std::llround(snv * snv));
      if (nv <= 0) {
        stop("Encountered non-positive block count in scaling factors.");
      }
      total += nv;
    }

    NumericVector det = detail_coeffs_by_level[lvl];
    if (det.size() != total) {
      stop("Detail coefficient length mismatch at level %d (expected %d, got %d).",
           lvl + 1, total, det.size());
    }

    NumericVector next_data(total);
    int idx_out = 0;
    int idx_det = 0;
    for (int b = 0; b < n_blocks; ++b) {
      double snv = sqrt_nvalid[b];
      int nv = static_cast<int>(std::llround(snv * snv));

      double avg = current[b] / snv;
      double sc_detail = sqrt_nvalid_div_8[b];
      for (int j = 0; j < nv; ++j) {
        next_data[idx_out + j] = det[idx_det + j] / sc_detail + avg;
      }
      idx_out += nv;
      idx_det += nv;
    }

    current = next_data;
  }

  return current;
}

// [[Rcpp::export]]
NumericMatrix inverse_lift_matrix_rcpp(NumericMatrix root_coeff,
                                      List detail_coeffs_by_level,
                                      LogicalVector mask_flat_morton_ordered,
                                      IntegerVector mask_dims,
                                      int levels,
                                      List scaling_factors_per_level) {
  std::vector<LevelInfo> infos = build_level_info(levels, scaling_factors_per_level);
  if (detail_coeffs_by_level.size() < levels) {
    stop("detail_coeffs_by_level must have length >= levels");
  }

  int n_time = root_coeff.nrow();
  if (n_time < 1) {
    stop("root_coeff must have at least one row");
  }
  if (root_coeff.ncol() != infos[levels - 1].n_blocks) {
    stop("root_coeff must have %d columns", infos[levels - 1].n_blocks);
  }

  std::vector<NumericMatrix> detail_mats;
  detail_mats.reserve(levels);
  for (int lvl = 0; lvl < levels; ++lvl) {
    NumericMatrix m = detail_coeffs_by_level[lvl];
    if (m.nrow() != n_time) {
      stop("detail_coeffs_by_level[[%d]] must have %d rows", lvl + 1, n_time);
    }
    if (m.ncol() != infos[lvl].total) {
      stop("detail_coeffs_by_level[[%d]] must have %d columns", lvl + 1, infos[lvl].total);
    }
    detail_mats.push_back(m);
  }

  int n_mask = infos[0].total;
  NumericMatrix out(n_time, n_mask);
  const double* root_ptr = root_coeff.begin();
  double* out_ptr = out.begin();

  std::vector<double> current;

  for (int tt = 0; tt < n_time; ++tt) {
    current.assign(infos[levels - 1].n_blocks, 0.0);
    std::vector<double> next;
    for (int k = 0; k < infos[levels - 1].n_blocks; ++k) {
      current[k] = root_ptr[tt + n_time * k];
    }

    for (int lvl = levels - 1; lvl >= 0; --lvl) {
      const LevelInfo& info = infos[lvl];
      next.assign(info.total, 0.0);

      const double* det_ptr = detail_mats[lvl].begin();
      int idx_out = 0;
      int idx_det = 0;
      for (int b = 0; b < info.n_blocks; ++b) {
        int nv = info.counts[b];
        double snv = info.sqrt_nvalid[b];
        double avg = current[b] / snv;
        double sc_detail = info.sqrt_nvalid_div_8[b];
        for (int j = 0; j < nv; ++j) {
          double det = det_ptr[tt + n_time * (idx_det + j)];
          next[idx_out + j] = det / sc_detail + avg;
        }
        idx_out += nv;
        idx_det += nv;
      }

      current.swap(next);
    }

    if (static_cast<int>(current.size()) != n_mask) {
      stop("Internal error: reconstruction width mismatch");
    }
    for (int k = 0; k < n_mask; ++k) {
      out_ptr[tt + n_time * k] = current[k];
    }
  }

  return out;
}
