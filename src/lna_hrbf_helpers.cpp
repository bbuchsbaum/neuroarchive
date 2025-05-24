#include <Rcpp.h>
#include <deque>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector label_components_6N_rcpp(LogicalVector flat_mask, IntegerVector dims) {
  const int nx = dims[0];
  const int ny = dims[1];
  const int nz = dims[2];
  const int total = flat_mask.size();
  IntegerVector labels(total, 0);

  auto to_1d = [nx, ny](int i, int j, int k) {
    return i + j * nx + k * nx * ny;
  };

  std::deque<int> queue;
  int current = 0;

  for (int idx = 0; idx < total; ++idx) {
    if (flat_mask[idx] && labels[idx] == 0) {
      current++;
      labels[idx] = current;
      queue.push_back(idx);

      while (!queue.empty()) {
        int p = queue.front();
        queue.pop_front();

        int i = p % nx;
        int tmp = p / nx;
        int j = tmp % ny;
        int k = tmp / ny;

        // six neighbors in 3D grid
        if (i > 0) {
          int n = to_1d(i - 1, j, k);
          if (flat_mask[n] && labels[n] == 0) {
            labels[n] = current;
            queue.push_back(n);
          }
        }
        if (i + 1 < nx) {
          int n = to_1d(i + 1, j, k);
          if (flat_mask[n] && labels[n] == 0) {
            labels[n] = current;
            queue.push_back(n);
          }
        }
        if (j > 0) {
          int n = to_1d(i, j - 1, k);
          if (flat_mask[n] && labels[n] == 0) {
            labels[n] = current;
            queue.push_back(n);
          }
        }
        if (j + 1 < ny) {
          int n = to_1d(i, j + 1, k);
          if (flat_mask[n] && labels[n] == 0) {
            labels[n] = current;
            queue.push_back(n);
          }
        }
        if (k > 0) {
          int n = to_1d(i, j, k - 1);
          if (flat_mask[n] && labels[n] == 0) {
            labels[n] = current;
            queue.push_back(n);
          }
        }
        if (k + 1 < nz) {
          int n = to_1d(i, j, k + 1);
          if (flat_mask[n] && labels[n] == 0) {
            labels[n] = current;
            queue.push_back(n);
          }
        }
      }
    }
  }

  return labels;
}



// [[Rcpp::export]]
IntegerMatrix poisson_disk_sample_component_rcpp(IntegerMatrix component_vox_coords_0based,
                                                 double radius_vox_sq,
                                                 int component_seed) {
    RNGScope scope; // ensure RNG gets properly set
    Function set_seed("set.seed");
    set_seed(component_seed);

    const int nvox = component_vox_coords_0based.nrow();
    if (nvox == 0) {
        return IntegerMatrix(0, 3);
    }

    // Shuffle candidate order for deterministic randomness
    Function sample_int("sample.int");
    IntegerVector perm = sample_int(nvox, nvox, false);

    // Determine bounding box of the component
    int min_i = component_vox_coords_0based(0, 0);
    int max_i = min_i;
    int min_j = component_vox_coords_0based(0, 1);
    int max_j = min_j;
    int min_k = component_vox_coords_0based(0, 2);
    int max_k = min_k;
    for (int r = 1; r < nvox; ++r) {
        int ii = component_vox_coords_0based(r, 0);
        int jj = component_vox_coords_0based(r, 1);
        int kk = component_vox_coords_0based(r, 2);
        if (ii < min_i) min_i = ii;
        if (ii > max_i) max_i = ii;
        if (jj < min_j) min_j = jj;
        if (jj > max_j) max_j = jj;
        if (kk < min_k) min_k = kk;
        if (kk > max_k) max_k = kk;
    }
    const int nx = max_i - min_i + 1;
    const int ny = max_j - min_j + 1;
    const int nz = max_k - min_k + 1;

    auto to_index = [nx, ny](int i, int j, int k) {
        return i + j * nx + k * nx * ny; // column-major
    };

    std::vector<char> available(nx * ny * nz, 0);
    for (int r = 0; r < nvox; ++r) {
        int idx = to_index(component_vox_coords_0based(r, 0) - min_i,
                           component_vox_coords_0based(r, 1) - min_j,
                           component_vox_coords_0based(r, 2) - min_k);
        available[idx] = 1;
    }

    std::vector< std::array<int,3> > selected;
    selected.reserve(nvox);

    std::deque< std::array<int,3> > q;

    for (int idp = 0; idp < nvox; ++idp) {
        int row = perm[idp] - 1; // convert 1-based to 0-based index
        int ci = component_vox_coords_0based(row, 0);
        int cj = component_vox_coords_0based(row, 1);
        int ck = component_vox_coords_0based(row, 2);

        int start_idx = to_index(ci - min_i, cj - min_j, ck - min_k);
        if (!available[start_idx]) continue;

        selected.push_back({ci, cj, ck});
        available[start_idx] = 0;
        q.clear();
        q.push_back({ci, cj, ck});

        while (!q.empty()) {
            auto p = q.front();
            q.pop_front();
            int pi = p[0];
            int pj = p[1];
            int pk = p[2];
            static const int off[6][3] = {
                {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
            };
            for (int nb = 0; nb < 6; ++nb) {
                int ni = pi + off[nb][0];
                int nj = pj + off[nb][1];
                int nk = pk + off[nb][2];
                if (ni < min_i || ni > max_i ||
                    nj < min_j || nj > max_j ||
                    nk < min_k || nk > max_k) continue;
                int nidx = to_index(ni - min_i, nj - min_j, nk - min_k);
                if (!available[nidx]) continue;
                double dx = static_cast<double>(ni - ci);
                double dy = static_cast<double>(nj - cj);
                double dz = static_cast<double>(nk - ck);
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < radius_vox_sq) {
                    available[nidx] = 0;
                    q.push_back({ni, nj, nk});
                }
            }
        }
    }

    const int nsel = selected.size();
    IntegerMatrix out(nsel, 3);
    for (int r = 0; r < nsel; ++r) {
        out(r, 0) = selected[r][0];
        out(r, 1) = selected[r][1];
        out(r, 2) = selected[r][2];
    }
    return out;
}

