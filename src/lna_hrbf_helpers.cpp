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
IntegerMatrix poisson_disk_sample_component_rcpp(IntegerMatrix component_vox_coords_0based, double radius_vox_sq, int component_seed) {
    RNGScope scope; // ensure RNG gets properly set
    Function set_seed("set.seed");
    set_seed(component_seed);

    int nvox = component_vox_coords_0based.nrow();
    if (nvox == 0) {
        return IntegerMatrix(0, 3);
    }

    // Create permutation of indices using R's sample.int
    Function sample_int("sample.int");
    IntegerVector perm = sample_int(nvox, nvox, false);

    std::vector< std::array<int,3> > selected;
    selected.reserve(nvox); // worst case

    for (int idx = 0; idx < nvox; ++idx) {
        int row = perm[idx] - 1; // convert 1-based to 0-based index
        int i = component_vox_coords_0based(row, 0);
        int j = component_vox_coords_0based(row, 1);
        int k = component_vox_coords_0based(row, 2);

        bool accept = true;
        for (const auto& centre : selected) {
            double dx = static_cast<double>(i - centre[0]);
            double dy = static_cast<double>(j - centre[1]);
            double dz = static_cast<double>(k - centre[2]);
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < radius_vox_sq) {
                accept = false;
                break; // early exit
            }
        }
        if (accept) {
            selected.push_back({i,j,k});
        }
    }

    int nsel = selected.size();
    IntegerMatrix out(nsel, 3);
    for (int r = 0; r < nsel; ++r) {
        out(r,0) = selected[r][0];
        out(r,1) = selected[r][1];
        out(r,2) = selected[r][2];
    }
    return out;

}

