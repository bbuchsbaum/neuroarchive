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

