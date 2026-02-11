# Rcpp Sprint: Performance Optimization Plan for neuroarchive

**Date**: 2025-10-27
**Goal**: Identify and implement high-impact C++ optimizations for compression pipeline
**Expected Overall Speedup**: 8-15x for full compression pipelines

---

## Executive Summary

Analysis of `vignettes/compression_methods.Rmd` and core transform implementations reveals significant opportunities for Rcpp acceleration. The package already has C++ implementations for HRBF, Haar wavelets, and block quantization. This document outlines remaining high-value targets.

**Key Finding**: Most bottlenecks are in data movement (copies, transposes) and element-wise operations that don't leverage BLAS optimization.

---

## Current C++ Implementation Status ✓

### Already Optimized (Do Not Re-implement)
- [x] **Spatial HRBF** (`lna_hrbf_helpers.cpp`, `lna_hrbf_rcpp.cpp`) - OpenMP parallelized
- [x] **Haar Octree Wavelet** (`haar_octwave.cpp`) - Uses Morton hashing
- [x] **Quantization Block-wise** (`quant_block_rcpp.cpp`) - Voxel-wise quantization for 4D blocks
- [x] **Morton Hashing** (`morton_hash_rcpp.cpp`) - Spatial indexing
- [x] **Sobel 3D** (`sobel3d.cpp`) - Edge detection for HRBF
- [x] **OMP Utilities** (`lna_omp_rcpp.cpp`) - OpenMP helpers

---

## HIGH PRIORITY Implementations 🔥

### [ ] 1. Global Quantization (3-5x speedup)

**File**: `R/transform_quant.R:415-480` (`.quantize_global()`)
**Target File**: `src/quantize_global_rcpp.cpp`
**Lines of Code**: ~150
**Used By**: Nearly all compression pipelines

#### Current Bottleneck
```r
# Multiple full passes over data
rng <- range(x)           # Pass 1: scan for min/max
m <- mean(x)              # Pass 2: compute mean
s <- sd(x)                # Pass 3: compute sd (requires mean)
q_raw <- round((x - offset) / scale)  # Pass 4: quantize
n_clipped <- sum(q_raw < 0 | q_raw > 2^bits - 1)  # Pass 5: count clips
q[q < 0] <- 0L           # Pass 6: clip low
q[q > qmax] <- qmax      # Pass 7: clip high
```

Each pass over 180M element array (typical fMRI) = 1.4GB scanned 7 times = **9.8GB memory traffic**

#### Proposed C++ Implementation
```cpp
// Single-pass: compute stats + quantize + clip
struct QuantStats {
  double sum, sumsq, min_val, max_val;
  int n_clipped_low, n_clipped_high;
};

// [[Rcpp::export]]
List quantize_global_rcpp(NumericVector x, int bits,
                          String method, bool center) {
  const int n = x.size();
  const int qmax = (1 << bits) - 1;

  // Single pass: compute all stats
  QuantStats stats = compute_stats_simd(x);  // SIMD-optimized

  // Compute scale/offset
  double scale, offset;
  compute_scale_offset(stats, bits, method, center, scale, offset);

  // Single pass: quantize + clip + count
  IntegerVector q(n);
  int n_clipped = quantize_and_clip(x, scale, offset, qmax, q);

  return List::create(
    _["q"] = q,
    _["scale"] = scale,
    _["offset"] = offset,
    _["n_clipped_total"] = n_clipped,
    _["clip_pct"] = 100.0 * n_clipped / n
  );
}
```

#### Performance Target
- **Small** (5K elements): 0.1ms → 0.05ms (2x)
- **Medium** (131K elements): 5ms → 1ms (5x)
- **Large** (3.6M elements): 150ms → 30ms (5x)

#### Testing Checklist
- [ ] Matches R output exactly for `method="range"`, `center=TRUE`
- [ ] Matches R output exactly for `method="range"`, `center=FALSE`
- [ ] Matches R output exactly for `method="sd"`, `center=TRUE`
- [ ] Matches R output exactly for `method="sd"`, `center=FALSE`
- [ ] Handles zero-variance input (all identical values)
- [ ] Handles extreme values (near machine precision)
- [ ] Handles 1-bit to 16-bit quantization
- [ ] Clip percentage matches R implementation
- [ ] Benchmark: 3x minimum speedup on 100M element array

---

### [ ] 2. Delta Encoding (2-4x speedup)

**File**: `R/transform_delta.R:78-105` (`forward_step.delta`)
**Target File**: `src/delta_encode_rcpp.cpp`
**Lines of Code**: ~200
**Used By**: Near-lossless compression pipelines

#### Current Bottleneck
```r
# Step 1: Permute axes (full array copy)
perm <- c(axis, setdiff(seq_along(dims), axis))
xp <- aperm(x, perm)  # 180M elements copied

# Step 2: Reshape (metadata operation, fast)
dim(xp) <- c(axis_len, prod(dims[-axis]))

# Step 3: Extract first values (copy)
first_vals <- xp[1, , drop = FALSE]  # Copy slice

# Step 4: Compute differences (allocate + subtract)
deltas <- xp[-1, , drop = FALSE] - xp[-nrow(xp), , drop = FALSE]
# Allocates TWO temporary matrices + result matrix
```

For 91×109×91×200 array (180M elements):
- `aperm()`: 720MB copy
- `xp[-1, ]`: 360MB temporary
- `xp[-nrow(xp), ]`: 360MB temporary
- `deltas`: 360MB result
**Total**: 1.8GB allocated

#### Proposed C++ Implementation
```cpp
// No copies - compute differences directly with axis-aware indexing
// [[Rcpp::export]]
List delta_encode_rcpp(NumericVector x, IntegerVector dims, int axis) {
  // Compute strides for axis-aware indexing
  std::vector<int> strides = compute_strides(dims);
  const int axis_len = dims[axis - 1];  // R is 1-indexed
  const int non_axis_size = compute_non_axis_size(dims, axis);

  // Allocate outputs
  NumericVector first_vals(non_axis_size);
  NumericVector deltas((axis_len - 1) * non_axis_size);

  // Single pass: extract first values + compute deltas
  #pragma omp parallel for if(non_axis_size > 10000)
  for (int i = 0; i < non_axis_size; ++i) {
    int idx_base = compute_base_index(i, axis, strides);

    // First value
    first_vals[i] = x[idx_base];

    // Deltas along axis
    for (int t = 1; t < axis_len; ++t) {
      int idx_curr = idx_base + t * strides[axis - 1];
      int idx_prev = idx_base + (t - 1) * strides[axis - 1];
      deltas[i + (t - 1) * non_axis_size] = x[idx_curr] - x[idx_prev];
    }
  }

  return List::create(
    _["first_values"] = first_vals,
    _["deltas"] = deltas
  );
}
```

#### Performance Target
- **Small** (5K elements): 0.5ms → 0.2ms (2.5x)
- **Medium** (131K elements): 20ms → 5ms (4x)
- **Large** (3.6M elements): 600ms → 150ms (4x)

#### Testing Checklist
- [ ] Matches R output for axis=1 (rows)
- [ ] Matches R output for axis=2 (columns)
- [ ] Matches R output for axis=3 (slices)
- [ ] Matches R output for axis=4 (time - most common)
- [ ] Matches R output for axis=-1 (last dimension)
- [ ] Handles 1D, 2D, 3D, 4D arrays correctly
- [ ] Handles edge case: 1 element along diff axis
- [ ] Handles edge case: 0 elements along diff axis
- [ ] Memory usage ≤ 50% of R implementation
- [ ] Benchmark: 2x minimum speedup on 100M element 4D array

---

### [ ] 3. Temporal Projection with Thresholding (2-3x speedup)

**File**: `R/temporal_project.R` (various `temporal_project.*` methods)
**Target File**: `src/temporal_project_rcpp.cpp`
**Lines of Code**: ~300
**Used By**: Temporal compression with shrinkage

#### Current Bottleneck
```r
# Projection is fast (BLAS-optimized)
coeff <- crossprod(basis, X)  # (n_basis x n_time) %*% (n_time x n_voxel)

# BUT: Thresholding loops over voxels in R
if (threshold_method != "none") {
  for (v in seq_len(ncol(coeff))) {
    coeff[, v] <- apply_threshold(coeff[, v], threshold_method, lambda)
  }
}
```

For 3.6M voxels with soft thresholding:
- R loop overhead: ~2ms per 1000 voxels = 7 seconds
- BLAS projection: 0.5 seconds
**Ratio**: Thresholding is 14x slower than projection!

#### Proposed C++ Implementation
```cpp
// Fused projection + thresholding
// [[Rcpp::export]]
NumericMatrix temporal_project_threshold_rcpp(
    NumericMatrix basis,     // n_time x n_basis
    NumericMatrix X,         // n_time x n_voxel
    String threshold_method,
    double lambda) {

  const int n_basis = basis.ncol();
  const int n_time = basis.nrow();
  const int n_voxel = X.ncol();

  // Use Eigen for BLAS-backed projection
  Eigen::Map<Eigen::MatrixXd> B(REAL(basis), n_time, n_basis);
  Eigen::Map<Eigen::MatrixXd> X_eigen(REAL(X), n_time, n_voxel);

  // Project: coeff = B' * X
  Eigen::MatrixXd coeff = B.transpose() * X_eigen;

  // Apply thresholding in parallel
  if (threshold_method != "none") {
    #pragma omp parallel for schedule(dynamic, 1000)
    for (int v = 0; v < n_voxel; ++v) {
      apply_threshold_inplace(coeff.col(v), threshold_method, lambda);
    }
  }

  return wrap(coeff);
}

// Threshold implementations
void apply_threshold_inplace(Eigen::VectorXd& x,
                             const String& method,
                             double lambda) {
  if (method == "soft") {
    // Soft thresholding: sign(x) * max(|x| - lambda, 0)
    x = x.unaryExpr([lambda](double v) {
      return std::copysign(std::max(std::abs(v) - lambda, 0.0), v);
    });
  } else if (method == "hard") {
    // Hard thresholding: x * I(|x| > lambda)
    x = x.unaryExpr([lambda](double v) {
      return std::abs(v) > lambda ? v : 0.0;
    });
  } else if (method == "bayes_shrink") {
    // Bayes shrinkage with automatic lambda estimation
    double sigma = estimate_noise_mad(x);
    double threshold = sigma * sigma / lambda;
    x = x.unaryExpr([threshold, sigma](double v) {
      return v / (1.0 + threshold / (v * v + 1e-10));
    });
  }
}
```

#### Performance Target
- **Small** (5K voxels): 5ms → 2ms (2.5x)
- **Medium** (131K voxels): 150ms → 60ms (2.5x)
- **Large** (3.6M voxels): 4000ms → 1500ms (2.7x)

#### Testing Checklist
- [ ] Matches R output for `threshold_method="none"` (projection only)
- [ ] Matches R output for `threshold_method="soft"`, various λ
- [ ] Matches R output for `threshold_method="hard"`, various λ
- [ ] Matches R output for `threshold_method="bayes_shrink"`
- [ ] Handles edge case: λ = 0 (no thresholding)
- [ ] Handles edge case: λ → ∞ (zeros everything)
- [ ] OpenMP parallelization works correctly (no race conditions)
- [ ] Benchmark: 2x minimum speedup on 1M voxel dataset with soft thresholding

---

## MEDIUM PRIORITY Implementations

### [ ] 4. Basis Reconstruction with Subsetting (1.5-2x speedup)

**File**: `R/transform_basis.R:13-82` (`invert_step.basis`)
**Target File**: `src/basis_reconstruct_rcpp.cpp`
**Lines of Code**: ~150
**Used By**: PCA/basis read path

#### Current Bottleneck
```r
# ROI subsetting creates copy
if (!is.null(roi_mask)) {
  vox_idx <- which(as.logical(roi_mask))
  basis <- basis[, vox_idx, drop = FALSE]  # Copy submatrix
}

# Time subsetting creates copy
if (!is.null(time_idx)) {
  coeff <- coeff[time_idx, , drop = FALSE]  # Copy submatrix
}

# Matrix multiply (BLAS-optimized but operates on copies)
dense <- coeff %*% basis
```

For 3.6M voxels with 1M ROI:
- ROI subsetting: Copy 80MB
- Time subsetting: Copy 8MB
- Matrix multiply: 8MB result
**Total**: 96MB allocated for 8MB result

#### Proposed C++ Implementation
```cpp
// Fused subsetting + reconstruction - no intermediate copies
// [[Rcpp::export]]
NumericMatrix basis_reconstruct_subset_rcpp(
    NumericMatrix coeff,     // n_time x n_basis
    NumericMatrix basis,     // n_basis x n_voxel
    Nullable<IntegerVector> time_idx_,
    Nullable<IntegerVector> vox_idx_) {

  // Handle optional subsetting
  std::vector<int> time_idx, vox_idx;
  bool subset_time = parse_index(time_idx_, time_idx, coeff.nrow());
  bool subset_vox = parse_index(vox_idx_, vox_idx, basis.ncol());

  const int n_time_out = subset_time ? time_idx.size() : coeff.nrow();
  const int n_vox_out = subset_vox ? vox_idx.size() : basis.ncol();

  NumericMatrix result(n_time_out, n_vox_out);

  // Multiply with on-the-fly subsetting
  #pragma omp parallel for collapse(2) if(n_time_out * n_vox_out > 100000)
  for (int t = 0; t < n_time_out; ++t) {
    for (int v = 0; v < n_vox_out; ++v) {
      int t_src = subset_time ? time_idx[t] - 1 : t;  // R is 1-indexed
      int v_src = subset_vox ? vox_idx[v] - 1 : v;

      double sum = 0.0;
      for (int k = 0; k < coeff.ncol(); ++k) {
        sum += coeff(t_src, k) * basis(k, v_src);
      }
      result(t, v) = sum;
    }
  }

  return result;
}
```

#### Performance Target
- **Full data**: 50ms → 30ms (1.7x) - modest gain due to BLAS efficiency
- **10% ROI + 50% time**: 25ms → 10ms (2.5x) - larger gain from avoiding copies

#### Testing Checklist
- [ ] Matches R output with no subsetting
- [ ] Matches R output with ROI subsetting only
- [ ] Matches R output with time subsetting only
- [ ] Matches R output with both ROI and time subsetting
- [ ] Handles edge case: empty ROI (returns empty matrix)
- [ ] Handles edge case: single voxel ROI
- [ ] Handles edge case: single timepoint
- [ ] Memory usage ≤ 50% of R implementation for 10% subsets
- [ ] Benchmark: 1.5x minimum speedup on 10% ROI + 50% time subset

---

### [ ] 5. RLE Encoding (3-5x speedup)

**File**: `R/utils_rle.R` (`.encode_rle`, `.decode_rle`)
**Target File**: `src/rle_encode_rcpp.cpp`
**Lines of Code**: ~100
**Used By**: Delta transform with `coding_method="rle"`

#### Current Bottleneck
```r
# R loop for run-length encoding
.encode_rle <- function(x) {
  if (length(x) == 0) return(list(values = numeric(0), lengths = integer(0)))

  runs <- rle(x)  # Base R rle is reasonably fast
  list(values = runs$values, lengths = runs$lengths)
}
```

Base R `rle()` is implemented in C, but:
- Not vectorized for certain data patterns
- Extra overhead from S3 dispatch
- Return structure requires conversion

#### Proposed C++ Implementation
```cpp
// [[Rcpp::export]]
List rle_encode_rcpp(NumericVector x) {
  const int n = x.size();
  if (n == 0) {
    return List::create(
      _["values"] = NumericVector(0),
      _["lengths"] = IntegerVector(0)
    );
  }

  // Pre-allocate worst case (no runs)
  std::vector<double> values;
  std::vector<int> lengths;
  values.reserve(n);
  lengths.reserve(n);

  double current_val = x[0];
  int current_len = 1;

  for (int i = 1; i < n; ++i) {
    if (x[i] == current_val) {
      current_len++;
    } else {
      values.push_back(current_val);
      lengths.push_back(current_len);
      current_val = x[i];
      current_len = 1;
    }
  }

  // Push final run
  values.push_back(current_val);
  lengths.push_back(current_len);

  return List::create(
    _["values"] = wrap(values),
    _["lengths"] = wrap(lengths)
  );
}

// [[Rcpp::export]]
NumericVector rle_decode_rcpp(NumericVector values, IntegerVector lengths) {
  int total_len = 0;
  for (int i = 0; i < lengths.size(); ++i) {
    total_len += lengths[i];
  }

  NumericVector result(total_len);
  int pos = 0;

  for (int i = 0; i < values.size(); ++i) {
    const double val = values[i];
    const int len = lengths[i];
    for (int j = 0; j < len; ++j) {
      result[pos++] = val;
    }
  }

  return result;
}
```

#### Performance Target
- **Low redundancy** (10% compression): 10ms → 5ms (2x)
- **Medium redundancy** (50% compression): 10ms → 2ms (5x)
- **High redundancy** (90% compression): 10ms → 1ms (10x)

#### Testing Checklist
- [ ] Matches R `rle()` output for random data
- [ ] Matches R `rle()` output for highly redundant data (many runs)
- [ ] Matches R `rle()` output for no redundancy (single-element runs)
- [ ] Handles edge case: empty input
- [ ] Handles edge case: single element
- [ ] Decode inverts encode perfectly (roundtrip test)
- [ ] Benchmark: 3x minimum speedup on typical delta-encoded fMRI data

---

## LOWER PRIORITY - Marginal Gains

### [ ] 6. Voxel-wise Quantization Enhancement

**Status**: Already partially implemented in `src/quant_block_rcpp.cpp`
**Current Usage**: Not fully integrated into R code
**Action**: Integration work rather than new implementation

#### Current Situation
- R code uses `matrixStats::rowMeans2()`, `rowSds()` which are already C-optimized
- `quant_block_rcpp.cpp` exists but only called for block processing
- Speed difference: 1.2-1.5x over matrixStats (not dramatic)

#### Integration Checklist
- [ ] Replace `.quantize_voxel_core()` R implementation with `quant_block_rcpp.cpp`
- [ ] Update `forward_step.quant()` to use Rcpp version by default
- [ ] Benchmark to verify speedup justifies code complexity
- [ ] Test: Matches matrixStats implementation exactly

**Priority**: LOW - matrixStats is already quite efficient

---

## Testing Infrastructure

### Benchmark Suite

**File**: `tests/testthat/test-rcpp_benchmarks.R`

```r
test_that("Rcpp quantization benchmarks", {
  skip_on_cran()

  sizes <- list(
    small = c(10, 10, 10, 50),      # 5K voxels
    medium = c(64, 64, 32, 200),    # 131K voxels
    large = c(91, 109, 91, 400)     # 3.6M voxels
  )

  for (size_name in names(sizes)) {
    dims <- sizes[[size_name]]
    x <- rnorm(prod(dims))
    dim(x) <- dims

    # Benchmark R version
    t_r <- system.time({
      r_result <- .quantize_global(as.vector(x), bits = 8,
                                   method = "range", center = TRUE)
    })

    # Benchmark C++ version
    t_cpp <- system.time({
      cpp_result <- quantize_global_rcpp(as.vector(x), bits = 8,
                                        method = "range", center = TRUE)
    })

    # Verify outputs match
    expect_equal(cpp_result$q, r_result$q)
    expect_equal(cpp_result$scale, r_result$scale, tolerance = 1e-10)
    expect_equal(cpp_result$offset, r_result$offset, tolerance = 1e-10)

    # Report speedup
    speedup <- t_r[3] / t_cpp[3]
    message(sprintf("%s: R=%.2fms, C++=%.2fms, speedup=%.1fx",
                   size_name, t_r[3]*1000, t_cpp[3]*1000, speedup))

    # Minimum speedup requirement
    if (size_name == "large") {
      expect_gte(speedup, 3.0,
                info = "Large array should achieve 3x+ speedup")
    }
  }
})
```

### Correctness Tests

**File**: `tests/testthat/test-rcpp_correctness.R`

Each Rcpp implementation should have:
- [ ] Edge case tests (empty, single element, extreme values)
- [ ] Exact match tests vs. R implementation
- [ ] Numerical stability tests (near-zero, near-infinity)
- [ ] Memory tests (no leaks, reasonable allocation)
- [ ] Thread safety tests (if using OpenMP)

---

## Development Guidelines

### Code Style
- Use Rcpp sugar for readable vector operations
- Use RcppEigen for matrix operations requiring BLAS
- Add OpenMP `#pragma omp parallel for` for loops > 10K iterations
- Document performance characteristics in comments

### Error Handling
```cpp
// Check dimensions early
if (dims.size() != 4) {
  stop("Expected 4D array, got %d dimensions", dims.size());
}

// Validate parameters
if (bits < 1 || bits > 16) {
  stop("bits must be between 1 and 16, got %d", bits);
}
```

### Benchmarking
```cpp
// Include timing macros for development
#ifdef BENCHMARK_MODE
  auto start = std::chrono::high_resolution_clock::now();
  // ... computation ...
  auto end = std::chrono::high_resolution_clock::now();
  Rcout << "Elapsed: "
        << std::chrono::duration<double, std::milli>(end - start).count()
        << " ms\n";
#endif
```

---

## Dependencies

### Required
- Rcpp (>= 1.0.0)
- RcppEigen (>= 0.3.3.0)

### Optional
- RcppParallel (for advanced parallelization)
- BH (Boost headers for advanced data structures)

### Build Configuration
```make
# In src/Makevars
PKG_CXXFLAGS = -DNDEBUG -O3 -march=native
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CXXFLAGS)
```

---

## Success Metrics

### Performance Targets
- [ ] Quantization global: 3x minimum on 100M element arrays
- [ ] Delta encoding: 2x minimum on 100M element 4D arrays
- [ ] Temporal projection: 2x minimum on 1M voxel datasets with thresholding
- [ ] Overall pipeline: 8-15x speedup on realistic workflows

### Quality Targets
- [ ] 100% numerical agreement with R implementations (within machine precision)
- [ ] Zero memory leaks (verified with valgrind)
- [ ] Thread-safe OpenMP usage (no race conditions)
- [ ] Passes R CMD check with no warnings

### Documentation Targets
- [ ] Each function has roxygen2 comments
- [ ] Each Rcpp function has usage examples
- [ ] Benchmark results documented in vignette
- [ ] Migration guide for users

---

## Timeline Estimate

**Phase 1** (Week 1): Quantization + Delta
- High-priority, high-impact
- Estimated: 20 hours (implementation + testing)

**Phase 2** (Week 2): Temporal Projection
- More complex due to thresholding variants
- Estimated: 25 hours (implementation + testing + OpenMP)

**Phase 3** (Week 3): Basis Reconstruction + RLE
- Medium priority, easier implementations
- Estimated: 15 hours (implementation + testing)

**Phase 4** (Week 4): Integration + Documentation
- Update vignettes, benchmarks, documentation
- Estimated: 10 hours

**Total**: ~70 hours of development time

---

## Risk Assessment

### Technical Risks
- **Numerical precision**: Floating-point differences between R and C++
  - *Mitigation*: Use same BLAS/LAPACK, extensive tolerance testing

- **Memory management**: Rcpp memory allocation issues
  - *Mitigation*: Use RAII patterns, test with large datasets

- **OpenMP overhead**: Parallelization slower for small data
  - *Mitigation*: Use `if(n > threshold)` guards on `#pragma omp`

### Integration Risks
- **API changes**: Modifying function signatures
  - *Mitigation*: Keep R wrappers, deprecate old functions gradually

- **Platform compatibility**: OpenMP not available on all systems
  - *Mitigation*: Make OpenMP optional, graceful degradation

### Maintenance Risks
- **Code duplication**: C++ and R versions out of sync
  - *Mitigation*: Remove R version after C++ proven stable, comprehensive tests

---

## References

- **Vignette**: `vignettes/compression_methods.Rmd`
- **Transform implementations**: `R/transform_*.R`
- **Existing C++**: `src/*.cpp`
- **Rcpp documentation**: https://www.rcpp.org/
- **RcppEigen**: https://cran.r-project.org/package=RcppEigen
- **OpenMP**: https://www.openmp.org/

---

## Change Log

- **2025-10-27**: Initial analysis and sprint plan created
