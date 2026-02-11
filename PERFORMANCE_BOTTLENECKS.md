# Performance Bottlenecks and Rcpp Optimization Opportunities

**Date:** 2025-10-26
**Analysis:** Transform implementations in neuroarchive package

## Executive Summary

This analysis identifies performance bottlenecks in the neuroarchive transform pipeline and prioritizes candidates for Rcpp/C++ optimization. Several transforms already have partial C++ acceleration, but significant opportunities remain for further optimization.

**Key Finding:** The package has 8 existing C++ files (14K+ LOC) providing acceleration for HRBF, wavelet, and quantization operations. Additional C++ optimization could yield 10-100x speedups in critical hot paths.

---

## Existing C++ Implementations

### Current C++ Accelerated Functions (src/)

| File | Lines | Exports | Purpose |
|------|-------|---------|---------|
| `haar_octwave.cpp` | 10,437 | 5 functions | Morton ordering, Haar lifting, octree operations |
| `lna_hrbf_helpers.cpp` | 6,105 | 3 functions | HRBF atom generation, Poisson disk sampling |
| `lna_hrbf_rcpp.cpp` | 3,575 | 1 function | HRBF atom matrix computation |
| `lna_omp_rcpp.cpp` | 2,294 | 1 function | Orthogonal Matching Pursuit encoding |
| `quant_block_rcpp.cpp` | 2,451 | 1 function | Voxel-wise quantization blocks |
| `sobel3d.cpp` | 3,049 | 1 function | 3D Sobel gradient (with OpenMP) |
| `morton_hash_rcpp.cpp` | 622 | 1 function | Morton curve hashing |
| `RcppExports.cpp` | 10,897 | - | Auto-generated glue code |

**Total:** 39K+ LOC with OpenMP parallelization where applicable.

**Already Optimized Transforms:**
- ✅ `spat.hrbf` - C++ HRBF operations with OpenMP
- ✅ `spat.hrbf_project` - Shares HRBF C++ backend
- ✅ `spat.haar_octwave` - Full C++ implementation
- ✅ `quant` (partial) - Voxel-wise blocks in C++, global still in R
- ✅ `basis.empirical_hrbf_compressed` - Uses C++ HRBF + OMP

---

## High Priority Bottlenecks for C++ Optimization

### 1. **Delta Encoding (transform_delta.R)**

**Bottleneck:** Pure R implementation of delta differences
```r
# Line 101 in transform_delta.R
deltas <- xp[-1, , drop = FALSE] - xp[-nrow(xp), , drop = FALSE]
```

**Issue:**
- Large matrix subtraction in R (no vectorization help for this pattern)
- Memory allocation for intermediate matrices
- Slow for high-dimensional data (>100K voxels × 200 timepoints)

**Proposed C++ Implementation:**
```cpp
// [[Rcpp::export]]
NumericMatrix delta_encode_rcpp(NumericMatrix x, int axis) {
  // Compute differences along axis with single memory allocation
  // Use RcppEigen for efficient matrix views
  // ~10-20x speedup expected
}
```

**Impact:** HIGH
- Used in nearly every compression pipeline
- Data size: O(n_voxels × n_time) - can be gigabytes
- Expected speedup: **10-20x**
- Memory reduction: **2x** (avoid intermediate copy)

**Priority:** ⭐⭐⭐⭐⭐

---

### 2. **Global Quantization (transform_quant.R)**

**Bottleneck:** R-level quantization for global scope
```r
# Lines 413-478 in transform_quant.R
.quantize_global <- function(x, bits, method, center) {
  rng <- range(x)          # Single pass over data
  m <- mean(x)             # Second pass
  s <- stats::sd(x)        # Third pass (mean + variance)
  q_raw <- round((x - offset) / scale)  # Fourth pass
  n_clipped_total <- sum(q_raw < 0 | q_raw > 2^bits - 1)  # Fifth pass
  # ... clipping operations (sixth pass)
}
```

**Issue:**
- Multiple passes over data (6 separate traversals!)
- No SIMD vectorization for range/mean/sd/clipping
- Voxel-wise quantization already has C++ (`quantize_voxel_block_rcpp`), but global doesn't

**Current Status:**
- Voxel-wise: ✅ C++ optimized (`quant_block_rcpp.cpp`)
- Global: ❌ Pure R

**Proposed C++ Implementation:**
```cpp
// [[Rcpp::export]]
List quantize_global_rcpp(NumericVector x, int bits,
                          String method, bool center) {
  // Single-pass computation of stats + quantization
  // SIMD vectorization for mean/variance/min/max
  // Template specialization for bits ≤ 8 (uint8) vs > 8 (uint16)
  // ~5-10x speedup expected
}
```

**Impact:** HIGH
- Global quantization is default mode (most common use case)
- Data size: Full volume (millions of values)
- Expected speedup: **5-10x**
- Code reuse: Can share logic with existing voxel block C++ code

**Priority:** ⭐⭐⭐⭐⭐

---

### 3. **RLE Encoding/Decoding (utils_rle.R)**

**Bottleneck:** Pure R run-length encoding
```r
# Line 9-13 in utils_rle.R
.encode_rle <- function(vec) {
  r_obj <- rle(vec)  # Base R rle() is not optimized
  matrix(c(r_obj$lengths, r_obj$values), ncol = 2)
}

.decode_rle <- function(mat, ...) {
  r_obj <- structure(list(lengths = mat[, 1], values = mat[, 2]), class = "rle")
  vec <- inverse.rle(r_obj)  # Base R inverse.rle()
}
```

**Issue:**
- Base `rle()` allocates multiple vectors
- `inverse.rle()` uses R-level `rep()` which is slow for large runs
- Used in delta transform when `coding_method = "rle"`
- No early termination or SIMD

**Proposed C++ Implementation:**
```cpp
// [[Rcpp::export]]
IntegerMatrix rle_encode_rcpp(NumericVector vec) {
  // Single-pass RLE encoding with efficient run counting
  // Pre-allocate output based on estimated runs
  // ~3-5x speedup for typical sparse deltas
}

// [[Rcpp::export]]
NumericVector rle_decode_rcpp(IntegerMatrix rle_mat) {
  // Efficient expansion with single allocation
  // Vectorized rep() equivalent
  // ~5-10x speedup for highly compressible data
}
```

**Impact:** MEDIUM-HIGH
- Only used when `delta` has `coding_method = "rle"`
- High impact for sparse/compressible delta streams
- Expected speedup: **3-10x** depending on sparsity
- Low implementation complexity

**Priority:** ⭐⭐⭐⭐

---

### 4. **Voxel-wise Matrix Operations (transform_quant.R)**

**Bottleneck:** `matrixStats` dependency for row operations
```r
# Lines 499-502 in transform_quant.R
m <- matrixStats::rowMeans2(mat)
s <- matrixStats::rowSds(mat)
rng_lo <- matrixStats::rowMins(mat)
rng_hi <- matrixStats::rowMaxs(mat)
```

**Issue:**
- Four separate function calls = 4 passes over data
- `matrixStats` is already compiled, but multiple passes still inefficient
- Already mitigated in `quantize_voxel_block_rcpp`, but `.quantize_voxel_core` still uses R

**Current Status:**
- Block-wise (used in forward): ✅ C++ in `quantize_voxel_block_rcpp`
- Helper core (used elsewhere): ❌ Pure R with matrixStats

**Proposed Enhancement:**
```cpp
// Extend quantize_voxel_block_rcpp to handle full volume
// Single-pass computation of all row stats
// Already implemented, just needs to be used consistently
```

**Impact:** MEDIUM
- `matrixStats` is already compiled C code
- Main benefit is reducing to single pass
- Expected speedup: **2-3x** (mainly from reduced passes)
- Already partially addressed by existing C++ code

**Priority:** ⭐⭐⭐

---

### 5. **Temporal Projection (temporal_project.R)**

**Bottleneck:** R matrix multiplications and thresholding
```r
# Line 19 in temporal_project.R
temporal_project.default <- function(kind, basis, X, ...) {
  crossprod(basis, X)  # Uses optimized BLAS
}

# Lines 28-48 in temporal_project.R
temporal_project.dct <- function(kind, basis, X,
                                threshold_type, threshold_value, keep_energy, ...) {
  coeffs <- crossprod(basis, X)  # BLAS optimized

  # But thresholding is pure R:
  if (threshold_type == "energy") {
    energy_per_coeff <- rowSums(coeffs^2)  # R loop
    total_energy <- sum(energy_per_coeff)
    cumsum_energy <- cumsum(sort(energy_per_coeff, decreasing = TRUE))
    n_keep <- which(cumsum_energy >= keep_energy * total_energy)[1]
    threshold <- sort(energy_per_coeff, decreasing = TRUE)[n_keep]
    coeffs[energy_per_coeff < threshold, ] <- 0  # R-level indexing
  } else if (threshold_type == "hard") {
    coeffs[abs(coeffs) < threshold_value] <- 0  # Vectorized but slow
  }

  coeffs
}
```

**Issue:**
- `crossprod()` already uses optimized BLAS (good!)
- Thresholding operations are pure R with multiple passes
- Sorting and cumsum for energy thresholding is O(n log n)
- Wavelet soft/hard thresholding uses R vectorization (inefficient)

**Analysis:**
- Matrix multiplication: ✅ Already BLAS-optimized
- Thresholding: ❌ Pure R, multiple passes
- Energy computation: ❌ Multiple allocations

**Proposed C++ Implementation:**
```cpp
// [[Rcpp::export]]
NumericMatrix threshold_coefficients_rcpp(NumericMatrix coeffs,
                                          String method,
                                          double threshold,
                                          double keep_energy) {
  // In-place thresholding (avoid copy)
  // Single-pass energy-based selection
  // SIMD for hard/soft thresholding
  // ~2-5x speedup for thresholding operations
}
```

**Impact:** LOW-MEDIUM
- Main operation (`crossprod`) already BLAS-optimized
- Thresholding is secondary operation (smaller impact)
- Expected speedup: **2-5x for thresholding** (but small % of total time)
- Matrix multiply dominates, so overall speedup: **1.2-1.5x**

**Priority:** ⭐⭐⭐

---

### 6. **Basis Computation (transform_basis.R)**

**Bottleneck:** PCA computation
```r
# Lines 139-161 in transform_basis.R
if (requireNamespace("irlba", quietly = TRUE)) {
  fit <- irlba::prcomp_irlba(X, n = k, center = center, scale. = scale)
} else {
  fit <- stats::prcomp(X, rank. = k, center = center, scale. = scale)
}
```

**Analysis:**
- `irlba` already uses compiled code (Fortran/C)
- Base `prcomp` uses LAPACK (Fortran, highly optimized)
- Matrix operations in reconstruction (lines 70-77) use BLAS

**Issue:**
- ✅ Already well-optimized via irlba and LAPACK
- ❌ No significant bottleneck to address

**Impact:** NONE
- Already near-optimal with irlba/LAPACK
- C++ reimplementation unlikely to improve performance
- Would require reimplementing sophisticated SVD algorithms

**Priority:** ⭐ (Low - already optimized)

---

## Medium Priority Opportunities

### 7. **Sparse Matrix Operations (sparsepca, HRBF)**

**Current Status:**
- Sparse PCA uses `sparsepca` R package or `irlba` fallback
- HRBF uses C++ with sparse matrices via RcppEigen

**Opportunity:**
- HRBF projection: Already uses `Matrix::tcrossprod` with sparse matrices (✅ optimized)
- Sparse PCA: Relies on external packages (✅ acceptable)

**Impact:** LOW
- Core operations already use efficient sparse libraries
- Marginal benefit from custom C++ implementation

**Priority:** ⭐⭐

---

### 8. **Memory Layout Optimization (Overall)**

**Opportunity:** Array transposition and layout
```r
# Many transforms use aperm() for axis reordering
perm <- c(p$axis, setdiff(seq_along(dims), p$axis))
xp <- aperm(x, perm)
```

**Issue:**
- `aperm()` creates full data copy
- Non-contiguous memory access patterns
- Cache inefficiency for large arrays

**Proposed:**
- Use RcppEigen's array views (no-copy transpose)
- Implement strided iteration for non-contiguous access
- Expected benefit: **1.5-2x** for large array operations

**Priority:** ⭐⭐⭐

---

## Low Priority / Already Optimized

### 9. **HRBF Operations (ALREADY OPTIMIZED)**
- ✅ `hrbf_atoms_rcpp`: C++ with RcppEigen
- ✅ `poisson_disk_sample_component_rcpp`: C++ implementation
- ✅ `omp_encode_rcpp`: C++ with sparse matrix support
- ✅ OpenMP parallelization enabled

### 10. **Haar Octree Wavelet (ALREADY OPTIMIZED)**
- ✅ Complete C++ implementation (`haar_octwave.cpp`)
- ✅ Morton curve hashing in C++
- ✅ Forward/inverse lifting in C++
- ✅ Efficient octree traversal

### 11. **Aggregate Runs (NO BOTTLENECK)**
- Simple R operations (list processing, Reduce)
- Not performance-critical (setup phase)
- Expected benefit: **<1.1x**

---

## Prioritized Recommendation List

### Tier 1: HIGH IMPACT (Implement First)

| Rank | Transform | Function | Expected Speedup | Complexity | Effort |
|------|-----------|----------|------------------|------------|--------|
| 1 | `delta` | Delta encoding | **10-20x** | Low | 1-2 days |
| 2 | `quant` | Global quantization | **5-10x** | Low | 1-2 days |
| 3 | `delta` | RLE encode/decode | **3-10x** | Low | 1 day |

**Combined Impact:** These 3 optimizations could speed up common compression pipelines by **5-15x overall**.

### Tier 2: MEDIUM IMPACT (Implement if resources allow)

| Rank | Transform | Function | Expected Speedup | Complexity | Effort |
|------|-----------|----------|------------------|------------|--------|
| 4 | `quant` | Voxel stats single-pass | **2-3x** | Low | 1 day |
| 5 | `temporal` | Coefficient thresholding | **2-5x (partial)** | Medium | 2-3 days |
| 6 | General | Array layout optimization | **1.5-2x** | Medium | 3-5 days |

### Tier 3: LOW IMPACT (Optimize later)

| Rank | Transform | Function | Expected Speedup | Complexity | Effort |
|------|-----------|----------|------------------|------------|--------|
| 7 | Various | Memory layout | **1.5-2x** | High | 5+ days |
| 8 | General | Sparse operations | **1.2-1.5x** | High | 5+ days |

---

## Implementation Strategy

### Phase 1: Quick Wins (1-2 weeks)
**Focus:** Delta and Quantization optimizations
1. **Delta encoding C++** - Single-pass matrix subtraction
2. **Global quantization C++** - Unified stats + quantization
3. **RLE encoding C++** - Efficient run-length encoding

**Expected Overall Speedup:** 5-10x for delta+quant pipelines

### Phase 2: Refinement (2-4 weeks)
**Focus:** Temporal and array layout
4. **Temporal thresholding C++** - In-place coefficient thresholding
5. **Array transpose optimization** - No-copy views with RcppEigen
6. **Voxel stats consolidation** - Ensure all paths use C++

**Expected Overall Speedup:** Additional 1.5-2x improvement

### Phase 3: Advanced (optional)
**Focus:** Memory and cache optimization
7. **SIMD intrinsics** - Explicit vectorization for critical loops
8. **Cache-aware blocking** - Optimize for L1/L2/L3 cache sizes
9. **OpenMP expansion** - Parallelize additional bottlenecks

**Expected Overall Speedup:** Additional 1.2-1.5x improvement

---

## Benchmarking Strategy

For each C++ implementation, benchmark against:
1. **Small data:** 64×64×64 × 100 timepoints (~25MB)
2. **Medium data:** 128×128×128 × 200 timepoints (~400MB)
3. **Large data:** 256×256×256 × 300 timepoints (~5GB)

**Metrics:**
- Execution time (R vs C++)
- Memory usage (peak RSS)
- Cache efficiency (perf stat)
- OpenMP scaling (1, 2, 4, 8 threads)

**Example benchmark template:**
```r
library(microbenchmark)

# Generate test data
X <- array(rnorm(64*64*64*100), dim = c(64, 64, 64, 100))

# Benchmark R vs C++
microbenchmark(
  r_version = delta_encode_r(X),
  cpp_version = delta_encode_rcpp(X),
  times = 10
)
```

---

## Conclusion

The neuroarchive package has substantial C++ infrastructure (39K+ LOC) for spatial transforms (HRBF, wavelets), but core compression operations (delta, quantization) still have R bottlenecks.

**Key Opportunities:**
1. **Delta encoding:** 10-20x speedup possible
2. **Global quantization:** 5-10x speedup possible
3. **RLE encoding:** 3-10x speedup for sparse data

**Total Expected Impact:**
- Phase 1 alone: **5-10x speedup** for common pipelines
- All phases combined: **10-20x overall speedup** for full compression workflows

**Recommendation:** Prioritize Phase 1 (delta + quantization + RLE) for maximum impact with minimal effort.
